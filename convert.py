import numpy
import os
import sys
import time
import hashlib
import h5py
import numpy.lib.recfunctions as rfn

def hashfile(afile, hasher, blocksize=65536):
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    return hasher.digest()

struct_snapshots = numpy.dtype([
    ('snapnum',numpy.int32,1),
    ('a', numpy.float32, 1),
    ('z', numpy.float32, 1),
    ('t', numpy.float32, 1),
    ('t0', numpy.float32, 1)
])
struct_lgalinput = numpy.dtype([
    ('Descendant',numpy.int32,1),
    ('FirstProgenitor',numpy.int32,1),
    ('NextProgenitor',numpy.int32,1),
    ('FirstHaloInFOFgroup',numpy.int32,1),
    ('NextHaloInFOFgroup',numpy.int32,1),
    ('Len',numpy.int32,1),
    ('M_Mean200',numpy.float32,1),
    ('M_Crit200',numpy.float32,1),
    ('M_TopHat',numpy.float32,1),
    ('Pos',numpy.float32,3),
    ('Vel',numpy.float32,3),
    ('VelDisp',numpy.float32,1),
    ('Vmax',numpy.float32,1),
    ('Spin',numpy.float32,3),
    ('MostBoundID',numpy.int64,1),
    ('SnapNum',numpy.int32,1),
    ('FileNr',numpy.int32,1),
    ('SubhaloIndex',numpy.int32,1),
    ('SubHalfMass',numpy.int32,1)
])
struct_lgaldbidsinput = numpy.dtype([
    ('HaloID',numpy.int64,1),
    ('FileTreeNr',numpy.int64,1),
    ('FirstProgenitorID',numpy.int64,1),
    ('LastProgenitorID',numpy.int64,1),
    ('NextProgenitorID',numpy.int64,1),
    ('DescendantID',numpy.int64,1),
    ('FirstHaloInFOFgroupID',numpy.int64,1),
    ('NextHaloInFOFgroupID',numpy.int64,1),
    ('Redshift',numpy.float64,1),
    ('PeanoKey',numpy.int32,1),
    ('dummy',numpy.int32,1)
])
def read_lgal_input_fulltrees_withids(folder,lastsnap,file,verbose):
    firstfile = file
    lastfile = file 
    nTrees = 0
    nHalos = 0
    nTreeHalos = numpy.array([],dtype=numpy.int32)
    output_Halos = numpy.array([],dtype=struct_lgalinput)
    output_HaloIDs = numpy.array([],dtype=struct_lgaldbidsinput)
    ifile = file
    filename = folder+'/trees_'+"%03d"%(lastsnap)+'.'+"%d"%(ifile)
    f = open(filename,"rb")
    this_nTrees = numpy.fromfile(f,numpy.int32,1)[0]
    nTrees += this_nTrees
    this_nHalos = numpy.fromfile(f,numpy.int32,1)[0]
    nHalos += this_nHalos
    if(verbose):
        print "File ", ifile," nHalos = ",this_nHalos
    nTreeHalos = numpy.fromfile(f,numpy.int32,this_nTrees)
    output_Halos = numpy.fromfile(f,struct_lgalinput,this_nHalos)
    f.close()
    filename = folder+'/tree_dbids_'+"%03d"%(lastsnap)+'.'+"%d"%(ifile)
    f = open(filename,"rb")
    output_HaloIDs = numpy.fromfile(f,struct_lgaldbidsinput,this_nHalos)
    f.close()
    return (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs)

def load_snapshot(alistfile):
    a = numpy.loadtxt(alistfile)
    nsnaps = len(a)
    print nsnaps
    print a
    return (nsnaps,a)
def convert(ifile):
    folder = "/lustre/scratch/astro/cs390/LGalaxies_Hen15_PublicRelease/MergerTrees/MR/treedata/"
    lastsnap = 63
    alistfile = "/lustre/scratch/astro/cs390/LGalaxies_Hen15_PublicRelease/input/zlists/zlist_MR.txt"
    f = h5py.File(folder+'/trees_'+str(ifile)+".hdf5", 'w')
    # Version
    f.attrs.create('Version', 0, dtype=numpy.int32)
    # Subversion
    f.attrs.create('Subversion', 1, dtype=numpy.int32)
    # Title
    f.attrs.create('Title', "The Mighty Peter")
    # Description
    f.attrs.create('Description', "This is for testing")
    # BoxsizeMpc -- I'm not convinced that we should use Mpc instead Mpc/h (It's quite difficult to remember)
    # so I will use Mpc/h to avoid the errors from myself
    f.attrs.create('BoxsizeMpc_h', 62.5, dtype=numpy.float32)
    # OmegaBaryon
    f.attrs.create('OmegaBaryon', 0.044, dtype=numpy.float32)
    # OmegaCDM
    f.attrs.create('OmegaCDM', 0.27-0.044, dtype=numpy.float32)
    # H100
    f.attrs.create('H100', 0.704, dtype=numpy.float32)
    # Sigma8
    f.attrs.create('Sigma8', 0.807, dtype=numpy.float32)
    
    #Group -- Snapshot
    snapshot_grp = f.create_group("Snapshots")
    (nsnaps,snapshot_data) = load_snapshot(alistfile)
    #NSnap
    print numpy.int32(nsnaps)
    snapshot_grp.attrs['NSnap'] = numpy.int32(nsnaps)
    #Snap
    snapshot_snap = snapshot_grp.create_dataset('Snap', data=snapshot_data)

    #Group -- MergerTrees
    mergertree_grp = f.create_group("MergerTrees")
    verbose = 1
    print "Reading tree",ifile
    (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,ifile,verbose)
    print "Done reading tree",ifile
    #TableFlag
    mergertree_grp.attrs['TableFlag'] = numpy.int32(1)
    #NTree
    mergertree_grp.attrs['NTrees'] = numpy.int32(nTrees)
    #NHalo
    mergertree_grp.attrs['NHalos'] = numpy.int32(nHalos)
    #NHalosInTree

    nhalosintree_data = mergertree_grp.create_dataset('NHalosInTree', data=nTreeHalos.astype(numpy.int32))
    #Halo
    print "Merging arrays"
    #halo = rfn.merge_arrays((output_Halos,output_HaloIDs), flatten = True, usemask = False)
    halo = join_struct_arrays((output_Halos,output_HaloIDs))
    print "Done merging arrays"
    halo = rfn.drop_fields(halo,['dummy','PeanoKey'])
    print "Outputting merger trees"
    nhalosintree_data = mergertree_grp.create_dataset('Halo', data=halo)
    print "Done"
def join_struct_arrays(arrays):
    sizes = numpy.array([a.itemsize for a in arrays])
    offsets = numpy.r_[0, sizes.cumsum()]
    n = len(arrays[0])
    joint = numpy.empty((n, offsets[-1]), dtype=numpy.uint8)
    for a, size, offset in zip(arrays, sizes, offsets):
        joint[:,offset:offset+size] = a.view(numpy.uint8).reshape(n,size)
    dtype = sum((a.dtype.descr for a in arrays), [])
    return joint.ravel().view(dtype)
def main():
    for i in range(512):
        convert(i)
if __name__ == "__main__":
    main()    
