import numpy
import os
import sys
import time
import hashlib
import h5py
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
    ('FirstProgenitor',numpy.int64,1),
    ('LastProgenitor',numpy.int64,1),
    ('NextProgenitor',numpy.int64,1),
    ('Descendant',numpy.int64,1),
    ('FirstHaloInFOFgroup',numpy.int64,1),
    ('NextHaloInFOFgroup',numpy.int64,1),
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
def convert():
    ifile = 100
    folder = "/lustre/scratch/astro/cs390/47Mpc/treedata/"
    lastsnap = 75
    alistfile = "/lustre/scratch/astro/cs390/47Mpc/snap.txt"
    f = h5py.File('trees_'+str(ifile)+".hdf5", 'w')
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
    (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,ifile,verbose)
    #TableFlag
    mergertree_grp.attrs['TableFlag'] = numpy.int32(1)
    #NTree
    mergertree_grp.attrs['NTrees'] = numpy.int32(nTrees)
    #NHalo
    mergertree_grp.attrs['NHalos'] numpy.int32(nHalos)
    #NHalosInTree
    nhalosintree_data = mergertree_grp.create_dataset('NHalosInTree', data=nTreeHalos.astype(numpy.int32))
    #Halo
    nhalosintree_data = mergertree_grp.create_dataset('Halo', data=output_Halos)

def main():
    convert()
if __name__ == "__main__":
    main()    
