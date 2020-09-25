
Get all PDB files:
```
wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/*
```
* there are 167,570 files (33 Gigabytes) (as of Sep 25, 2020)


Convert pdb file to DSSP file, which contains phi and psi for each residue for each chain:
```
./mkdssp -i ~/pdb-files/pdb1a00.ent.gz -o ~/pdb-files/pdb1a00.dssp
```

Compile mkdssp program:

1. download dssp-2.3.0 from https://github.com/cmbi/dssp/releases
2. unzip
3. cd to directory with source
4.
```
./autogen.sh
./configure
make
```

see also https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html