ComputeCoulOnAtom
=================


Welcome to the ComputeCoulOnAtom code!
This code will compute the overall Coulomb energy between one chosen atom and *all* other atoms of the sample. 
A single configuration (pdb) and its topology (psf) are needed. 
The code assumes dielectric constant to be that of water at 298K (i.e. 78.4), if not otherwise specified.


BUILD (Linux)
-------------
g++ -lm compute_coulombEne_on1atom.cpp -o ComputeCoulOnAtom


USAGE
-----
./ComputeCoulOnAtom <pdbfile> <psffile> <atomid> (<dielectric>)


EXAMPLE
-------
As an example a double stranded DNA (5 base pairs, sense strand:
cyt-cyt-cyt-cyt-cyt) has been generated using Avogadro [1] and
solvated + ionized with VMD [2] plugins. To compute the Coulomb energy
of the Phosphate with id = 33 (on a cytosine residue) with all other
atoms of the sample, the following command should be run:

./ComputeCoulOnAtom ionized.pdb ionized.psf 33

Not specifying the last parameter, the dielectric constant is assumed
to be that of water at 298K. The electrostatic energy of this
Phosphate atom is -0.513721 eV.