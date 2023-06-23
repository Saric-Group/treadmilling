09 Nov 2022, Ivan Palaia.

The new fix_bond_react.cpp contains:
- new custom function rxndiffIvan that takes as input a variable and a fragment index. It computes the difference between the variable for the first atom of the fragment and for the second atom of the fragment: var(1) - var(2).
- correction for bug that occurred when more reactions happen at once and one of them did not have angles or dihedral defined in the final conigutation.

The new fix_gcmc :
- works in 2d as well, for insertion and deletion of individual atoms.

----------------------------------------------------------------------

20 Feb 2023, Christian Vanhille Campos. 

A new flag vector for the reactions is introduced: modify_create_nucrand

This vector has size nreacts (number of reactions implemented in the fix) and stores and integer that defines whether the newly nucleated filament dimer takes a random position or not
  -1 = no random position, the new dimer is just a copy of the nucleator dimer following the map in the reaction files
   1 = random position. Instead of following the map in the reaction files, if this option is implemented the position of the inserted particles is defined as follows:

	The first template particle position is redefined to take a random value within the box (uniform sampling using the local number generator of the particular reaction).
	The rest of the template particles positions are redefined from this new random position using the original displacements of the template nucleator.
	The new positions are now fit to this new template using the map defined in the reaction files.

Note that:

	1) A nucleator molecule is still needed in the system, but its position will not affect the nucleation
	2) The orientation of the nucleator is still used to define the orientation of the random template

----------------------------------------------------------------------

05 Apr 2023, Christian Vanhille Campos.

THERE HAVE BEEN NO NEW MODIFICATIONS (according to my notes) SINCE 20/02/2023. HOWEVER, THE IMPLEMENTATION OF CREATION TIME STORAGE (for the hydrolysis part of treadmilling, when computing poff) HAS CHANGED. NOTE THAT THIS PERTAINS TO THE INPUT SCRIPT AND NOT THE LAMMPS FILES SO NOTHING ELSE WAS CHANGED HERE. The modifications to the input script were implemented to solve bugs when nucleating new filaments, as the old implementation would not store the nucleation time for tail monomers of the new dimers. AS IT IS NOW IT WORKS! To see the input script please refer to the examples in MD/examples and to the simulation-generating Python script in MD/generation_files

-----------------------------------------------------------------------

23 Jun 2023, Ivan Palaia.

I am merging our version of bond/react, with the one from the last distributed version of LAMMPS (15 Jun 2023).