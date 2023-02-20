README FILE EXPLAINING MODIFICATIONS IMPLEMENTED ON 20/03/2022 TO:
  fix_bond_react.cpp
  fix_bond_react.h

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

Chris - 20/02/2023