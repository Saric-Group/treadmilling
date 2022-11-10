Ivan Palaia, 9 Nov 2022.

The new fix_bond_react.cpp contains:
- new custom function rxndiffIvan that takes as input a variable and a fragment index. It computes the difference between the variable for the first atom of the fragment and for the second atom of the fragment: var(1) - var(2).
- correction for bug that occurred when more reactions happen at once and one of them did not have angles or dihedral defined in the final conigutation.

The new fix_gcmc :
- works in 2d as well, for insertion and deletion of individual atoms.