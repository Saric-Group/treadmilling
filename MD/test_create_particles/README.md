This is a temporary directory for development purposes

The aim here is to implement treadmilling dynamics in MD (within LAMMPS) using the bond/react fix but creating particles as we polymerise and nucleate, and deleting particles as we depolymerise and dissolve, instead of explicitly considering the free monomers via GCMC like we've done so far. This should allow to play with the growth rates more freely.
