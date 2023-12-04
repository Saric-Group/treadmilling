Molecular dynamics model for treadmilling filaments

The model is implemented in LAMMPS using the REACTION packages to execute the polymerisation kinetics.

To run it you need to install a modified version of LAMMPS on your laptop. You can do this easily by executing "source ./build_lammps.sh" from the command line within this directory.

To generate a treadmilling simulation run write_simulation.py in ./generation_files. This will initialise a new folder for the simulation, where an input file (in.local) and a configuration file (configuration.txt) for LAMMPS will be created. This code will also copy all the reaction templates (in /Reactions) to the simulation folder. Please refer to the README file in ./generation_files for details on how to execute write_simulation.py.
