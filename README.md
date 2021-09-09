# GenEvaPa
Python wrapper to model evaporation process using molecular dynamics simulations

# Requirements
* Python3 or greater
* MDAnalysis 1.1.1 or greater

# Installation
Download and place GenEvaPa.py in a directory included in your $PATH variable.

# Arguments
Execute `GenEvaPa.py -h` for a brief description of arguments.

Required Arguments:
* `-s filename` Structure input file (e.g. input.gro) 
* `-f filename` Trajectory input file (e.g. traj.xtc) -can also be the same as structure file
* `-tol ##` Minimum distance between reference and solvent (Angstrom)
* `-nd ##` Number of solvent molecules to delete
* `-o filename` Output structure file (e.g. out.gro)
* `-w 'MDAnalysis selection'` Selection group for solvent to be removed (e.g. 'name OW')
* `-ref 'MDAnalysis selection'` Selection group for reference for solvation shells (e.g. 'resid 1-5')

Optional Arguments:
* `-sol 'MDAnalysis selection'` Selection group for solvent if different from w (e.g. 'resname SOL')
* `-misc 'MDAnalysis selection'` Selection group for other groups in system (e.g. 'resname Na or resname Cl'

# Instructions
First generate molecular dynamics simulation output structure file. Then run `GenEvaPa.py` to generate new molecular dynamics structure file. Finally this structure file can be used as input for further molecular dynamics simulation.

# Additional Information
GenEvaPa does not automatically update topology due to differences between molecular dynamics suites. An example of automating this process can be seen in the included Example.
GenEvaPa defaults to orthorhombic simulation boxes. This can be updated at the bottom of the script, or handled post deletion as shown in the included Example.
