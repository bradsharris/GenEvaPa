# Examples
There are two sample systems included, sample commands below are constructed to be 10% per deletion step and both folders have sample automation scripts included. 

* Sample system 1: `Toy/` contains 4 sugar molecules, 3849 water molecules, and 8 sodium ions to balance charges in a ~5 nm^3 box.
    This system is useful for getting a feel for how the script and automation processes function without requiring long simulation times. 
    
* Sample system 2: `LargeSystem/` contains a starting configuration from the manuscript, with 27 sugar molecules, 162623 water molecules, and 54            neutralizing sodium ions in a roughly ~17 nm^3 box.

# Requirements 
* Python3 or greater
* MDAnalysis 1.1.1 or greater
* Gromacs 5.1 or Greater

# How to Use

There are two primary ways to use GenEvaPa.py, manually adjusting a single topology file by running the script once directly, and in an automated fashion to simulate an evaporation process.

Below are example commands to be ran in a linux terminal within each sample folder without addition of scripts to $PATH

# Manually adjusting a single topology file

To Generate a modified topology file `Evaporated1.gro` in the `LargeSystem/` example using GenEvaPa run the script with corresponding flag options outlined previously or with the -h flag. If the script has been placed into $PATH or has been made executable by running a command such as `chmod a+x GenEvaPa.py`. The script can then be executed as follows:

`./GenEvaPa.py -f md0.gro -s md0.gro -tol 10 -nd 16262 -w 'name OW' -ref 'resid 1-5' -sol 'resname SOL' -misc 'resname Na' -o Evaporated1.gro`

If the script has not been made executable then the script can be run using:

`python3 GenEvaPa.py -f md0.gro -s md0.gro -tol 10 -nd 16262 -w 'name OW' -ref 'resid 1-5' -sol 'resname SOL' -misc 'resname Na' -o Evaporated1.gro`

Similarly, in the `Toy/` system a command such as `./GenEvaPa.py -f md0.gro -s md0.gro -tol 10 -nd 385 -w 'name OW' -ref 'resid 1-5' -sol 'resname SOL' -misc 'resname Na' -o Evaporated1.gro`can be run, with the number of solvent to be deleted set lower to account for the smaller system. 

Manually running GenEvaPa like this comes with 2 drawbacks
1) Topology file needs manually updated to remove the deleted solvents (Renumber SOL in topology file to the correct number)
2) Output structure file (Evaporated1.gro) has the wrong box dimensions (Replace the last line with the box dimensions in md0.gro)

These drawbacks can be corrected manually, or can be scripted such as in the `prepAutomation.py` script in the automated workflow.

# Automated workflow 
To create the first step in the evaporation loop `Evaporated1.gro` a sample script `prepAutomation.py` is included in both folders. This script runs the command from the above section on manual adjustment. The script as written maintains variable adjustment on the tolerance and number of solvent to delete. The script as written assumes an input structure file `md0.gro` and the corresponding flags for the sugar, water, and sodium ions. An example of running this for each sample system:

* LargeSystem: `./prepAutomation.py -tol 10 -nd 16262`
* Toy:  `python3 prepAutomation.py -tol 10 -nd 385`

After Evaporated1.gro is generated with correct box dimensions and matching topology file for the given step (`FmocHep_step1.top`) the deletion process can be automated using `AutoEvap.py`.

* Using AutoEvap.py for other systems requires modification of the input groups in line 168 and the topology to be deleted from in line 133 & 166 (Lines 118, 115, and 123 for `prepAutomation.py` respectively)

To run the included automation locally on a device the following commands can be run:
* LargeSystem: `./AutoEvap.py -tol 10 -nd 16262`
* Toy:  `./AutoEvap.py -tol 10 -nd 385`

An example bash script for executing `AutoEvap.py` using Slurm on an HPC cluster is provided in run_evap.sh. The script is executed using `sbatch run_evap.sh`

If the scripts end due to running out of available solvent to delete they can be restarted from the last step by adjusting the i and j parameters in lines 104 and 105 respectively of both `prepAutomation.py` and `AutoEvap.py` and reducing the tolerance as low as 0 angstrom and/or changing the number of removed solvent molecules in the execution commands. 
