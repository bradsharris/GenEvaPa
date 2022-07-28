# Example
10% per deletion step automation sample workflow and scripts

# Requirements 
* Python3 or greater
* MDAnalysis 1.1.1 or greater
* Gromacs 5.1 or Greater

# How to Use
Example commands in linux terminal within example folder without addition of scripts to $PATH

To Generate Evaporated1.gro using GenEvaPa

`./GenEvaPa.py -f md0.gro -s md0.gro -tol 10 -nd 16262 -w 'name OW' -ref 'resid 1-5' -sol 'resname SOL' -misc 'resname Na' -o Evaporated1.gro`

Manually running GenEvaPa like this comes with 2 drawbacks
1) Topology file needs manually updated to remove the deleted solvents (Renumber SOL in topology file to the correct number)
2) Output structure file (Evaporated1.gro) has the wrong box dimensions (Replace the last line with the box dimensions in md0.gro)

These drawbacks can be corrected manually, or can be handled by the `prepAutomation.py` script

`./prepAutomation.py -tol 10 -nd 16262`

After Evaporated1.gro is generated with the correct box dimensions and the topology properly updated the deletion process can be automated using `AutoEvap.py`.

* Using AutoEvap.py for other systems requires modification of the input groups in line 169 and the topology to be deleted from in line 175

An example script for executing `AutoEvap.py` using Slurm on an HPC cluster is provided in run_evap.sh. The script is executed using `sbatch run_evap.sh`
