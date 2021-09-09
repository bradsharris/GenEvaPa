#!/usr/bin/env python3

import argparse
import math
import os
import subprocess
import re
import fileinput
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-tol', help='Minimum distance between sugar and water (Angstrom)', required=True)
parser.add_argument('-nd', help='Number of Waters to delete (integer)', required =True)
args = parser.parse_args()

tolerance = float(args.tol)
number_delete = int(args.nd)
## pseudo-Code

def run_gromacs(argument):
	myexec1 = 'gmx'
	myarg = [myexec1] + argument
	mycwd = os.getcwd()
	p = subprocess.Popen(myarg,
							executable=myexec1,
							cwd=mycwd,
							# stdin=subprocess.PIPE,
							# stdout=subprocess.PIPE,
							# stderr=subprocess.PIPE
							)

	p.wait()

def run_sbatch(argument):
	myexec2 = 'sbatch'
	myarg2 = [myexec2] + argument
	mycwd = os.getcwd()
	p = subprocess.Popen(myarg2,
							executable=myexec2,
							cwd=mycwd,
							# stdin=subprocess.PIPE,
							# stdout=subprocess.PIPE,
							# stderr=subprocess.PIPE
							)	

	p.wait()

def run_evaporation(argument):
	myexec3 = './GenEvaPa.py'
	myarg3 = [myexec3] + argument
	mycwd = os.getcwd()
	p = subprocess.Popen(myarg3,
							executable=myexec3,
							cwd=mycwd,
							# stdin=subprocess.PIPE,
							# stdout=subprocess.PIPE,
							# stderr=subprocess.PIPE
							)
	p.wait()

def replaceAll(file,search,replace):
	for line in fileinput.input(file, inplace=1):
		if search in line:
			line = line.replace(search,replace)
		sys.stdout.write(line)

def create_new_file(oldfile,newfile):
	cwd=os.getcwd()
	if os.path.isfile(newfile):
		print("file exists")
	else:
		with open(oldfile, 'r') as ofile:
			lines = ofile.readlines()
		with open(newfile, 'w') as nfile:
			nfile.writelines(lines)

def modify_evap(paste_file,cp_file):
	source_dir=os.getcwd()
	with open(cp_file, 'r') as cfile:
		fulltext = cfile.readlines()
	lastline = fulltext[-1]

	with open(paste_file, 'r') as pfile:
		fullertext = pfile.readlines()
	with open(paste_file, 'w')	as tfile:
		tfile.writelines(fullertext[:-1] + [lastline])

def topology_SOL_delete(topologyfile,number_delete):
	with open(topologyfile, 'r') as tfile:
		lines = tfile.readlines()

	for p in range(-1,-10, -1):
		print(lines[p])
		if "SOL" in lines[p]:
			break


	numSOLS = int(lines[p][3:])
	numSOLS = numSOLS - (number_delete)
	print(numSOLS)
	print(p)
	newline = str(lines[p][:3]) + '         ' + str((numSOLS)) + "\n"
	print(newline)
	with open(topologyfile, 'w') as wfile:
		wfile.writelines(lines[:p] + [newline] + lines[p+1:])

i = 1
j = 1
i_max = 9
j_max = 9
while True:
	k = i + 1
	bn_evap = "Evaporated{:d}.gro".format(i)
	bn_em = "em{:d}".format(i)
	bn_nvt = "nvt{:d}".format(i)
	bn_npt = "npt{:d}".format(i)
	bn_md = "md{:d}".format(i)
	bn_srun = "srun_md{:d}.sh".format(i)
	new_npt = "npt{:d}".format(k)
	new_md = "md{:d}".format(k)
	new_srun = "srun_md{:d}.sh".format(k)

	em_tpr = bn_em + ".tpr"
	em_gro = bn_em + ".gro"
	nvt_tpr = bn_nvt + ".tpr"
	nvt_gro = bn_nvt + ".gro"
	npt_tpr = bn_npt + ".tpr"
	npt_gro = bn_npt + ".gro"
	md_xtc = bn_md + ".xtc"
	md_gro = bn_md + ".gro"
	md_edr = "em{:d}.edr".format(i)
	md_tpr = bn_md + ".tpr"
	energy_out = "energy{:d}.xvg".format(i)
	visc_out = "visco{:d}.xvg".format(i)

	step1 = ['grompp','-f','minim.mdp','-c',bn_evap,'-p','FmocHep.top','-o',em_tpr]
	step2 = ['mdrun','-deffnm',bn_em,'-nt','1']
	step3 = ['grompp','-f','nvt.mdp','-c',em_gro,'-p','FmocHep.top','-o',nvt_tpr]
	step4 = ['mdrun','-deffnm',bn_nvt,'-nt','1']
	# step5 = ['grompp','-f','npt.mdp','-c',nvt_gro,'-p','FmocHep.top','-o',npt_tpr]
	# step6 = ['mdrun','-deffnm',bn_npt]
	step7 = ['grompp','-f','md.mdp','-c',nvt_gro,'-p','FmocHep.top','-o',md_tpr]
	step8 = ['mdrun','-deffnm',bn_md,'-nt','1']


	run_gromacs(step1)
	run_gromacs(step2)

	run_gromacs(step3)
	run_gromacs(step4)

	# run_gromacs(step5)
	# run_gromacs(step6)

	run_gromacs(step7)
	run_gromacs(step8)
	

	# create_new_file(bn_srun,new_srun)
	# replaceAll(new_srun,bn_npt,new_npt)
	# replaceAll(new_srun,bn_md,new_md)

	# step7 = [new_srun]
	# run_sbatch(step7)

	j = j + 1

	new_evap = "Evaporated{:d}.gro".format(j)
	step9 = ['-f',md_gro,'-s',md_gro,'-tol',str(tolerance),'-nd',str(number_delete),'-o',new_evap,'-w',"'name OW'",'-ref',"'resid 1-5'",'-sol',"'resname SOL'",'-misc',"'resname Na'"]
	print(j)
	print(new_evap)
	run_evaporation(step9)
	modify_evap(new_evap,md_gro)

	topology_SOL_delete("FmocHep.top",number_delete)
	
	i = i + 1
	if i > i_max:
		break
	if j > j_max:
		break
