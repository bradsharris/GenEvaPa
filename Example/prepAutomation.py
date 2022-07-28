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
		# print(lines[p])
		if "SOL" in lines[p]:
			break


	numSOLS = int(lines[p][3:])
	print("Previous step solvent #:",numSOLS)
	numSOLS = numSOLS - (number_delete)
	print("Current step solvent #:",numSOLS)
	newline = str(lines[p][:3]) + '         ' + str((numSOLS)) + "\n"
	# print(newline)
	with open(topologyfile, 'w') as wfile:
		wfile.writelines(lines[:p] + [newline] + lines[p+1:])

i = 0
j = 0
i_max = 1
j_max = 1
while i <= i_max:
	bn_evap = "Evaporated{:d}.gro".format(i)
	bn_md = "md{:d}".format(i)
	md_xtc = bn_md + ".xtc"
	md_gro = bn_md + ".gro"
	md_edr = "em{:d}.edr".format(i)
	md_tpr = bn_md + ".tpr"
	cur_top = "FmocHep_step{:d}.top".format(i)
	j = j + 1
	new_evap = "Evaporated{:d}.gro".format(j)
	step9 = ['-f',md_gro,'-s',md_gro,'-tol',str(tolerance),'-nd',str(number_delete),'-o',new_evap,'-w','name OW','-ref','resid 1-5','-sol','resname SOL','-misc','resname Na']
	print("preparing step #",j)
	print("new filename is", new_evap)
	run_evaporation(step9)
	modify_evap(new_evap,md_gro)
	new_top = "FmocHep_step{:d}.top".format(j)
	create_new_file(cur_top,new_top)
	topology_SOL_delete(new_top,number_delete)
	quit()
	i = i + 1
	if i > i_max:
		break
	if j > j_max:
		break
