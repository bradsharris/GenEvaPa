#!/usr/bin/env python3
import argparse
import warnings
import math
import MDAnalysis as mda
import random
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import AnalysisFromFunction

class my_bin():
	def __init__(self):
		self.sugar_atoms = []
		self.OW_atoms = []


class my_system():
	def __init__(self,u_in,min_dist):
		self.u = u_in
		self.n = []
		self.unit_cell = [0,0,0]
		self.min_dist = min_dist
		if min_dist < 1.0:
			bin_width_local = 1.0
		else:
			bin_width_local = min_dist
		for j in range(3):
			self.unit_cell[j] = self.u.trajectory[0]._unitcell[j]
			self.n.append(math.ceil(self.unit_cell[j]/bin_width_local))
		self.bins = [ [ [ my_bin() for i in range(self.n[2])] for j in range(self.n[1]) ] for k in range(self.n[0]) ]
		self.init_bins()
		self.deletable_resids = []
		self.deletable_atoms = []
		self.non_deleted_resids = []
		self.non_deleted_waters = []
		self.remaining_waters = []
		self.final_non_deleted_waters = []
		self.final_system_sel = []
		self.all_sugar_shells = []
	def init_bins(self):
		for atom in water_sel:
			ind = []
			for j in range(3):
				if atom.position[j] > self.unit_cell[j]:
					ind.append(math.ceil(atom.position[j]/self.unit_cell[j] * self.n[j])-self.n[j])
				else: 
					ind.append(math.ceil(atom.position[j]/self.unit_cell[j] * self.n[j])-1)
			self.bins[ind[0]][ind[1]][ind[2]].OW_atoms.append(atom)

		for atom in ref_sel:
			ind = []
			for j in range(3):
				if atom.position[j] > self.unit_cell[j]:
					ind.append(math.ceil(atom.position[j]/self.unit_cell[j] * self.n[j])-self.n[j])
				else: 
					ind.append(math.ceil(atom.position[j]/self.unit_cell[j] * self.n[j])-1)
			self.bins[ind[0]][ind[1]][ind[2]].sugar_atoms.append(atom)

	def find_deletable_waters(self):
		for i in range(self.n[0]):
			for j in range(self.n[1]):
				for k in range(self.n[2]):
					for OW_atom in self.bins[i][j][k].OW_atoms:
						if self.is_deletable_OW(OW_atom,i,j,k):
							self.deletable_resids.append(OW_atom.resid)
						else:
							self.non_deleted_resids.append(OW_atom.resid)

	# i is the index, j is the cartesian index
	def fix_pbc(self,i,j):
		image = 0
		if i < 0:
			i += self.n[j]
			image -= 1
		elif i >= self.n[j]:
			i -= self.n[j]
			image += 1
		return i,image

	def find_nearby_sugars(self,i,j,k):
		nearby_sugars = []
		nearby_sugar_images = []
		for i2 in range(i-1,i+2):
			image = [0,0,0]
			[i2p, image[0]] = self.fix_pbc(i2,0)
			for j2 in range(j-1,j+2):
				[j2p, image[1]] = self.fix_pbc(j2,1)
				for k2 in range(k-1,k+2):
					[k2p, image[2]] = self.fix_pbc(k2,2)
					for sugar in self.bins[i2p][j2p][k2p].sugar_atoms:
						nearby_sugars.append(sugar)
						nearby_sugar_images.append(image)
		return [nearby_sugars, nearby_sugar_images]

	def is_deletable_OW(self,OW,i,j,k):
		[nearby_sugars, nearby_sugar_images] = self.find_nearby_sugars(i,j,k)
		for j,sugar in enumerate(nearby_sugars):
			distance = self.distance_calc(OW,sugar,nearby_sugar_images[j])
			if distance <= self.min_dist:
				return False
		return True	

	def distance_calc(self,atom1,atom2,atom2_image):
		distance2 = 0
		atom2_position = []
		for j in range(3):
			if atom2_image[j] == 0:
				atom2_position.append(atom2.position[j])
			else:
				atom2_position.append(atom2.position[j] + atom2_image[j]*self.unit_cell[j])
		for j in range(3):
			temp = atom2_position[j] - atom1.position[j]
			distance2 += temp*temp
		return math.sqrt(distance2)

	def create_deletable_atomgroup(self):
		self.deletable_resids = sorted(self.deletable_resids)
	 

		atom_ind = 0
		deletable_ind = 0
		cur_deletable_resid = self.deletable_resids[deletable_ind]
		n_atoms = SOL_sel.n_atoms
		n_deletable_resids = len(self.deletable_resids)
		while(True):
			cur_atom = SOL_sel[atom_ind]
			cur_resid = cur_atom.resid
			if cur_resid > cur_deletable_resid:
				deletable_ind = deletable_ind + 1
				cur_deletable_resid = self.deletable_resids[deletable_ind]
			if cur_resid == cur_deletable_resid:
				self.deletable_atoms.append(cur_atom)
			atom_ind = atom_ind + 1
			if atom_ind >= n_atoms or deletable_ind >= n_deletable_resids:
				break


		return mda.core.groups.AtomGroup(self.deletable_atoms)

	def random_deletion(self):
		k = len(self.deletable_resids) - n_delete
		self.non_deleted_waters = random.sample(self.deletable_resids, k)



	def create_post_deletion_atomgroup(self):
		self.random_deletion()
		self.final_non_deleted_waters = self.non_deleted_waters + self.non_deleted_resids
		self.final_non_deleted_waters = sorted(self.final_non_deleted_waters)
	 

		atom_ind2 = 0
		deletable_ind2 = 0
		if len(self.final_non_deleted_waters) !=0:
			cur_deletable_resid2 = self.final_non_deleted_waters[deletable_ind2]
			n_atoms2 = SOL_sel.n_atoms
			n_deletable_resids2 = len(self.final_non_deleted_waters)
			while(True):
				cur_atom2 = SOL_sel[atom_ind2]
				cur_resid2 = cur_atom2.resid
				if cur_resid2 > cur_deletable_resid2:
					deletable_ind2 = deletable_ind2 + 1
					if atom_ind2 >= n_atoms2 or deletable_ind2 >= n_deletable_resids2:
						break
					cur_deletable_resid2 = self.final_non_deleted_waters[deletable_ind2]
				if cur_resid2 == cur_deletable_resid2:
					self.remaining_waters.append(cur_atom2)
				atom_ind2 = atom_ind2 + 1
				if atom_ind2 >= n_atoms2 or deletable_ind2 >= n_deletable_resids2:
					break

		self.append_misc()
		self.append_sugars()
		return mda.core.groups.AtomGroup(self.final_system_sel)

	def append_sugars(self):
		for atom in ref_sel:
			self.remaining_waters.append(atom)
		self.final_system_sel = sorted(self.remaining_waters)

	def append_misc(self):
		for atom in misc_sel:
			self.remaining_waters.append(atom)

	def sugar_shell(self):
		self.non_deleted_resids = sorted(self.non_deleted_resids)
	 

		atom_ind3 = 0
		deletable_ind3 = 0
		cur_deletable_resid3 = self.non_deleted_resids[deletable_ind3]
		n_atoms3 = SOL_sel.n_atoms
		n_deletable_resids3 = len(self.non_deleted_resids)
		while(True):
			cur_atom3 = SOL_sel[atom_ind3]
			cur_resid3 = cur_atom3.resid
			if cur_resid3 > cur_deletable_resid3:
				deletable_ind3 = deletable_ind3 + 1
				cur_deletable_resid3 = self.non_deleted_resids[deletable_ind3]
			if cur_resid3 == cur_deletable_resid3:
				self.all_sugar_shells.append(cur_atom3)
			atom_ind3 = atom_ind3 + 1
			if atom_ind3 >= n_atoms3 or deletable_ind3 >= n_deletable_resids3:
				break

		for atom in ref_sel:
			self.all_sugar_shells.append(atom)
		self.all_sugar_shells = sorted(self.all_sugar_shells)
		return mda.core.groups.AtomGroup(self.all_sugar_shells)

parser = argparse.ArgumentParser()

parser.add_argument('-tol', help='Minimum distance between sugar and water (Angstrom)', required=True)
parser.add_argument('-s', help='Structure input file (.gro)', required=True)
parser.add_argument('-f', help='Trajectory input file (.xtc)', required=True)
parser.add_argument('-o', help='Output basename', required=True)
parser.add_argument('-nd', help='Number of Waters to delete', required =True)
parser.add_argument('-w', help='MDAnalysis group selection for solvent to be removed', required=True)
parser.add_argument('-ref', help='MDAnalysis group selection for reference for hydration shells', required=True)
parser.add_argument('-sol', help='DAnalysis group selection for solvent (if different from watername)',required=False)
parser.add_argument('-misc',help='DAnalysis group selection for other groups in system', required=False)

args = parser.parse_args()

## MDAnalysis Universe Declarations and Atom selection
u = mda.Universe(args.s,args.f)
# watername = 'resname ' + args.w
# refname = 'resname ' + args.ref
water_sel = u.select_atoms(args.w)
ref_sel = u.select_atoms(args.ref)
if args.sol is not None:
	SOL_sel = u.select_atoms(args.sol)
else:
	SOL_sel = u.select_atoms(args.w)
if args.misc is not None:
	misc_sel = u.select_atoms(args.misc)
else:
	misc_sel = u.select_atoms('')
# misc_sel = u.select_atoms('resname GRAP or resname GLY')

## Trajectory Selection
n_delete = int(args.nd)
my_test = my_system(u,float(args.tol))
my_test.find_deletable_waters()
print('number of deletable waters: ' + str(len(my_test.deletable_resids)))
ag = my_test.create_post_deletion_atomgroup()

############# IO
u2atoms = ag
# coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(),u2atoms).run().results
### DEPRECATED COMMANDS MAY BE NEEDED IN MDAnalysis versions below 2.0 ###
u2 = mda.Merge(u2atoms)
# u2.load_new(coordinates, format=MemoryReader)
with mda.Writer(args.o, reindex=True) as w:
	w.fmt['box_orthorhombic'].format(box=(u.trajectory[0]._unitcell[0],u.trajectory[0]._unitcell[1],u.trajectory[0]._unitcell[2]))
	w.write(u2.atoms)
