#!/usr/bin/env python

import sys
import math
from os import listdir
from os.path import isfile, join

def get_xyz():
	"""
	Returns the number of .xyz files in the current directory
	"""
	onlyfiles = [f for f in listdir(".") if isfile(join(".", f))]
	xyz = 0
	for f in onlyfiles:
		if ".xyz" in f:
			xyz += 1
	return xyz

def get_atoms_in_species(n):
	"""
	Returns a dictionary with the number of atoms in each species: {species_n : x atoms}
	"""
	all_species = {}
	for i in range(1, n+1):
		with open(str(i) + ".xyz") as f:
			number_atoms = int(f.readline())
			all_species[i] = number_atoms
	return all_species

def get_molecules_from_species(lines, start, first_atom, n):
	"""
	Returns the number of molecules of each species given in the input file and the number of the line where the next species begins
	"""
	num_molecs = 0
	other_species_index = 0
	for i in range(start, len(lines), n):
		if first_atom in lines[i].split()[2]:
			num_molecs += 1
		else:
			first_other_species = i
			break
	return num_molecs, first_other_species


def remove_lines(lst, unwanted_str):
	"""
	Given a list, it removes any lines that contain the substring given as the second argument of the function.
	> Used to remove lines with ANISOU
	"""
	new_lst = []
	for l in lst:
		if not unwanted_str in l:
			new_lst.append(l)
	return new_lst

def get_box_params(crystal):
	"""
	It has a lot of equations and prints, prettily and with the nice spaces in between, three lines with the cell vectors.
	Returns a string
	"""
	mult_a = math.ceil(32.2 / float(crystal[1]))
	mult_b = math.ceil(32.2 / float(crystal[2]))
	mult_c = math.ceil(32.2 / float(crystal[3]))

	a = float(crystal[1]) * mult_a
	b = float(crystal[2]) * mult_b
	c = float(crystal[3]) * mult_c

	alpha = math.radians(float(crystal[4]))
	beta = math.radians(float(crystal[5]))
	gamma = math.radians(float(crystal[6]))

	ba = b * math.cos(gamma)
	bb = b * math.sin(gamma)
	ca = c * math.cos(beta)
	cb =( c / math.sin(gamma) ) * ( math.cos(alpha) - math.cos(beta) * math.cos(gamma) )
	cc =( c / math.sin(gamma) ) * math.sqrt( 1 - math.cos(alpha)**2 - math.cos(beta)**2 - math.cos(gamma)**2 + 2 * math.cos(alpha)* math.cos(beta) * math.cos(gamma) )

	box_params = "{:20.12f}{:20.12f}{:20.12f}\n{:20.12f}{:20.12f}{:20.12f}\n{:20.12f}{:20.12f}{:20.12f}\n".format(a, 0, 0, ba, bb, 0, ca, cb, cc)

	return box_params

def change_tag(n, i, sizes):
	"""
	Returns the tag for a certain atom, as it is in the n.xyz file
	"""
	with open(str(n) + ".xyz") as f:
		lines = f.read().splitlines()
		atom = lines[i + 2]
	return atom.split()[0]

def get_atoms(lines, sizes, molecules):
	"""
	It changes the original tag with the one in the xyz files.
	It returns information related to it.
	Returns ordered list with strings with wanted format
	"""
	atoms = {}
	n = 1
	current_s = sizes[n] * molecules[n] 
	current_a = 0

	for i in range(0, len(lines)):

		if current_s <= 0:
			n = n + 1
			if not n in atoms:
				break
			current_s = sizes[n] * molecules[n]
			current_a = 0
		if current_a >= sizes[n]:
			current_a = 0

		l_sep = lines[i].split()
		atom_id = int(l_sep[1])
		x = float(l_sep[5])
		y = float(l_sep[6])
		z = float(l_sep[7])
		
		atom_tag = change_tag(n, current_a, sizes[n])
		atoms[atom_id] = {"x": x, "y": y, "z": z, "tag": atom_tag}
		current_s = current_s - 1
		current_a = current_a + 1
	return atoms

def get_easy_atoms(lines):
	atoms = {}
	for i in range(3, len(lines)):
		l_sep = lines[i].split()

		atom_id = int(l_sep[1])
		atom_tag = l_sep[2]
		x = float(l_sep[3])
		y = float(l_sep[4])
		z = float(l_sep[5])

		atoms[atom_id] = {"x": x, "y": y, "z": z, "tag": atom_tag}
	return atoms

def write_config(box_params, atoms):
	"""
	Creates / overwrites a CONFIG file in the same folder with the expected content. 
	It respects the given format.
	"""
	f = open("CONFIG", "w")
	f.write(" PBD\n")
	f.write("         0         3\n")
	f.write(box_params)
	for i in range(1, len(atoms)+1):
		f.write("{:<9}{:>9}\n".format(atoms[i]["tag"], i))
		f.write("{:20.12f}{:20.12f}{:20.12f}\n".format(atoms[i]["x"], atoms[i]["y"], atoms[i]["z"]))
	f.close()

def main():

	input_name = ""
	if len(sys.argv) == 1:
		onlyfiles = [f for f in listdir(".") if isfile(join(".", f))]
		for f in onlyfiles:
			if ".pdb" in f:
				input_name = f
	else:
		input_name =  sys.argv[1]
		if not ".pdb" in input_name:
			input_name = input_name + ".pdb"

	with open(input_name) as f:

		file = f.read()
		lines = file.splitlines()

		if not "CRYST" in file:
			box_params = "{:20.12f}{:20.12f}{:20.12f}\n{:20.12f}{:20.12f}{:20.12f}\n{:20.12f}{:20.12f}{:20.12f}\n".format(0, 0, 0, 0, 0, 0, 0, 0, 0)
			atoms = get_easy_atoms(lines)

		else:
			#Get number of species and respective number of atoms
			num_species = get_xyz()
			species_size = get_atoms_in_species(num_species)
			
			#Clean input
			aux_lines = remove_lines(lines, "ANISOU")
			input_lines = remove_lines(aux_lines, "CONECT")

			#Number of molecules of each species, stored in dictionary. 
			start = 5
			species_molecs = {}
			for s in range(1, num_species + 1):
				first_atom = input_lines[start].split()[2]
				species_molecs[s], start = get_molecules_from_species(input_lines, start, first_atom, species_size[s])

			#Box parameters
			crystal = input_lines[1].split()
			box_params = get_box_params(crystal)

			#Atom tags
			len_aux = len(input_lines) - 2
			atoms = get_atoms(input_lines[5:len_aux], species_size, species_molecs)
			
		#Writes the CONFIG file
		write_config(box_params, atoms)
		print("Success: CONFIG file completed, based on {}".format(input_name)) 

if __name__ == '__main__':
	main()