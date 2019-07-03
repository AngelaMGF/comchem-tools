# comchem-tools
*Toolkit for Computational Chemistry*

Useful tools to ease Computational Chemistry tasks.

 - You will need to have Python installed to work with these (developed with Python 3.7; also tested in Python 2.7)

# Files

### PDB to CONFIG
Given an input file with `.pdb` extension, converts content into a CONFIG file. Two types of pdb can be given:
	1 - RASMOL file created from the tool cfg2pdb or rvc2pdb (not yet ready)
	2 - pdb file generated by Mercury from a CIF file.

Usage:

* Option 1:

`python pdb2cfg.py input_file_name`

* Option 2:
`chmod +x pdb2cfg.py ` (Make it executable)
`./pdb2cfg.py input_file_name`

Where `input_file_name` is the `.pdb` (with or without extension).

It requires an arbitrary number of `.xyz` files in the same folder.


Usage:

`./pdb2cfg.py input_file_name`

Where `input_file_name` is the `.pdb` (with or without extension).


Requires, in the same folder, a number of `.xyz` files equal to the number of species in the .pdb file:
 - Must contain the OPLS-AA tags of each species, respecting the standart .xyz format (XMOL) and the order of the .pdb
 - Each file must be named like `1.xyz` ; `2.xyz` ; .... with the order in which the species are presented in the .pdb file
