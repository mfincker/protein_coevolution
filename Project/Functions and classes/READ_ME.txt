————————————————————————————————————————————————
06/09/14 - Functions and classes folder Read_me
————————————————————————————————————————————————

This folder contains the functions and classes used in making the sector database or running or analysis. 

Content:
	- @Sector: contains the class description of the Sector class, used to represent a sectors and store important information.

	- Coordinates PCA: function to run PCA on the x,y,z coordinates of the centroids of the residues in a single sector as well as running PCA on all sectors contained in a sector database (cell array of sector objects).

	- mat2PDB_v1: contains the function mapping sector coordinates to its protein structure.

	- MISC wrapper: folder containing all functions pertaining to the SCA analysis

	- Phylogenetic distance: contains the function that, given a MSA and a reference sequence, will calculate the phylogenetic distance between all sequences in the MSA and the reference sequence.

	- Residue_PCA: contains functions counting amino acid counts and frequencies in a sector database as well as a function carrying PCA on a 20xn matrix containing the amino acid counts or frequencies in the n sectors from the database.

	- Sector_helper_functions: contains all functions necessary to run the Sector class constructor.