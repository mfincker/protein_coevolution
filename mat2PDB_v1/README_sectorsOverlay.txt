READ ME for sectorsOverlay

sectorsOverlay will take a sector database divided by protein and will produce a pdb file of all the sectors.  The individual sectors will be colored by B-factor and can be separated in pyMOL by coloring by Bfactor.

Steps to using sectorsOverlay:
1. Get sector database of proteins
	a) In sectorhelperfunctions there is a function called “getSectorsByProtein”, apply that to a sector database (from a group or subgroup).  This will give a new database dividided by protein

2. Use sectorsOverlay
	!!!! be sure not to use sectorOverlay this is v1 !!!!
	a) sectorsOverlay will take a database of a given protein
	* e.g.) from subgroup of toxins we get a 1x6 database of proteins
	sectorsOverlay(toxins{1}) will give a pdb of the sectors from the 	first protein in the database. 

3. Visualize in pymol
	1. type command ‘fetch’  pdbid
	e.g) fetch 7HAL
	2. open the file of sectors 
	e.g.) 7HAL_sectors
	
	2. type command: hide everything
	3. type command: show cartoon
	4. on the right there will be 2 files: one will be the pdb and the sectors
	5. color sector file: C (click on the right hand panel) and then click spectrum >> by factor
	6. show sectors : S in GUI dots


—— EXAMPLE code—

sectorsbyProtein = getSectorsByProtein(subgroupsectorDB);

pdboutput from sectorsOverlay(sectorsbyProtein{1})

*open pymol*

see step 3


	