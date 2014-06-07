########################################################
# 06/06/2014 - Maeva
# Python script to parse the membraneDB.xml files
# containing the pdb IDs and information of all crystallized
# membrane protein.
# The script creates a nested dictionnary with the protein 
# group as first key, protein subgroup as second key and protein
# pdb ID as third one. The fourth key is then 'pdb', 'group',
# 'subgroup', 'memberProteins' or 'relatedProteins'.
# Once the dictionnary has been created, it is written to file
# where every group of 5 lines corresponds to:
# 	- pdbID
#	- group
#	- subgroup
#	- member proteins
#	- related proteins
# See ./uniqueMembraneProteinsDB.txt for an example of the output.

import re

#############
# FUNCTIONS #
#############

# Extract the name of the xml level currently parsed 
# and passed as input argument
def extractName(text):
	pattern = '<name>.*</name>'
	name = re.search(pattern, text)
	if name == None:
		name = 'None'
	else:
		name = name.group()[6:-7]
	return name

# Extract the pdb code of the xml level currently parsed 
# and passed as input argument
def extractPdbCode(text):
	pattern = '<pdbCode>.*</pdbCode>'
	pdbId = re.search(pattern, text)
	if pdbId == None:
		pdbId = 'None'
	else:
		pdbId = pdbId.group()[9:-10]
	return pdbId

# Extract all pdb code present at the xml level currently parsed,
# and passed as input argument
def extractAllPdbCode(text):
	pdbIdList = []
	pattern = '<pdbCode>.*</pdbCode>'
	pdbId = re.findall(pattern, text)
	if pdbId == None:
		pdbId = 'None'
	else:
		for pdb in pdbId:
			pdbIdList.append(pdb[9:-10])
	return pdbIdList


##########
# SCRIPT #
##########

# Open file
filename = './membraneDB.xml';


f = open(filename, 'r')

lines = f.read()

membraneDB = {}
uniqueProteinPdbs = []
allProteinPdbs = []

# Extract groups
groups = lines.split('<group>')
for group in groups:
	# Extract group name
	groupName = extractName(group)
	print(groupName)


	if groupName != 'None':
		membraneDB[groupName] = {}

		# Extract subgroup name
		subgroups = group.split('<subgroup>')
		subgroups = subgroups[1:] # Remove the group name
		for subgroup in subgroups:
			subgroupName = extractName(subgroup)
			print('\t' + subgroupName)

			if subgroupName != 'None':
				membraneDB[groupName][subgroupName] = {}

				# Extract proteins
				proteins = subgroup.split('<protein>')
				proteins = proteins[1:] # Remove subgroup name split part
				for protein in proteins:
					proteinPdb = extractPdbCode(protein)

					if proteinPdb != 'None' and not(proteinPdb in uniqueProteinPdbs):
						membraneDB[groupName][subgroupName][proteinPdb] = {}
						membraneDB[groupName][subgroupName][proteinPdb]['pdbId'] = proteinPdb
						membraneDB[groupName][subgroupName][proteinPdb]['memberProteins'] = []
						membraneDB[groupName][subgroupName][proteinPdb]['relatedProteins'] = []
						uniqueProteinPdbs.append(proteinPdb)
						allProteinPdbs.append(proteinPdb)

						# Extract member proteins
						memberProteins = protein.split('<memberProtein>')

						# Extract related proteins in the protein def
						relatedProteins = memberProteins[0].split('<relatedPdbEntries>')
						if len(relatedProteins) > 1:
							relatedPdbEntries = extractAllPdbCode(relatedProteins[1])
							for pdb in relatedPdbEntries:
								if not (pdb in allProteinPdbs):
									allProteinPdbs.append(pdb)
									membraneDB[groupName][subgroupName][pdb] = {}
									membraneDB[groupName][subgroupName][pdb]['pdbId'] = pdb
									membraneDB[groupName][subgroupName][pdb]['memberProteins'] = []
									membraneDB[groupName][subgroupName][pdb]['relatedProteins'] = []
								membraneDB[groupName][subgroupName][proteinPdb]['relatedProteins'].append(pdb)




						memberProteins = memberProteins[1:]
						for memberProtein in memberProteins:
							memberProteinPdb = extractPdbCode(memberProtein)
							if memberProteinPdb != 'None':
								if not (memberProteinPdb in allProteinPdbs):
									allProteinPdbs.append(memberProteinPdb)
									membraneDB[groupName][subgroupName][memberProteinPdb] = {}
									membraneDB[groupName][subgroupName][memberProteinPdb]['pdbId'] = memberProteinPdb
									membraneDB[groupName][subgroupName][memberProteinPdb]['memberProteins'] = []
									membraneDB[groupName][subgroupName][memberProteinPdb]['relatedProteins'] = []
								membraneDB[groupName][subgroupName][proteinPdb]['memberProteins'].append(memberProteinPdb)
						
							# Extract related proteins in the memberProtein
							relatedProteins = memberProtein.split('<relatedPdbEntries>')
							if len(relatedProteins) > 1:
								relatedPdbEntries = extractAllPdbCode(relatedProteins[1])
								for pdb in relatedPdbEntries:
									if not (pdb in allProteinPdbs):
										allProteinPdbs.append(pdb)
										membraneDB[groupName][subgroupName][pdb] = {}
										membraneDB[groupName][subgroupName][pdb]['pdbId'] = pdb
										membraneDB[groupName][subgroupName][pdb]['memberProteins'] = []
										membraneDB[groupName][subgroupName][pdb]['relatedProteins'] = []
									membraneDB[groupName][subgroupName][proteinPdb]['relatedProteins'].append(pdb)


print('Unique proteins: %i', len(uniqueProteinPdbs))
print('All proteins: %i', len(allProteinPdbs))

f.close()


###############
# WRITE FILES #
###############

# Write down all proteins
fall = open('./allMembraneProteinsDB.txt','w')
# Write down unique protein (no doulbe entry of a pdb ID)
funique = open('./uniqueMembraneProteinsDB.txt','w')
# Group files

groupNames = membraneDB.keys()
for group in groupNames:
	filename = './' + group + 'DB.txt'
	fgroup = open(filename,'a')

	subgroupNames = membraneDB[group].keys()
	for subgroup in subgroupNames:
		proteinIds = membraneDB[group][subgroup].keys()
		for ID in proteinIds:
			
			if ID in uniqueProteinPdbs:
				funique.write(membraneDB[group][subgroup][ID]['pdbId'] + '\n')
				funique.write(group + '\n')
				funique.write(subgroup + '\n')
				funique.write('\t'.join(membraneDB[group][subgroup][ID]['memberProteins']) + '\n')
				funique.write('\t'.join(membraneDB[group][subgroup][ID]['relatedProteins']) + '\n')


			fall.write(membraneDB[group][subgroup][ID]['pdbId'] + '\n')
			fall.write(group + '\n')
			fall.write(subgroup + '\n')
			fall.write('\t'.join(membraneDB[group][subgroup][ID]['memberProteins']) + '\n')
			fall.write('\t'.join(membraneDB[group][subgroup][ID]['relatedProteins']) + '\n')

			fgroup.write(membraneDB[group][subgroup][ID]['pdbId'] + '\n')
			fgroup.write(group + '\n')
			fgroup.write(subgroup + '\n')
			fgroup.write('\t'.join(membraneDB[group][subgroup][ID]['memberProteins']) + '\n')
			fgroup.write('\t'.join(membraneDB[group][subgroup][ID]['relatedProteins']) + '\n')
fall.close()
funique.close()





	