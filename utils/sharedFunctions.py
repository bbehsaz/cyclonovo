import sys, string


def getSpecializedBuildingBlocks(monomers,cyclonovoPath): 
	if monomers == "ripps":
		building_blocks_main = getBuildingBlocks(cyclonovoPath+"/configs/aa_polymer_cyclomasses_ripps.txt")
	elif monomers == "top25":
		building_blocks_main = getBuildingBlocks(cyclonovoPath+"/configs/aa_polymer_masses_25.txt")
	else:	
		building_blocks_main = getBuildingBlocks(cyclonovoPath+"/configs/aa_polymer_masses.txt")
	polymer_repeat_units = getBuildingBlocks(cyclonovoPath+"/configs/polymer_repeat_masses.txt")
	return building_blocks_main, polymer_repeat_units
def getBuildingBlocks(path2bbFile): 
# The function gets the building block files and return the required datastructures 
	digits = 3
	aminoMAsses = {}
	standardMasses = set()
	mass2name = {}
	with open(path2bbFile) as bbFile:
		for line in bbFile:
				linesplit = line.strip().split()
				aminoMAsses[linesplit[0]] = round(float(linesplit[3]),digits)
				standardMasses.add(round(float(linesplit[3]),digits))

				mass2name[round(float(linesplit[3]),digits)] = linesplit[0]
	return standardMasses


def writeOriginalSpectra(lines,outputMGF):
	for line in lines:
		outputMGF.write(line)




