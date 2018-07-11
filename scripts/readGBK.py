import sys, string

#THIS is right now only for one peptide

def readGBK(gbkFileAddress):
	aaMassesFile = open("MASSES.txt","r")
	aaMassesFile.readline()
	aaMasses = {}
	while(True):
		line = aaMassesFile.readline()
		if not line:
			break
		aaMasses[line.split()[0]] = float(line.split()[5])
	antiSmashPredicts ={}
	antiSmashPredictsScore = {}
	pos = 0
	gbkFile = open(gbkFileAddress,"r")
	while(True):
		line =gbkFile.readline().strip()
		if not line:
			break
		predictions = line.split()[2]
		aminoAcids = predictions.split(";")
		antiSmashPredicts[pos] = []
		antiSmashPredictsScore[pos] = []
		for i in range(5):
			aminoScore = aminoAcids[i]
			aminoName = aminoScore.split("(")[0]
			score = float(aminoScore.split("(")[1].split(")")[0])
			antiSmashPredicts[pos].append(aaMasses[aminoName])
			antiSmashPredictsScore[pos].append(score)
		pos += 1

	return antiSmashPredicts,antiSmashPredictsScore


#print readGBK(sys.argv[1])[0]
