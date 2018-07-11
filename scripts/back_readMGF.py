import sys, string

def readMGF(mgfFile):
	peaksFile = open(mgfFile, "r")
	cyclePath={}

	peptides = {}
	peaks = {}
	peaksnIntensity = {}
	distances ={}
	nline= int(40)
	peptideMass = {}
	peptidePrecursorMass = {}
	antimartinIDs= {}
	peptideNames = {}
	n=0
	condition ={}
	anti2pep = {}

	antimartin2pep={}
	nline= int(379)
	peptideMass = {}
	peptidePrecursorMass = {}
	peptideNames = {}
	n=0
	condition ={}
	pepMasses={}
	mapantimartin2pep={}

	 
	theoreticalSpectrum={}

	peaksFile = open(sys.argv[1], "r")

	#Reading the peaks and intensities.
	fileLines = {}
	charges= {}
	antimartinIDs= {}
	n= -1
	pval = 0
	antimartin = 1
	fileLines={}
	pvalues = {}
	while(True):
		line = peaksFile.readline()
		if not line:
			break
		line = line.strip()
		if line == "BEGIN IONS":
			while(True):
				line = peaksFile.readline()
				if not line[0].isdigit():
					if line[0:6]=="CHARGE":
						charge=int(line[7:8])
					if line[0:6]=="PEPMAS":
						pepMass=float(line.strip().split()[0][8:])
					if line[0:5]=="DBID=":
						line= line.strip()
						antimartin=line[15:]
						fileLines[antimartin]="BEGIN IONS\n"
						fileLines[antimartin]=fileLines[antimartin]+line+"\n"
					if line[0:15]=="DEREPLICATORPV=":
						line = line.strip().split("=")
						pval = line[1]
					if line[0:5]=="NAME=":
						peptide=line.strip()[5:]
						
						fileLines[antimartin]=fileLines[antimartin]+line

				else:
					break
			n=n+1
			fileLines[antimartin]=fileLines[antimartin]+"CHARGE="+str(charge)+"\n"
			fileLines[antimartin]=fileLines[antimartin]+"NAME="+str(peptide)+"\n"
			fileLines[antimartin]=fileLines[antimartin]+"PEPMASS="+str(pepMass)+"\n"
			#this if statement is to only include cyclic peptides in cyclic peptide IDs.
			peptides[peptide]=[]

			peaksnIntensity[peptide]={}
			peaks[peptide]=[]
			charges[peptide]=charge
			pepMasses[peptide]=pepMass
			antimartinIDs[peptide]=antimartin
			anti2pep[antimartin] =peptide
			pvalues[peptide] = pval
						# if antimartinIDs[peptide][0]=="0":
						# 	antimartinIDs[peptide]=line[16:]
			antimartin2pep[peptide]=antimartinIDs[peptide]

			antimartinID= antimartinIDs[peptide]
			mapantimartin2pep[antimartinID]= peptide
			while line!="END IONS":
				peakLine = line.strip().split()
				peakMass = float(peakLine[0])
				intensity = float(peakLine[1])
				#peaksUnfiltered[peptide].append(peakMass)
				peaks[peptide].append(peakMass)
				peaksnIntensity[peptide][peakMass]=intensity
				line = peaksFile.readline().strip()
			if line == "END IONS":
				if n == nline:
					break
				continue
	# print peaks
	# print peaksnIntensity
	# return peaks, peaksnIntensity, pepMasses
	return peaksnIntensity,peaks,pepMasses,antimartinIDs,anti2pep,pvalues
# print readMGF(sys.argv[1])
	# molDir = "/Users/bahar/workspace/antibiotic-sequencing/code/data/data_high/frags/"
	# nAA = {}
	# cyclePathID = {}
	# cyclepath= {}
	# #Read the real masses for each peptide and store in aaMasses[antimartinID]
	# for peptide in antimartinIDs:
	# 	#antimartinID = 28811 ##### THIS must iterate
	# 	# if condition[antimartinID]!="cyclic":
	# 	# 	continue
	# 	antimartinID = antimartinIDs[peptide]
	# 	structureFile = open(molDir+"antimartin_2012_"+antimartinID+".fragment.txt","r")
	# 	while(True):
	# 		line = structureFile.readline()
	# 		if not line:
	# 			break 
	# 		line = line.strip().split()
	# 		if line[2]=="components":
	# 			nNodes = int(line[4])
	# 			nAA[antimartinID]=nNodes
	# 			aaMasses[antimartinID] = [0]*nNodes
	# 			edges = {}
	# 			neighbors = {}
	# 			inneighbors={}
	# 			for i in range(nNodes):
	# 				line = structureFile.readline().strip().split()
	# 				#aaMasses[antimartinID][i] =float("%.1f" %  float(line[2]))
	# 				aaMasses[antimartinID][i] =float("%.2f" %  float(line[2]))
	# 				edges[i] = [0]*nNodes
	# 				neighbors[i]=[]
	# 				inneighbors[i]=[]
	# 		if line[2]=="bonds":
	# 			nEdges = int(line[4])
	# 			nCycles = nEdges - int(nNodes-1)
	# 			for i in range(nEdges):
	# 				line= structureFile.readline().strip().split()
	# 				edges[int(line[0])][int(line[2])]= 1
	# 				neighbors[int(line[0])].append(int(line[2]))
	# 				inneighbors[int(line[2])].append(int(line[0]))
	# 	outDegree = [0]*nNodes
	# 	inDegree = [0]*nNodes
	# 	outEnds = []
	# 	inEnds = []
	# 	inBranch = []
	# 	outBranch = []
	# 	nodeStatus=[""]*nNodes
	# 	condition[antimartinID] = "cyclic"
	# 	for i in range(nNodes):
	# 		if len(neighbors[i])!=1:
	# 			condition[antimartinID] = "branch-cyclic"

	# 	for i in range(nNodes):
	# 	# 	outDegree[i]= sum(edges[i])
	# 	# 	if outDegree[i]==0:
	# 	# 		outEnds.append(i)
	# 	# 		nodeStatus[i] = "outEnd"
	# 	# 		condition[antimartinID] = "branch-cyclic"
	# 	# 	if outDegree[i]>1:
	# 	# 		outBranch.append(i)
	# 	# 		nodeStatus[i]="outBranch"
	# 	# 		condition[antimartinID] = "branch-cyclic"
	# 	# 	for j in range(nNodes):
	# 	# 		x=inDegree[i]
	# 	# 		inDegree[i]=x+edges[j][i]
	# 	# 	if inDegree[i]==0:
	# 	# 		inEnds.append(i)
	# 	# 		nodeStatus="inEnd"
	# 	# 		condition[antimartinID] = "branch-cyclic"
	# 	# 	if inDegree[i]>1:
	# 	# 		inBranch.append(i)
	# 	# 		nodeStatus[i]="inBranch"
	# 	# 		condition[antimartinID] = "branch-cyclic"
	# 		# at this point we start generating linearize stuff
	# 		thSpectra = {}
	# 		longpath=[0]*nNodes
	# 		longpathID=[0]*nNodes
	# 		if condition[antimartinID]=="cyclic":
	# 			cyclePath[antimartinID]=[]
	# 			cyclePathID[antimartinID]=[]
	# 			cyclePath[antimartinID].append(0)
	# 			cyclePathID[antimartinID].append(0)
	# 			nextNode= neighbors[0][0]
	# 			i=0
	# 			while i<nNodes:
	# 				currenNode = nextNode
	# 				cyclePath[antimartinID].append(aaMasses[antimartinID][currenNode])
	# 				cyclePathID[antimartinID].append(currenNode)
	# 				nextNode=neighbors[currenNode][0]
	# 				i=i+1

