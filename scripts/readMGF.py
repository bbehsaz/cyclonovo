import sys, string

def readMGF(mgfFile):
	peaksFile = open(mgfFile, "r")
	cyclePath={}
	protonMass = 1.00728
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
	nline= int(200000000)
	peptideMass = {}
	peptidePrecursorMass = {}
	peptideNames = {}
	n=0
	condition ={}
	pepMasses={}
	mapantimartin2pep={}

	 
	theoreticalSpectrum={}

	#Reading the peaks and intensities.
	fileLines = {}
	charges= {}
	antimartinIDs= {}
	n= -1
	pval = 0
	antimartin = 1
	fileLines={}
	pvalues = {}
	retentions = {}
	inchiIDS= {}
	while(True):
		line = peaksFile.readline()
		if not line:
			break
		line = line.strip()
		if line == "BEGIN IONS":
			peptide= ""
			inchi = ""
			while(True):
				line = peaksFile.readline()
				if not line[0].isdigit():
					if line[0:6]=="CHARGE":
						charge = int(''.join(c for c in line.split("=")[1] if c.isdigit()))

					if line[0:6]=="PEPMAS":
						pepMass=float(line.strip().split()[0][8:])
					if line[0:5]=="DBID=":
						line= line.strip()
						if line[5:] !="":
							antimartin = line[15:]
						# fileLines[antimartin]="BEGIN IONS\n"
						# fileLines[antimartin]=fileLines[antimartin]+line+"\n"
					if line[0:15]=="DEREPLICATORPV=":
						line = line.strip().split("=")
						pval = line[1]
					if "TITLE=" in line:
						peptide=line.strip().split("=")[1]
					if "RTINSECONDS" in line:
						retention = line.split("=")[1].strip()
					elif "INCHI=InChI" in line:
						inchi = line.split("INCHI=InChI=")[1]

						# print newline
						# print newline.split("INCHI=InChI=")[1].count("N")
						
						# fileLines[antimartin]=fileLines[antimartin]+line
					# if line[0:6]=="TITLE=":
						# if peptide!= "":

						# 	peptide=line.strip()[6:]			
						# fileLines[antimartin]=fileLines[antimartin]+line
				else:
					break
			if antimartin == n-1:
				antimartin = str(n)
			n=n+1
			# print peptide
			if peptide == "":
				peptide=str(n)
			else:
				peptide = peptide + "_"+str(n) 

			# fileLines[antimartin]=fileLines[antimartin]+"CHARGE="+str(charge)+"\n"
			# fileLines[antimartin]=fileLines[antimartin]+"NAME="+str(peptide)+"\n"
			# fileLines[antimartin]=fileLines[antimartin]+"PEPMASS="+str(pepMass)+"\n"
			#this if statement is to only include cyclic peptides in cyclic peptide IDs.
			peptides[peptide]=[]
			retentions[peptide] = retention

			inchiIDS[peptide] = inchi
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

			# antimartinID= antimartinIDs[peptide]
			mapantimartin2pep[antimartin]= peptide
			while line!="END IONS":
				peakLine = line.strip().split()
				peakMass = round(float(peakLine[0]),3)
				intensity = float(peakLine[1])
				#peaksUnfiltered[peptide].append(peakMass)
				peaks[peptide].append(peakMass)
				peaksnIntensity[peptide][peakMass]=intensity
				
				if charge == 2:
					charge1PeakMass = (2*peakMass) -2 + protonMass
					peaks[peptide].append(charge1PeakMass)
					peaksnIntensity[peptide][peakMass]=intensity
					print pepmass
					print (2*pepMass) - 2+ protonMass
					pepMasses[peptide] = (2*pepMass) - 2+ protonMass
				line = peaksFile.readline().strip()
			if line == "END IONS":
				#print len(peaks[peptide])
				# if n == nline:
				# 	break
				continue

	# print peaks
	# print s
	# return peaks, peaksnIntensity, pepMasses
	return peaksnIntensity,peaks,pepMasses,antimartinIDs,anti2pep,pvalues,charges, retentions
