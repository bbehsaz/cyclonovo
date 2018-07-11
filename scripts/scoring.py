import sys,string


#this function generally scores a spectra against a genome mining result
#arguments: experimentalSpectrum, theoreticalSpectrum, error value
def scoreSpecAa(exSpec,thSpec,e):
	score = 0
	for exPeak in exSpec:
		for thPeak in thSpec:
			adjustedTHPeak = thPeak + 1.00728
			if abs(exPeak- adjustedTHPeak )< e:
				score += 1
	return score


def scorePeptideFoundFrags(possiblePeptide,exSpec,thspec,e): #possiblePeptide, is the peptide sequence as a list of masses
	foundTradnSatMasses={}

	foundTradMasses={}
	scoreTrad = 0
	scoreSat = 0
	numThSpec = len(set(thspec.values()))
	matrixThCheck = {}
	# startends = cyclicTheoreticalSpectrum(possiblePeptide,peptide)[1]
	correctDistances = 0
	if thspec == 0:
			return 0
	for i in range(min(200,len(exSpec))):
		peak = exSpec[i]
		check = 0
		e = 0.015

		if check==0:
			for mass in thspec.values():
			# for mass in thSpecWRep:
				adjustedMass = mass +1.00728
				if abs(peak-adjustedMass)<e:
					correctDistances=correctDistances+1
					scoreTrad = scoreTrad+1
					# print mass
					#peaks_annotations.write(str(i)+"\t"+ str(distance)+"\t"+str(distancesCleaned[peptide][distance])+"\t"+"traditional"+"\t"+str(mass)+"\t"+combinations[peptide][mass] + "\n")
					#peaks_annotations.write(str(i)+"\t"+str(sorted_peaks[peptide][i])+"\t"+str(peaksnIntensity[peptide][sorted_peaks[peptide][i]])+"\t"+str(mass)+"\t"+combinations[peptide][mass]+"\t"+"Traditional"+"\n")
					if mass not in foundTradnSatMasses:
						foundTradnSatMasses[mass] = 1
					if mass not in foundTradMasses:
						foundTradMasses[mass] = i+1
					check=1
					break
		# e = 0.015
		# peak = pepMasses[peptide]-sorted_peaks[peptide][i]
		# if check==0:
		# 	for mass in set(thspec.values()):
		# 		# adjustedMass = mass +1.00728
		# 		adjustedMass = mass
		# 		# print mass
		# 		if abs(peak-adjustedMass)<e:
		# 			correctDistances=correctDistances+1
		# 			scoreTrad = scoreTrad+1
		# 			# print mass
		# 			#peaks_annotations.write(str(i)+"\t"+ str(distance)+"\t"+str(distancesCleaned[peptide][distance])+"\t"+"traditional"+"\t"+str(mass)+"\t"+combinations[peptide][mass] + "\n")
		# 			#peaks_annotations.write(str(i)+"\t"+str(sorted_peaks[peptide][i])+"\t"+str(peaksnIntensity[peptide][sorted_peaks[peptide][i]])+"\t"+str(mass)+"\t"+combinations[peptide][mass]+"\t"+"Traditional"+"\n")
		# 			if mass not in foundTradnSatMasses:
		# 				foundTradnSatMasses[mass] = 1
		# 			check=1
		# 			break
		# peak = sorted_peaks[peptide][i] +18.01528
		# # check = 0
		# e= 0.015
		# if check==1:
		# 	for mass in list(thspec.values()):
		# 		adjustedMass = mass +1.00728
		# 		# adjustedMass = mass
		# 		if abs(peak-adjustedMass)<e:
		# 			correctDistances=correctDistances+1
		# 			scoreTrad = scoreTrad+1
		# 			# print mass
		# 			#peaks_annotations.write(str(i)+"\t"+ str(distance)+"\t"+str(distancesCleaned[peptide][distance])+"\t"+"traditional"+"\t"+str(mass)+"\t"+combinations[peptide][mass] + "\n")
		# 			#peaks_annotations.write(str(i)+"\t"+str(sorted_peaks[peptide][i])+"\t"+str(peaksnIntensity[peptide][sorted_peaks[peptide][i]])+"\t"+str(mass)+"\t"+combinations[peptide][mass]+"\t"+"Traditional"+"\n")
		# 			if mass not in foundTradnSatMasses:
		# 				foundTradnSatMasses[mass] = 1

		# 			check=1
		# 			break
		# peak = sorted_peaks[peptide][i] + 17.02655
		# # check = 0
		# if check==0:
		# 	for mass in list(thspec.values()):
		# 		adjustedMass = mass +1.00728
		# 		# adjustedMass = mass
		# 		if abs(peak-adjustedMass)<e:
		# 			correctDistances=correctDistances+1
		# 			scoreTrad = scoreTrad+1
		# 			# print mass
		# 			#peaks_annotations.write(str(i)+"\t"+ str(distance)+"\t"+str(distancesCleaned[peptide][distance])+"\t"+"traditional"+"\t"+str(mass)+"\t"+combinations[peptide][mass] + "\n")
		# 			#peaks_annotations.write(str(i)+"\t"+str(sorted_peaks[peptide][i])+"\t"+str(peaksnIntensity[peptide][sorted_peaks[peptide][i]])+"\t"+str(mass)+"\t"+combinations[peptide][mass]+"\t"+"Traditional"+"\n")
		# 			if mass not in foundTradnSatMasses:
		# 				foundTradnSatMasses[mass] = 1

		# 			check=1
		# 			break
		# peak = sorted_peaks[peptide][i] + 28.00615
		# if check==0:
		# 	for mass in list(thspec.values()):
		# 		adjustedMass = mass +1.00728
		# 		# adjustedMass = mass
		# 		if abs(peak-adjustedMass)<e:
		# 			correctDistances=correctDistances+1
		# 			scoreTrad = scoreTrad+1
		# 			# print mass
		# 			#peaks_annotations.write(str(i)+"\t"+ str(distance)+"\t"+str(distancesCleaned[peptide][distance])+"\t"+"traditional"+"\t"+str(mass)+"\t"+combinations[peptide][mass] + "\n")
		# 			#peaks_annotations.write(str(i)+"\t"+str(sorted_peaks[peptide][i])+"\t"+str(peaksnIntensity[peptide][sorted_peaks[peptide][i]])+"\t"+str(mass)+"\t"+combinations[peptide][mass]+"\t"+"Traditional"+"\n")
		# 			if mass not in foundTradnSatMasses:
		# 				foundTradnSatMasses[mass] = 1

		# 			check=1
		# 			break
		# peak = sorted_peaks[peptide][i] + 1.00728
		# check = 0
		# if check==0:
		# 	for mass in list(thspec.values()):
		# 		adjustedMass = mass +1.00728
		# 		# adjustedMass = mass
		# 		if abs(peak-adjustedMass)<e:
		# 			correctDistances=correctDistances+1
		# 			scoreTrad = scoreTrad+1
		# 			# print mass
		# 			#peaks_annotations.write(str(i)+"\t"+ str(distance)+"\t"+str(distancesCleaned[peptide][distance])+"\t"+"traditional"+"\t"+str(mass)+"\t"+combinations[peptide][mass] + "\n")
		# 			#peaks_annotations.write(str(i)+"\t"+str(sorted_peaks[peptide][i])+"\t"+str(peaksnIntensity[peptide][sorted_peaks[peptide][i]])+"\t"+str(mass)+"\t"+combinations[peptide][mass]+"\t"+"Traditional"+"\n")
		# 			if mass not in foundTradnSatMasses:
		# 				foundTradnSatMasses[mass] = 1

		# 			check=1
		# 			break
	matrixThCheck = {}
		# print "----"
	matrixIntensity = {}
	for k in range(len(possiblePeptide)):
		matrixThCheck[k] = [0]*len(possiblePeptide)
		for j in range(len(possiblePeptide)):
			# print j
			# if k == j:
			# 	matrixThCheck[k].append(0)
			# 	continue
			if thspec[(k,j)] in foundTradMasses:
				matrixThCheck[k][j] = thspec[(k,j)]
				# matrixThCheck[k][j] = foundTradMasses[thspec[(k,j)]]
	for k in range(len(possiblePeptide)):
		matrixIntensity[k] = [0]*len(possiblePeptide)
		for j in range(len(possiblePeptide)):
			# print j
			# if k == j:
			# 	matrixThCheck[k].append(0)
			# 	continue
			if thspec[(k,j)] in foundTradMasses:
				# matrixThCheck[k][j] = thspec[(k,j)]
				matrixIntensity[k][j] = foundTradMasses[thspec[(k,j)]]
	tradMatches = foundTradMasses
	if scoreTrad/(len(set(thspec.values()))*1.000) == 1:
		print len(set(thspec.values()))
		print scoreTrad
	return scoreTrad/(len(set(thspec.values()))*1.000), tradMatches, numThSpec,matrixThCheck,matrixIntensity