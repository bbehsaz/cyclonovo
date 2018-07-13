#!/usr/bin/env python

import sys
import sys, string, random
from scripts.readMGF import readMGF
from scripts.readStream import readStream
from operator import itemgetter
from itertools import combinations 



def generateAllAutconv(
	peaks,e,thresh,representative,standardAminoMasses,precursorMass,charge): 
###################################################################
#The function generate all autoconvolutiosn from peaks. 
#Arguments: 1) peaks: a dictionary with peak masses as keys and intensities as values. 2) error value. 3) thresh: The maximum number of peaks to be considered. 4) 
###################################################################
	sorted_peaks = sorted(peaks, key=lambda k: peaks[k], reverse=True ) #sort peaks based on intensity
	tempstnd = list(standardAminoMasses)[:]
	standardAminoMasses = set([round(x,3) for x in tempstnd])
	# minIntensity = peaks[sorted_peaks[len(sorted_peaks)-int(round(len(sorted_peaks)/4.0,0))]]
	minIntensity =0
	digits = 3
	# print minIntensity
	# print([x for x in sorted_peaks if peaks[x]<minIntensity])
	distances={} #includes the distances with their multiplicities
	tuples = {} # pairs of peaks resulted in a specific autoconvolution (using the rank of the peak when the peaks are sortd)
	tuplesMasses = {} # pairs of peaks resulted in a specific autoconvolution (using the mass of the peaks)

	allExistingDistances= []
	for i in range(min(thresh,len(peaks))):

		for j in range(min(thresh,len(peaks))):
			if  peaks[sorted_peaks[j]] < minIntensity:
				break
			distance = round(sorted_peaks[j]-sorted_peaks[i],3)
			if distance>300 or distance<40:
				continue
			distanceRounded = round(distance, digits)
			if distanceRounded >20:
				allExistingDistances.append(distanceRounded)
			if distanceRounded > 20 and  distanceRounded not in distances:
				tuples[distanceRounded]= []
				tuplesMasses[distanceRounded]= []
				tuples[distanceRounded].append((i,j))
				tuplesMasses[distanceRounded].append((sorted_peaks[j],sorted_peaks[i]))
				distances[distanceRounded] = 1
			elif distanceRounded > 20:
				tuples[distanceRounded].append((i,j))
				tuplesMasses[distanceRounded].append((sorted_peaks[j],sorted_peaks[i]))
				distances[distanceRounded] = distances[distanceRounded] + 1

	distancesSorted = {}
	standardAutconvPeaks = {}
	standardAutconvCleaned = {}
	# print sorted(distances.items(), key = itemgetter(1), reverse=True)
	autconvGroups = {}
	groupNum = 0
	clusterPeaks = {}
	for distance1 in sorted(distances.keys()): #Find clusters among distances
		distance = round(distance1,digits)
		groupFound = False
		for mass in standardAminoMasses:
			if abs(mass-distance)<(5)*e:
				for group in autconvGroups:
					if abs(autconvGroups[group][0][0]-distance)<e or abs(autconvGroups[group][0][1]-distance)<e:
						if autconvGroups[group][0][0] > distance:
							autconvGroups[group][0][0] = distance
						if autconvGroups[group][0][1] < distance:
							autconvGroups[group][0][1] = distance
						clusterPeaks[group] += tuplesMasses[distance] 
						autconvGroups[group][1].append(distance)
						groupFound=True
				if not groupFound:
					autconvGroups[groupNum] = ([distance,distance],[distance])
					clusterPeaks[groupNum] = tuplesMasses[distance] 
					#this dictionary has numbers as keys (starting from 0 up)
					# and each group has a tuple, first component is a list of size with minimum number being the first entry and the highest entry is the second component (basically the range of the group)
					#the second component of the value tuple is the list of all distances that appeared in that range.
					groupNum += 1
				break
	averageDistances = {}
	medianDistances = {}
	clusterPeaksCleaned = {}
	for group in autconvGroups:

		totalSum = 0
		numElements = 0
		allPeaksGroups = []
		for distance in autconvGroups[group][1]:
			totalSum += distances[distance]*distance
			numElements += distances[distance]

		# averageDist = round(totalSum/(numElements*1.0000),3)
		# averageDistances[averageDist] = numElements
		medianDist = autconvGroups[group][1][int(len(autconvGroups[group][1])/2)]
		medianDistances[medianDist] = numElements
		clusterPeaksCleaned[medianDist] = clusterPeaks[group]
	if representative == "average":
	#### GENERATING THE CLUSTER REPRESENTATIVES USING AVERAGE OF THE MEMBERS OF EACH CLUSTER
		averageDistancesSorted = sorted(averageDistances.items(), key=itemgetter(1), reverse=True)
		if averageDistancesSorted[0][1]<7 or averageDistancesSorted[1][1]<7 :
			return False

		filteredTuples = []
		# print averageDistances
		filteredDistances = {}
		for tuplevalue in averageDistancesSorted:
			if tuplevalue[1]>2:
				filteredTuples.append(tuplevalue)
				filteredDistances[tuplevalue[0]] = tuplevalue[1]
	else:
		## GENERATING THE CLUSTER REPRESENTATIVES USING MEDIAN OF THE MEMBERS OF EACH CLUSTER
		medianDistancesSorted = sorted(medianDistances.items(), key=itemgetter(1), reverse=True)
		filteredTuples = []
		filteredDistances = {} # THIS DICTIONARY HOLDS THE FINAL CLUSTERS
		filteredPeaks = {}
		for tuplevalue in medianDistancesSorted:
			if tuplevalue[1]>1:
				filteredTuples.append(tuplevalue)
				filteredPeaks[tuplevalue[0]] = clusterPeaksCleaned[tuplevalue[0]]
				filteredDistances[tuplevalue[0]] = tuplevalue[1]

	ALLclusterPeaksCleaned = {}
	ALLcleanedAutoconvolutions = {}
	clustersPeaksFiltered = filteredPeaks
	ALLcleanedAutoconvolutionsPeaks = {}
	for mass in clustersPeaksFiltered:
		clusterNum = 0
		ALLclusterPeaksCleaned = {}
		for i in range(0,len(clustersPeaksFiltered[mass])):
			peakPair1 = sorted(clustersPeaksFiltered[mass])[i]
			matched = False
			for j in range(0,len(clustersPeaksFiltered[mass])):
			# for j in range(len(standardAutconvPeaks[mass])):
				peakPair2 = sorted(clustersPeaksFiltered[mass])[j]
				if abs(peakPair2[1]- peakPair1[1])<29 and abs(peakPair2[0]- peakPair1[0])<29:
					found = False
					matched = True
					
					for cluster in ALLclusterPeaksCleaned:
						if peakPair1 in ALLclusterPeaksCleaned[cluster]:
							ALLclusterPeaksCleaned[cluster][peakPair2] = 1
							found = True
						if peakPair2 in ALLclusterPeaksCleaned[cluster]:
							ALLclusterPeaksCleaned[cluster][peakPair1] = 1
							found = True			
					if not found:
						clusterNum +=1
						ALLclusterPeaksCleaned[clusterNum] = {}
						ALLclusterPeaksCleaned[clusterNum][peakPair1] =1
						ALLclusterPeaksCleaned[clusterNum][peakPair2] =1
			if not matched:
				clusterNum +=1
				ALLclusterPeaksCleaned[clusterNum] = {}
				ALLclusterPeaksCleaned[clusterNum][peakPair1] =1

		ALLcleanedAutoconvolutionsPeaks[mass] = ALLclusterPeaksCleaned
		ALLcleanedAutoconvolutions[mass] = len(ALLclusterPeaksCleaned)

	finalClusters = {} #will contain all the clusters after satellites are removed with their max and min distances and their diamaeters

	for dis in sorted(ALLcleanedAutoconvolutions.items(), key= itemgetter(1), reverse=True):
		if dis[1]>0: 
			all_pair_peaks = []
			peaksPairs = ALLcleanedAutoconvolutionsPeaks[dis[0]]
			allclusterdistances = []
			all_pair_peaks = []
			for i in peaksPairs:
				for cluster in peaksPairs:
					for pair in peaksPairs[cluster]:
						all_pair_peaks.append(pair)
						allclusterdistances.append(round(abs(pair[0]-pair[1]),3))

			minDis = min(allclusterdistances)
			maxDis = max(allclusterdistances)
			distance = dis[0]
			
			finalClusters[distance] = (dis[1], round(maxDis- minDis,3)/2.0, minDis, maxDis, all_pair_peaks) #tuple with the frequency as the first component and the diameter of the cluster as second
			# print str(dis[0])+ "\t" + str(dis[1]) + "\t"+ str(round(maxDis- minDis,3) ) + "\t" + str(minDis) + "\t" + str(maxDis) 

	return finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge

def findSpectraPassingThresh(
	finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge,standardAminoMasses,e,alpha,betta):
	standardAutconvPeaks = {}
	standardAutconvCleaned = {}
	stand_aa_peak_pairs = {}
	for mass in standardAminoMasses:
		standardAutconvCleaned[mass] = 0
		# peaksIndexAminoDistances[mass] = []
		standardAutconvPeaks[mass] = []
		stand_aa_peak_pairs[mass] = []
	N=30
	# minFreq = sorted(ALLcleanedAutoconvolutions.values(), reverse=True)[3]
	totalIncluded = 0
	numSTDFound = 0
	prevValue =0
	topNconvolutions = {}
	relatedClusters= {}
	peaks_in_cluster = {}
	for mass in standardAminoMasses:
		relatedClusters[mass] = []
	sortedConvolutios = sorted(ALLcleanedAutoconvolutions.items(), key= itemgetter(1), reverse=True)
	# thresholdValue = 0
	protonMass = 1.00728
	pepMass = precursorMass*charge - protonMass
	import math
	thresholdValue = round(math.ceil(alpha *pepMass + betta))
	# thresholdValue = 2
	foundMasses = set()
	allSTNDclusters = []
	for distance1 in sortedConvolutios:
		distance = distance1[0]
		min_distance = finalClusters[distance][2]
		max_distance = finalClusters[distance][3]
		peak_pairs = finalClusters[distance][4]
		currentVal = ALLcleanedAutoconvolutions[distance]
		if currentVal <1: 
			continue
		groupFound = False
		if currentVal >= thresholdValue:
			for mass in set(standardAminoMasses):
				# if abs(mass-distance)<= finalClusters[distance][1] + e:
				if (mass-distance)<(e) and ((abs(mass - min_distance) < e) or (abs(mass - max_distance) < e) or (min_distance<mass and mass<max_distance) ):
					relatedClusters[mass].append(ALLcleanedAutoconvolutions[distance]) 
					standardAutconvCleaned[mass] += ALLcleanedAutoconvolutions[distance]
					stand_aa_peak_pairs[mass] += peak_pairs
					if not groupFound:
						numSTDFound+=1
						groupFound = True
					if mass not in foundMasses:
						foundMasses.add(mass)
		totalIncluded +=1

	numSTDFound = len(foundMasses)
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# COUNTING NUMBER OF CORRECT AMINO ACIDS (monomers) IN THE CONVOLUTIONS
	numCorrectConvolutions = 0
	totalstd = 0
	totalstd = 0

	for x in finalClusters:
		apeak = []
		z=  [finalClusters[x][0],finalClusters[x][1],finalClusters[x][2],finalClusters[x][3]]
		apeak = [a[0] for a in finalClusters[x][4]]
		apeak += [a[1] for a in finalClusters[x][4]]
	import math
	return numSTDFound, N , topNconvolutions ,standardAutconvCleaned, int(math.ceil(thresholdValue)), stand_aa_peak_pairs


def find_repeat_chain(clusters): 
	#finds a chain of peaks equal to the representative of the cluster
	chains = {}
	n =0 
	for pair in clusters:
		chains[n] = [sorted(pair)[0], sorted(pair)[1]]
		n+=1
	for pair in clusters:
		for i in chains:
			chain = chains[i]
			if sorted(pair)[0] == chain[-1]:
				chain.append(sorted(pair)[1])
	max_chain_size = max([len(chains[x]) for x in chains])

	return max_chain_size


def generate_convolutions(
	standardAminoMasses,peaks,e,thresh,charge,representative, precursorMass):
	constituentMonoMers = []
	finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge = generateAllAutconv(peaks,e,thresh,representative,standardAminoMasses,precursorMass,charge)
	return finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge, standardAminoMasses, e


def find_proteinogenic_clusters(
	finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge,standardAminoMasses,e,alpha,betta):
	numStnd, N , topNconvolutions ,standardAutconvCleaned, thresholdValue, stand_aa_peak_pairs= findSpectraPassingThresh(
		finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge,standardAminoMasses,e,alpha,betta)
	return numStnd, N, topNconvolutions, standardAutconvCleaned, thresholdValue,stand_aa_peak_pairs

def output_cyclopeptide_polymers(
	numStnd, N, topconvolutions, distancesCleaned,polymerMasses,charge,output, 
	precursor, retention,lines,writeOriginalSpectra,num_peaks,peptide,alpha,betta,thresholdValue,stand_aa_peak_pairs):
	
	totalCyclopeptide = 0

	compound = "nonpeptidic"
	polymer= False
	if numStnd >-1:
		sortedtop30Convolutions = sorted(distancesCleaned.items(), key=itemgetter(1), reverse=True)
		founds = [d for d in sortedtop30Convolutions if d[1]>0]

		final_STNDconvolutions = [] # this list will hold all the final convolutions that will contribute to the final convolutions with filtered convoltuions
		
		checkMasses = []
		numPolymers = []

		for mass in [x[0] for x in sortedtop30Convolutions if x[1]>0]:
			if mass in polymerMasses:
				if mass>30:
					numPolymers.append(mass)
			elif int(mass) not in set(checkMasses):
				if int(mass)==87:
					continue
				if charge==2 and (2*int(mass) in set(checkMasses) or int(mass)/2 in set(checkMasses)):
					continue
				final_STNDconvolutions.append(mass)
				# checkMasses += [int(mass)-1,int(mass),int(mass)+1, int(mass)+2, int(mass)-2]
				checkMasses += [int(mass)-1,int(mass),int(mass)+1]
		if len(numPolymers)>1:
			if len(final_STNDconvolutions)< 3*len(numPolymers)+1:
				polymer= True
		import math
		
		if len(final_STNDconvolutions)>1 and not polymer:
			totalCyclopeptide +=1
			compound = "cyclopeptide"
		elif len(final_STNDconvolutions)>1 and polymer:
			found_chains = 0
			n=0
			chains = []
			for mass in set(numPolymers):
				n+=1
				max_chain =  find_repeat_chain(stand_aa_peak_pairs[mass])
				chains.append(max_chain)
				if max_chain>3:
					found_chains +=1
			if found_chains >1:
				compound = "polymer"
			else:
				compound = "cyclopeptide"

		elif len(numPolymers)>1:
			compound = "polymer"

		else:
			compoud = "unclassified"

		return compound, str([x for x in sortedtop30Convolutions if x[1]>0]).replace(" ", ""), numStnd


