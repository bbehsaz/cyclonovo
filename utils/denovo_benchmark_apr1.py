#!/usr/bin/env python
# import sys
import sys, string, random
import subprocess

def output_denovo_results(candidatePeptides, precursorMass, retention, charge, peptide, benchmark_file, kmerSize, kmerThreshold,reconstructions_file):
	#This function just the final reconstructions in the format wanted! This one is for benchamrking CYCLOLIBRARY
	from operator import itemgetter
	
	reconstruction_sorted = sorted(candidatePeptides.items(), key= itemgetter(1), reverse=True)
	scores = set()
	
	peptidesUptoCorrectScore = []
	correctFound = False

	for (reconstruction,score) in reconstruction_sorted:
		# scores.add(score)
		if len(scores)>4:
			break
		elif score >12:
			peptidesUptoCorrectScore.append(reconstruction)
			scores.add(score)
		if score>11:
			reconstructions_file.write("{}\t{}\t{}\t{}\t{}\n".format(",".join([str(x) for x in reconstruction]),precursorMass, retention, charge, score))
	if len(scores) == 0:
		maxScore = 0
	else:
		maxScore = max(scores)

	# arguments_to_print = [(kmerSize,kmerThreshold), precursorMass, retention, charge, len(peptidesUptoCorrectScore), len(scores), maxScore]
	# benchmark_file.write("\t".join([str(x) for x in arguments_to_print])+"\n")

	if maxScore>12:
		return 1
	else:
		return 0


def denovo_sequence(peaks, pepMass, building_blocks,kmerSize,e,delta,kmerThreshold,verboseprint):
# De novo sequence spectrum with mass pepMass using amino acids in building blocks 
	lengths = choose_length(pepMass)
	allcompos = []
	allFeasibleCombinations, failedLenths = generate_feasible_combinations(
		building_blocks, lengths, peaks,pepMass, delta)

	allMatches = {}
	for k in range(1,kmerSize):
			currentMatches = matchAllKmers(building_blocks, k, peaks, e)
			allMatches.update(currentMatches)
	#} CHanging back
	for l in allFeasibleCombinations:
		if l < kmerSize+1:
			continue
		allcompos += list([sorted(x) for x in allFeasibleCombinations[l]])
	verboseprint("Checking feasible combinations")	
	kmerScores, kmerFrequencies, kmerSequences = {}, {}, {}
	
	candidatePeptides = {} #tuples (peptide, score)
	for l in [x for x in allFeasibleCombinations if len(allFeasibleCombinations[x])>0]:
		for combination in allFeasibleCombinations[l]:
			verboseprint("Checking amino acid combination:\t{}".format(combination))
			all_kmerFrequencies = {}
			kmers = {}
			all_kmerSequences= {}
			kmers,all_kmerFrequencies[kmerSize],all_kmerSequences[kmerSize] = generateAllKmerLinearScore(allMatches,list(combination),kmerSize, pepMass)
			chosenKmers = { key:value for key, value in kmers.items() if value> kmerThreshold}

			# verboseprint("#Putative Kmers: {}".format(len(chosenKmers)))
			adjacencies, nodes, edges = constructDeBruijn(chosenKmers, all_kmerFrequencies, all_kmerSequences,kmerSize) # Build the graph
			adjacencies, nodes, edges = pruneGraph(adjacencies,nodes,edges)
			listAllPaths = []
			for node in nodes:
				if (pepMass + protonMass) - sum(nodes[node][1])>-e:
					listAllPaths.extend(findWalksLengthK(combination, adjacencies,node,l,nodes,building_blocks))
			allPaths= set(listAllPaths)

			if len(allPaths) < 1:
				verboseprint("No feasible cycle found!")
				continue
			checkedSeqs = set()

			for path in allPaths:
				for x in combination:
					if combination.count(x)!= path.count(x):
						continue
				if not pathIsClosed(list(path),chosenKmers,kmerSize):
					continue
				checked = False

				for permutation in relatedPermutations(list(path)):
					if permutation in checkedSeqs:
						checked = True
						break
				if checked:
					continue
				else:
					checkedSeqs.add(path)
			numfeasibleCycles=0
			for path in checkedSeqs:
				matchedPeaksFragments = findCyclicScore(path,kmerSize,peaks,e)
				score =  len(set([x[0] for x in matchedPeaksFragments]))
				candidatePeptides[tuple(path)] = score
				if score>11:
					# verboseprint("{}\t{}".format(path, score))
					numfeasibleCycles +=1
			if numfeasibleCycles>0:
				verboseprint("Number of feasible cycles: {}".format(numfeasibleCycles))
			else:
				verboseprint("No feasible cycle found!")
	return candidatePeptides
					


def matchAllKmers(foundAutoConv, k, peaks, e): 
#this function takes the list of selected autoconvolutions (FoundAutoConv) and returns 
#a dictionary containing all possible k-mers and if they match (with error e) a peak (1 if yes, 0 if no)
	import itertools
	allMatches = {}
	protonMass = 1.00728
	for x in itertools.combinations_with_replacement(foundAutoConv, k):
		elements = x
		for kmer in itertools.permutations(elements):
			if kmer in allMatches:
				continue
			allMatches[kmer] = 0
			for peak in peaks:
				if abs(sum(kmer) - (peak - protonMass) ) < e:
					allMatches[kmer] = 1
					break
	return allMatches

def generateAllKmerLinearScore(allMatches,combination,k,precursorMass): 
#this function generates all the kmers using the amino acids in combination and calculaets their linear 
#score using allMatches 
	import itertools

	aaFreq = {}
	aminoFrequencies = {}
	kmerFrequencies = {}
	kmerSequences = {}
	allKmersScore = {}
	list_combination = list(combination)
	aminos_set = set(list_combination)
	for x in aminos_set:
		aminoFrequencies[x] = list_combination.count(x)
	for compos in itertools.combinations_with_replacement(combination, k):
		for x in itertools.permutations(list(compos), k):
			# kmerName = x
			kmer = list(x)
			kmerName = tuple(kmer)
			if kmerName in allKmersScore:
				continue
			notPossible = False			
			kmerFrequencies[kmerName] = {a:kmer.count(a) for a in aminoFrequencies.keys() }
			kmerSequences[kmerName] = kmer
			for a in set(kmer):
				if aminoFrequencies[a] < kmerFrequencies[kmerName][a]:
					notPossible = True
					break
			if notPossible:
				# allKmersScore[kmerName] = -1 # kmer 
				continue
			if sum(kmer) > precursorMass:
				# allKmersScore[kmerName] = -1 
				continue
			

			allKmersScore[kmerName] = 0
			for i in range(k-1):
				for j in range(i+1,k):
					linearFragment = kmer[i:j]
					# print j
					# print linearFragment
					# if len(linearFragment) = 1:
					# 	linearFragment = 
					allKmersScore[kmerName] += allMatches[tuple(linearFragment)]
	return allKmersScore,kmerFrequencies, kmerSequences
	

def linear_score_fragments(fragments, peaks, e):
	#matchedFragments = []
	score = 0
	protonMass = 1.00728
	matched = []
	for fragment in fragments:
		for peak in sorted(peaks):
			if abs(peak-(sum(fragment)+protonMass))<e:
				score += 1 
				matched.append((peak,round(sum(fragment),2)))
				# matched[peak] = tuple(fragment)
				break
	return matched
def findCyclicScore(peptide,kmerSize,peaks,e):
	cThSpect = {}
	listCThSpect = []
	seThSpect = {}
	pepname = peptide
	n=len(pepname)
	revpeptide = list(pepname)
	revpeptide.reverse()
	matches = []
	allfragments = []
	for i in range(len(pepname)):
		linearVersion = peptide[i:] + peptide[:i]	
		# matches+=kmerMatches["-".join(str(x) for x in peptide[0:kmerSize])]		
		# matches += linear_score_fragments(
		allfragments +=	[linearVersion[0:j] for j in range(len(pepname))]
	allfragments += [tuple(peptide)]
	matches = linear_score_fragments(list(set(allfragments)), peaks, e  ) 
	# print (len(matches))
	return matches

def generate_kmers_score(combination,kmerSize,pepMass,peaks,e): 
#this function generates all the kmers using the amino acids in combination and calculaets their linear 
#score using an iterative method 
	aaFreq = {}
	aminoFrequencies = {}

	list_combination = list(combination)
	aminos_set = set(list_combination)
	for x in aminos_set:
		aminoFrequencies[x] = list_combination.count(x)
	def extend_to_k(k,kmerScores, kmerFrequencies, kmerSequences): #all-kmers with their linear scores using the k-1-mers in the kmers_scores
		kmerScores[k] = {}
		kmerFrequencies[k] = {} #Frequency of all possible aa's in kmers with size k
		kmerSequences[k] = {}

		for aa in aminos_set:
			for k1mer in kmerScores[k-1]:
				k1merSeq = kmerSequences[k-1][k1mer][:]
				newKmer_sequence = k1merSeq + [aa]
				if len(newKmer_sequence)>0:
					if sum(newKmer_sequence) >pepMass:
						continue
				if kmerFrequencies[k-1][k1mer][aa]+1 > aminoFrequencies[aa]:
					continue
				if k1mer == '':
					newKmer = str(aa)
				else:
					newKmer = k1mer + "-" + str(aa)
				newScore =kmerScores[k-1][k1mer][:]
				newScore += linear_score_fragments(
					[newKmer_sequence[i:] for i in range(1,k)], peaks, e)
				kmerScores[k][newKmer] = newScore
				kmerFrequencies[k][newKmer] = {}
				for x in kmerFrequencies[k-1][k1mer]:
					kmerFrequencies[k][newKmer][x] = kmerFrequencies[k-1][k1mer][x]
				kmerFrequencies[k][newKmer][aa] = kmerFrequencies[k-1][k1mer][aa] + 1
				kmerSequences[k][newKmer] = newKmer_sequence
		newk = k+1

		if newk== kmerSize + 1:
			return kmerScores, kmerFrequencies , kmerSequences
		return extend_to_k(newk, kmerScores, kmerFrequencies, kmerSequences)

	kmerScores = {}
	kmerFrequencies = {}
	kmerSequences = {}
	kmerScores[0] = {'':[]}
	kmerFrequencies[0] = {'':{aa:0 for aa in aminos_set}}
	kmerSequences[0] = {'':[]}
	all_kmerScores, all_kmerFrequencies , all_kmerSequences  = extend_to_k(
		1, kmerScores, kmerFrequencies, kmerSequences)
	final_allKmerScores = {}
	kmerMatches = all_kmerScores[kmerSize].copy()

	for size in all_kmerScores:
		final_allKmerScores[size] = {}
		for kmer in all_kmerScores[size]:
			final_allKmerScores[size][kmer] = len(( [x for x in all_kmerScores[size][kmer]] ))
	return final_allKmerScores, all_kmerFrequencies , all_kmerSequences, kmerMatches



def generate_All_Kmers(combination,kmerSize,pepMass,peaks,e): 
#this function generates all the kmers using the amino acids in combination and calculaets their linear 
#score using an iterative method 
	aaFreq = {}
	aminoFrequencies = {}

	list_combination = list(combination)
	aminos_set = set(list_combination)
	for x in aminos_set:
		aminoFrequencies[x] = list_combination.count(x)
	def extend_to_k(k,kmerScores, kmerFrequencies, kmerSequences): #all-kmers with their linear scores using the k-1-mers in the kmers_scores
		kmerScores[k] = {}
		kmerFrequencies[k] = {} #Frequency of all possible aa's in kmers with size k
		kmerSequences[k] = {}

		for aa in aminos_set:
			for k1mer in kmerScores[k-1]:
				k1merSeq = kmerSequences[k-1][k1mer][:]
				newKmer_sequence = k1merSeq + [aa]
				if len(newKmer_sequence)>0:
					if sum(newKmer_sequence) >pepMass:
						continue
				if k1mer == '':
					newKmer = str(aa)
				else:
					newKmer = k1mer + "-" + str(aa)
				newScore =kmerScores[k-1][k1mer][:]
				newScore += linear_score_fragments(
					[newKmer_sequence[i:] for i in range(1,k)], peaks, e)
				kmerScores[k][newKmer] = newScore
				kmerFrequencies[k][newKmer] = {}
				for x in kmerFrequencies[k-1][k1mer]:
					kmerFrequencies[k][newKmer][x] = kmerFrequencies[k-1][k1mer][x]
				kmerFrequencies[k][newKmer][aa] = kmerFrequencies[k-1][k1mer][aa] + 1
				kmerSequences[k][newKmer] = newKmer_sequence
		newk = k+1

		if newk== kmerSize + 1:
			return kmerScores, kmerFrequencies , kmerSequences
		return extend_to_k(newk, kmerScores, kmerFrequencies, kmerSequences)

	kmerScores = {}
	kmerFrequencies = {}
	kmerSequences = {}
	kmerScores[0] = {'':[]}
	kmerFrequencies[0] = {'':{aa:0 for aa in aminos_set}}
	kmerSequences[0] = {'':[]}
	all_kmerScores, all_kmerFrequencies , all_kmerSequences  = extend_to_k(
		1, kmerScores, kmerFrequencies, kmerSequences)
	final_allKmerScores = {}
	kmerMatches = all_kmerScores[kmerSize].copy()

	for size in all_kmerScores:
		final_allKmerScores[size] = {}
		for kmer in all_kmerScores[size]:
			final_allKmerScores[size][kmer] = len(( [x for x in all_kmerScores[size][kmer]] ))
	return final_allKmerScores, all_kmerFrequencies , all_kmerSequences, kmerMatches



def constructDeBruijn(kmers, all_kmerFrequencies, all_kmerSequences,kmerSize):
# Constructs de Bruijn graph from k-mers with k-1
	adjacencies = {}
	nodes2aminos = {}
	n = 0
	edges = []
	nodes = {}
	def suffix(string):
		return string[1:]
	def prefix(string):
		return string[:len(string)-1]
	#To get the deburijn graph that actually nodes are the k-mer size ... only for plotting purposes!
	for kmer in kmers:
		# add prefix and suffix nodes
		prefixK1merSeq = prefix(all_kmerSequences[kmerSize][kmer])
		prefixK1merName = "-".join(str(aa) for aa in prefixK1merSeq)
		# prefixK1merFrequencies = all_kmerFrequencies[kmerSize-1][prefixK1merName]
		prefixK1merFrequencies = {a:prefixK1merSeq.count(a) for a in set(prefixK1merSeq)}
		nodes[prefixK1merName] = (prefixK1merSeq,prefixK1merFrequencies)
		
		suffixK1merSeq = suffix(all_kmerSequences[kmerSize][kmer])
		suffixK1merName = "-".join(str(aa) for aa in suffixK1merSeq)
		# suffixK1merFrequencies = all_kmerFrequencies[kmerSize-1][suffixK1merName]
		suffixK1merFrequencies = {a:suffixK1merSeq.count(a) for a in set(suffixK1merSeq)}
		nodes[suffixK1merName] = (suffixK1merSeq,suffixK1merFrequencies)

		if prefixK1merName in adjacencies:
			adjacencies[prefixK1merName].append(suffixK1merName)
			edges.append((prefixK1merName,suffixK1merName))
		else:
			adjacencies[prefixK1merName] = [suffixK1merName]
			edges.append((prefixK1merName,suffixK1merName))

	return adjacencies, nodes, edges




def pathIsClosed(path, chosenKmers, kmerSize):
	closed = True
	for i in range(0,kmerSize):
		kmerSeq = path[-(kmerSize-i):]+path[:i]
		kmerName = "-".join([str(x) for x in kmerSeq])
		kmerName=tuple(kmerSeq)
		if kmerName in chosenKmers:
			continue
		else:
			closed = False
			return closed
	return closed



def relatedPermutations(permutation):
	cyclicPermutations = [permutation[i:]+permutation[:i] for i in range(len(permutation))]
	permutation.reverse()
	reverseCyclicPermutations = [permutation[i:]+permutation[:i] for i in range(len(permutation))]
	allRelated = cyclicPermutations+reverseCyclicPermutations
	allRelatedNames=[]
	for x in allRelated:
		name = tuple([y for y in x])
		# name = '-'.join([str(y) for y in x])		
		allRelatedNames.append(name)
		x.reverse()
		name = tuple([y for y in x])
		# name = '-'.join([str(y) for y in x])
		allRelatedNames.append(name)

	return allRelatedNames

protonMass = 1.00728


def findWalksLengthK(combination, graph, start, l, nodes,building_blocks):
#Finds all Walks of length k from the node start in a directed graph. Uses modified BFS for traversal.
	combinationFreq = {aa:combination.count(aa) for aa in set(building_blocks)}
	allPaths = []
	pathAAfrequencies = {}
	for aa in building_blocks:
		if aa in nodes[start][1]:
			pathAAfrequencies[aa] = nodes[start][1][aa]
		else:
			pathAAfrequencies[aa] = 0
	startSequence = nodes[start][0]
	queue = [(start,pathAAfrequencies, startSequence)]
	while True:
		if not queue:
			break
		(vertex,frequencies, path) = queue.pop(0)
		if vertex in graph:
			for nextName in set(graph[vertex]):
				next = nodes[nextName][0][:]
				nextAA = next[-1]
				newFrequencies = frequencies.copy()
				if nextAA in frequencies:
					newFrequencies[nextAA]= frequencies[nextAA]+1
				else:
					newFrequencies[nextAA]= 1
				if newFrequencies[nextAA] > combinationFreq[nextAA]:
					continue
				newPath = path + [nextAA]
				if len(newPath) == l:
					feasiblePath = True
					for aa in building_blocks:
						if newFrequencies[aa] != combinationFreq[aa]:
							feasiblePath = False
							break
					if feasiblePath:
						allPaths.append(newPath)
				if len(newPath) < l:
					queue.append((nextName,newFrequencies,newPath))
				if len(newPath) > l+1:
					break
	return [tuple(path1) for path1 in allPaths]

def spellPath(path):
	spelledSeq = list(path[0])
	for i in range(1,len(path)):
		# pointer = len(kmers[0])-1
		spelledSeq = spelledSeq + [path[i][-1]]
	return spelledSeq
def findLeaves(adjacencies1): 
	newleaves= []
	for node in adjacencies1:
		if len(adjacencies1[node]) == 0:
			newleaves.append(node)
	return newleaves
def pruneGraph(adjacencies,nodes,edges):
	backwardEdges = {}
	leaves = {}
	newleaves= []
	singletons = []
	for node in nodes:
		backwardEdges[node] =[]
	for node in nodes:
		if node not in adjacencies:
			adjacencies[node] = []
			continue
		for tail in adjacencies[node]:
			
			backwardEdges[tail].append(node)

	for node in nodes:
		if len(adjacencies[node])==0:
			if node not in backwardEdges:
				singletons.append(node)
	[nodes.remove(x) for x in singletons]
	for node in nodes:
		if node not in adjacencies:
			adjacencies[node] = []
	# singletons = set(nodes) - set(set(nodes) - set(adjacencies.keys()) )

	newleaves= findLeaves(adjacencies)
	while len(newleaves)>0:

		leaves = newleaves[:]
		newleaves = []
		
		for leaf in leaves:

			for parent in backwardEdges[leaf]:
			# parent = backwardEdges[leaf][0]
			# if len(adjacencies[parent]) ==1:
			# 	newleaves.append(parent)
				adjacencies[parent].remove(leaf)
				edges.remove((parent,leaf))				
			# backwardEdges[leaf].remove(parent)
			del backwardEdges[leaf]
			del adjacencies[leaf]
			del nodes[leaf]
		newleaves = findLeaves(adjacencies)

	return adjacencies, nodes , edges

def generate_feasible_combinations(foundAutoConv, pepLengths, peaks, pepMass, delta):
# Generates all feasible combinations for lengths that 
	failedLenths  = []
	allFeasibleCombinations = {}
	for l in pepLengths:
		allFeasibleCombinations[l] = set()		
		allFeasibleCombinations[l]= find_feasible_combintions(
			foundAutoConv, l, pepMass, delta)
		if len(allFeasibleCombinations) == 0:
			failedLenths.append(l)
	return allFeasibleCombinations, failedLenths

def find_feasible_combintions(foundAutoConv, l, peptideMass, delta): 
# Among all possible combinations of size l of the selected autconvolutions (foundAutoConv) returns the ones that
# their sum matches the precursor mass allowing for delta error
	listfeasibleCombinations = set()
	from operator import itemgetter
	import itertools
	protonMass = 1.00728
	for combination in itertools.combinations_with_replacement(foundAutoConv, l):
		if abs(sum(combination) - ( peptideMass) ) < delta:
			listfeasibleCombinations.add(tuple(sorted([float(x) for x in combination])))
	return listfeasibleCombinations

def choose_length(pepMass):
	if pepMass > 1100:
		pepLengths = [9,10,11]
	if pepMass > 950:
		pepLengths = [8,9,10]
	elif pepMass > 800:
		pepLengths = [7,8,9]
	elif pepMass > 700:
		pepLengths = [6,7,8,9]
	elif pepMass > 600:
		pepLengths = [5,6,7]
	elif pepMass > 500:
		pepLengths = [4,5,6,7]
	elif pepMass > 400:
		pepLengths = [3,4,5]
	return pepLengths



def find_aa_for_denovo(
	finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge,standardAminoMasses,e,polymer_repeat_units):
# find proteinogenic clusters with convolution 2 or higher as building blocks for de novo sequencing 
	standardAutconvPeaks = {}
	standardAutconvCleaned = {}
	from operator import itemgetter
	stand_aa_peak_pairs = {}
	for mass in standardAminoMasses:
		standardAutconvCleaned[mass] = 0
		# peaksIndexAminoDistances[mass] = []
		standardAutconvPeaks[mass] = []
		stand_aa_peak_pairs[mass] = []
	N=30
	polymerPeaks = {}
	for mass in polymer_repeat_units:
		polymerPeaks[mass] = []
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
	protonMass = 1.00728
	pepMass = precursorMass*charge - protonMass
	import math
	thresholdValue = 1
	foundMasses = set()
	allSTNDclusters = []
	polymerPeaks = set()
	for distance1 in sortedConvolutios:
		distance = distance1[0]
		min_distance = finalClusters[distance][2]
		max_distance = finalClusters[distance][3]
		peak_pairs = finalClusters[distance][4]
		currentVal = ALLcleanedAutoconvolutions[distance]
		if currentVal <1: 
			continue
		groupFound = False
		if currentVal > thresholdValue:
			for mass in set(standardAminoMasses):
				if mass in polymer_repeat_units:
					if abs(mass-distance)<e and ((abs(mass - min_distance) < e) or (abs(mass - max_distance) < e) or (min_distance<mass and mass<max_distance) ):
						for pair in peak_pairs:
							polymerPeaks.add(pair[0])
							polymerPeaks.add(pair[1])
					continue

				if abs(mass-distance)<=2*e and ((abs(mass - min_distance) < e) or (abs(mass - max_distance) < e) or (min_distance<mass and mass<max_distance) ):
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
	import math
	return numSTDFound, N , topNconvolutions , standardAutconvCleaned, int(math.ceil(thresholdValue)), stand_aa_peak_pairs,polymerPeaks





