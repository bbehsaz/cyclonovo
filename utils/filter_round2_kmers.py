
def filter_intensity_kmerScore(
	numStnd, N , topNconvolutions ,standardAutconvCleaned, thresholdValue, stand_aa_peak_pairs, polymerPeaks,intensities,generate_All_Kmers,realPepMass,e):
	from operator import itemgetter
	import itertools
	polymerIntensity = 0
	final_compound = "unclassified"
	for peak in polymerPeaks:
		polymerIntensity += intensities[peak]
	allintensity = sum(intensities.values())
	usedIntensity = 0
	allPeaksDescribed= set()
	for x in stand_aa_peak_pairs:
		if standardAutconvCleaned[x]>1:
			for pair in stand_aa_peak_pairs[x]:
				allPeaksDescribed.add(pair[0])
				allPeaksDescribed.add(pair[1])		
	usedIntensity = sum([intensities[p] for p in allPeaksDescribed])
	temppercentProteinoPeaks= usedIntensity/float(allintensity)
	percentPolymerPeaks = polymerIntensity/float(allintensity)
	if temppercentProteinoPeaks > percentPolymerPeaks:
		percentProteinoPeaks = temppercentProteinoPeaks
		final_compound = 'polymer'
	else:
		percentProteinoPeaks = min(percentPolymerPeaks,temppercentProteinoPeaks)
	building_blocks = [key for key, value in standardAutconvCleaned.iteritems() if value>1]
	kmerSize = 5
	all_kmerScores, all_kmerFrequencies , all_kmerSequences, all_kmerMatches = generate_All_Kmers(
		building_blocks,kmerSize,realPepMass,intensities.keys(),e)
	highetScoredKmer = ("",0)
	sorted_kmers = sorted(all_kmerScores[5].items(), key=itemgetter(1),reverse=True)
	if len(all_kmerScores[5])>0:
		highetScoredKmer = sorted_kmers[0]
	topkmer = ""
	if highetScoredKmer[1]>3 and percentProteinoPeaks>0.6:
		final_compound = "cyclopeptide"
		topkmer = highetScoredKmer[0]
		if percentProteinoPeaks <0.7 and percentPolymerPeaks>0.35:
						final_compound="polymer"
	if topkmer != "":
		roundedKmer = "-".join([str(round(float(km),2)) for km in topkmer.split("-")])
	else:
		roundedKmer = "-"
	return final_compound, highetScoredKmer[1],percentPolymerPeaks,percentProteinoPeaks, all_kmerScores, all_kmerFrequencies , all_kmerSequences, all_kmerMatches, roundedKmer


# def reportResults(cyclopeptide_spectra,cyclopeptide_convolutions):
# 	cyclopeptide_spectra = open(output+"_cyclopeptides.mgf","a")
# 	cyclopeptide_convolutions= open(output+"_cyclopeptides_convolutions.txt","a")

# 	info = [precursorMass, peptide, percentAnnotatedPeaks, highetScoredKmer[1],charges[peptide], percentPolymerPeaks, compound_typ]
# 	# info = ["("+str(alpha)+","+str(betta)+")",pepID, pval, score, compound_type, str(len(correctSeq))]
# 	compound_type_reports_file[(alpha,betta)].write("\t".join(str(z) for z in info)+"\n")




