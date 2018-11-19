#!/usr/bin/env python
#################################################
# CycloNovo: algorithm for high-throughput de novo cyclopeptide sequencing using tandem mass spectrometry
# Pevzner lab 
# University of California San Diego, La Jolla, CA, USA
#################################################
import sys
import os
import traceback
import logging
from os.path import basename

from os.path import join, abspath, isdir
from datetime import datetime
from site import addsitedir
from optparse import OptionParser
import argparse
import os
from utils.filter_apr1 import *
from utils.denovo_benchmark_apr1 import *
from utils.filter_round2_kmers import *
from scripts.readMGF import readMGF
from scripts.readStream import readStream
from scripts.readStringMGF import readStringMGF

from utils.sharedFunctions import getSpecializedBuildingBlocks
from utils.sharedFunctions import writeOriginalSpectra
import math
from operator import itemgetter
import itertools
from subprocess import Popen,PIPE

def is_valid_file(parser, arg):
	"""
	Check if arg is a valid file that already exists on the file system.
	Parameters
	"""
	arg = os.path.abspath(arg)
	if not os.path.exists(arg):
		parser.error("The file %s does not exist!" % arg)
	else:
		return arg


def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc:  # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass
		else:
			raise

def get_original_fpath(fpath, file_mapping=None):
	return file_mapping[basename(fpath)] if file_mapping and basename(fpath) in file_mapping else fpath


total = 0
def get_parser():
	"""Parse arguments and check if the spectrum file exists. """
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(prog="cyclonovo.py",description= "Algorith for analyzing cyclopeptides using Tandem Mass Spectra ...", 
		formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("-s", "--spectrum", dest="mgfFile", type=lambda x: is_valid_file(parser, x), help="Input spectra (mgf/mzXML)", metavar="FILE", required=True)
	parser.add_argument("-o", "--output", dest="output", help="Output dir",required=True)
	parser.add_argument('--denovo', dest='denovo', help="De novo sequence cyclopeptide spectra", action='store_true')
	# parser.add_argument("-d", "--precursor_ion_thresh", dest='delta', help='Precursor Ion Thresh for custom Running mode', default='0.02')
	parser.add_argument("-e","--ion_thresh", dest='e', help='Fragment/Precursor Ion Thresh for custom Running mode', default='0.02')
	parser.add_argument("--pname", dest="pname",
						help="Prefix for generated files", default="cyclonovo")
	parser.add_argument("--monomers", dest="monomers", action="store",
						help="Monomers to be considered ('standard','top25' or'ripps')",
						choices=['standard', 'top25', 'ripps'], default='standard')

	# parser.add_argument("--preprocess", dest="preprocess", action="store_true",
	# 					help='Filter and Preprocess the spectra (otherwise it\'s assumed it\'s preprocessed\')')
	# parser.add_argument("--reprsentative", dest='represent', help='Median or Average as cluster representatives', default='median')
	parser.add_argument("-k", "--kmer_size", dest='kmerSize', help='k-mer size for building de Bruijn graph', default=int(2),type=int)
	parser.add_argument("-t", "--kmer_threshold", dest='kmerThreshold', help='Threshold value for selecting k-mers for graph', default=int(1),type=int)
	parser.add_argument("-b", "--betta", dest='betta', help='Betta value for calculating multiplicity threshold', default=float(-1),type=float)
	parser.add_argument("-a", "--alpha", dest='alpha', help='Alpha value for calculating multiplicity thresholds', default=float(150.0),type=float)
	parser.add_argument("-v","--verbosity", dest="verbose", action="store_true")
	# parser.add_argument("--topcyclic", dest="topcyclic", action="store_true",help="Consider most common cyclopeptidic amino acids", default=False)
	parser.add_argument("--preprocess", dest="filter", action="store_true",help='Filter and Preprocess the spectra (otherwise it\'s assumed it\'s preprocessed\')')
	# parser.add_argument("--ripps", dest="ripps", action="store_true",help="Consider modifications common in RiPPs", default=False)
	parser.add_argument("--kmer_score", dest='kmerScoreIdent', help='k-mer score threshold for identifying cyclopeptidic spectra', default=int(4),type=int)
	parser.add_argument("--cyclointensity", dest='cyclointensityThresh', help='cycloIntesity threshold for identifying cyclopeptidic spectra', default=int(60),type=int)
	parser.add_argument("--num_frequent_clusters", dest='cyclointensityThresh', help='cycloIntesity threshold for identifying cyclopeptidic spectra', default=int(2),type=int)
	parser.add_argument("--aminoacid_multiplicity", dest='aminoThresh', help='Amino acid multiplicity threshold', default=int(1),type=int)
	return parser

if __name__ == "__main__":
	args = get_parser().parse_args()
	cyclonovoPath = os.path.dirname(os.path.realpath(__file__))
	
	if args.verbose:
		def verboseprint(stringtoprint):
			print(stringtoprint)
	else:
		verboseprint = lambda *a, **k: None # do-nothing function
	verboseprint("++ Arguments: " + str(args))
	verboseprint("++ Running CycloNovo from directory: " + cyclonovoPath)
	verboseprint("++ Reading input spectrum: " + args.mgfFile)
	protonMass = 1.00728 #This is simply mass of hydrogen --- do not remove!
	digits = 2
	# building_blocks_main = getBuildingBlocks(cyclonovoPath+"/configs/aminoacidMasses.txt")
	# if args.ripps:
	# 	building_blocks_main = getBuildingBlocks(cyclonovoPath+"/configs/aa_polymer_cyclomasses_ripps.txt")
	# if args.topcyclic:
	# 	building_blocks_main = getBuildingBlocks(cyclonovoPath+"/configs/aa_polymer_masses_25.txt")
	# else:	
	# 	building_blocks_main = getBuildingBlocks(cyclonovoPath+"/configs/aa_polymer_masses.txt")
	# print building_blocks_main
	# polymer_repeat_units = getBuildingBlocks(cyclonovoPath+"/configs/polymer_repeat_masses.txt")
	building_blocks_main, polymer_repeat_units = getSpecializedBuildingBlocks(args.monomers, cyclonovoPath)
	outputdir = args.output
	if args.filter:
		Process=Popen([cyclonovoPath+'/scripts/print_spectrum '+ str(args.mgfFile) +' --scan_num -1 --print_spectrum '],shell=True,stdin=PIPE, stdout=PIPE)
	else:
		Process=Popen([cyclonovoPath+'/scripts/print_spectrum '+ str(args.mgfFile) +' --scan_num -1 --print_spectrum --no_filter --no_merge'],shell=True,stdin=PIPE, stdout=PIPE)
		# Process=Popen(['cat '+ str(args.mgfFile) ],shell=True,stdin=PIPE, stdout=PIPE)
	peaksfile = Process.communicate()[0]
	mkdir_p(outputdir)
	output = os.path.join(outputdir, args.pname)
	peaksnIntensity,pepMasses,charges,retentions,fileLines= readStringMGF(peaksfile)
	

	verboseprint("# of spectra in the file:\t{}".format(len(peaksnIntensity)))
	# peaksnIntensity,pepMasses,charges,retentions,fileLines= readStream(streamFile)
	representative = 'median'
	peptides_considered = 0
	num_raw_spectra = 0
	num_spectra_analyzed = 0
	num_cyclopeptide_spectra_dic = {}
	num_polymer_spectra_dic = {}
	compound_type_reports_file = {}
	betta_values = [args.betta]
	alpha_values = [args.alpha]

	for x in alpha_values:
			alpha = round(1/float(x),3)
			for betta in betta_values:
				num_cyclopeptide_spectra_dic[(alpha,betta)] = 0
				num_polymer_spectra_dic[(alpha,betta)] = 0
				# compound_type_reports_file[(alpha,betta)] = open(output+"_compoundType_"+str(x)+"_"+str(betta)+".txt","w")
	# summary_file = open(output+"_cyclonovo_summary.txt","a")
	kValues = [int(args.kmerSize)]
	benchmark_file = {}
	reconstructions_file = {}
	kthresholdValues = [int(args.kmerThreshold)]
	for kmerSize in kValues: 
		for kmerThreshold in kthresholdValues:
				# benchmark_file[(kmerSize,kmerThreshold)] = open(output+"_sequencing_summary.txt","w")
				benchmark_file[(kmerSize,kmerThreshold)] = 1
				reconstructions_file[(kmerSize,kmerThreshold)] = open(output+"_sequencing_reconstructions.txt","w")
	allidsconsidered = []
	allsequenced = 0
	for peptide in peaksnIntensity:
		verboseprint("=======================================================")
		verboseprint("{}\t{}\t{}".format(pepMasses[peptide], retentions[peptide], charges[peptide]))
		total +=1
		addAminos = []
		thisbuildingblock = building_blocks_main.copy()
		standardAutconvCleaned = []
		precursorMass = pepMasses[peptide]
		num_raw_spectra +=1
		realPepMass = pepMasses[peptide]*charges[peptide] - protonMass
		if realPepMass < 500 or realPepMass > 2000:
			verboseprint("Peptide mass out of the acceptable range")
			continue
		if len(peaksnIntensity[peptide])<20:
			verboseprint("There's not enough peaks in MS/MS spectrum")
			continue
		num_spectra_analyzed +=1
		e= float(args.e)
		finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge, standardAminoMasses, e = generate_convolutions(
					thisbuildingblock,peaksnIntensity[peptide],float(args.e),int(precursorMass/10.0)+10,charges[peptide],representative,pepMasses[peptide])
		
		for x in alpha_values:
					alpha = round(1/float(x),3)
					for betta in betta_values:						
						suffix = "_"+str(x)+"_"+str(betta)
						num_protenogenic,N, foundConvolutions, distancesCleaned, thresholdValue, stand_aa_peak_pairs = find_proteinogenic_clusters(
							finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge,standardAminoMasses,e,alpha,betta)						
						compound_type,convolutions, numStnd = output_cyclopeptide_polymers(
							num_protenogenic,N, foundConvolutions, distancesCleaned,polymer_repeat_units,charge,output+suffix, 
							pepMasses[peptide], retentions[peptide],fileLines[peptide],writeOriginalSpectra,len(peaksnIntensity[peptide]),peptide, alpha, betta,int(thresholdValue),stand_aa_peak_pairs)
		final_compound = "unclassified"

		if compound_type == 'polymer':
			final_compound = 'notcyclic'
		# if compound_type == "cyclopeptide":
		if True:
			intensities = peaksnIntensity[peptide]
			numStnd, N , topNconvolutions ,standardAutconvCleaned, thresholdValue, stand_aa_peak_pairs, polymerPeaks= find_aa_for_denovo(
						finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge,standardAminoMasses,e, polymer_repeat_units)
			final_compound, highetScoredKmer,percentPolymerPeaks,percentProteinoPeaks, all_kmerScores, all_kmerFrequencies , all_kmerSequences, all_kmerMatches, topkmer = filter_intensity_kmerScore(
						numStnd, N , topNconvolutions ,standardAutconvCleaned, thresholdValue, stand_aa_peak_pairs, polymerPeaks,intensities,generate_All_Kmers,realPepMass,e)

		else: 
			percentPolymerPeaks = 0
			percentProteinoPeaks = 0
			highetScoredKmer = 0
		verboseprint("id\tkmerscore\tcyclointensity\t{}\t{}\t{}".format(peptide,highetScoredKmer,percentProteinoPeaks))
		if final_compound == "cyclopeptide":
			num_cyclopeptide_spectra_dic[(alpha,betta)] += 1
			cyclopeptide_spectra_file = open(output+"_cyclospectra.mgf","a")
			writeOriginalSpectra(fileLines[peptide],cyclopeptide_spectra_file)
		if final_compound == 'polymer':
			num_polymer_spectra_dic[(alpha,betta)] += 1

		if len(standardAutconvCleaned)>2:
			building_blocks = [key for key, value in standardAutconvCleaned.iteritems() if value>1]
			sorted_building_blockes = [x for x in standardAutconvCleaned.iteritems() ]
		# info = [precursorMass, retentions[peptide],charges[peptide], peptide, final_compound, percentProteinoPeaks, highetScoredKmer, percentPolymerPeaks]
		if final_compound != "cyclopeptide":
			final_compound = "notcyclic"
		building_blocks_multiplicity =  "["+ ",".join([ "("+str(round(key,2))+","+str(value)+")" for key, value in sorted(standardAutconvCleaned.items(), key=itemgetter(1), reverse=True) if value>2]) +"]"
		info = [precursorMass, retentions[peptide],charges[peptide], peptide, final_compound, highetScoredKmer, 
				percentProteinoPeaks, topkmer, building_blocks_multiplicity ] #reports classification measure for each spectrum
		compound_type_reports_file = open(output+"_classification_report.txt","a")
		compound_type_reports_file.write("\t".join(str(z) for z in info)+"\n")

		if final_compound == "cyclopeptide":
			verboseprint("classified: " + str("cyclospectrum"))
			if not args.denovo:
				continue
			for kmerSize in kValues: 
				
				for kmerThreshold in kthresholdValues: 
					verboseprint("@@@")
					verboseprint("De novo sequence ... \nkmer: size {}, threshold {}".format(kmerSize,kmerThreshold))
					if len(building_blocks)<3:
						verboseprint("Number of predicted cyclopeptidic amino acids is low")
						break
					verboseprint("Predicted Amino Acids:")
					building_blocks = []
					for p in sorted(standardAutconvCleaned.iteritems(), key=itemgetter(1), reverse=True):
						if p[1]>args.aminoThresh:
							verboseprint("{}\t{}".format(round(p[0],2),p[1]))
							building_blocks.append(p[0])
					candidateSequences = denovo_sequence(
						peaksnIntensity[peptide].keys(), realPepMass, building_blocks,kmerSize,e,e,kmerThreshold,verboseprint)
					allsequenced += output_denovo_results(candidateSequences, precursorMass, retentions[peptide], charges[peptide], peptide,
						benchmark_file[(kmerSize,kmerThreshold)], kmerSize, kmerThreshold,reconstructions_file[(kmerSize,kmerThreshold)])
					nameofrecontfile = output+"_sequencing_reconstructions.txt"
					nameofcyclospecfile = output+"_cyclospectra.mgf"
	reconstructions_file[(kmerSize,kmerThreshold)].close()
	verboseprint("calculating P-values and cleaning up ...")
	if args.denovo:
		Popen(["sh " +cyclonovoPath+'/scripts/generate_pvalus_recontfile_mgf.sh '+nameofrecontfile+ " " + nameofcyclospecfile + " " + str(output) + " "+ str(cyclonovoPath) +" " +str(e)]
					,shell=True)
	
	cyclopeptide_spectra_file.close()
	#report the number of cyclopeptides identified in output_summary.txt
	finalMGF = args.mgfFile.replace(" ", "_")
	argstoprint = [finalMGF, str(len(peaksnIntensity)), num_spectra_analyzed, num_cyclopeptide_spectra_dic[(alpha,betta)], allsequenced]
	# for (alpha,betta) in num_cyclopeptide_spectra_dic:
	# 	summary_file.write("\t".join([str(ap) for ap in argstoprint])+"\n")
	verboseprint("CycloNovo succesfully finished. Find the results in "+str(output))






