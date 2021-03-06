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
import errno
import os
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
from scripts.identiy_cyclospectra import *
from scripts.denovo_seq_spectrum import *
from utils.sharedFunctions import getSpecializedBuildingBlocks
from utils.sharedFunctions import writeOriginalSpectra
from utils.commons import *
import math
from operator import itemgetter
import itertools
from subprocess import Popen,PIPE

total = 0
def get_parser():
    """Parse arguments and check if the spectrum file exists. """
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(prog="cyclonovo.py",
        description= "Algorith for analyzing cyclopeptides using Tandem Mass Spectra ...", 
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--spectrum", dest="mgfFile", type=lambda x: is_valid_file(parser, x), help="Input spectra (mgf/mzXML)", metavar="FILE", required=True)
    parser.add_argument("-o", "--output", dest="output", help="Output dir",required=True)
    parser.add_argument('--denovo', dest='denovo', help="De novo sequence cyclopeptide spectra", action='store_true')
    # parser.add_argument("-d", "--precursor_ion_thresh", dest='delta', help='Precursor Ion Thresh for custom Running mode', default='0.02')
    parser.add_argument("-e","--ion_thresh", dest='e', help='Fragment/Precursor Ion Thresh for custom Running mode', default='0.02')
    parser.add_argument("-p", "--params", dest="params", help="params.xml with all settings (GNPS running mode)")

    parser.add_argument("--pname", dest="pname",
                        help="Prefix for generated files", default="cyclonovo")
    parser.add_argument("--monomers", dest="monomers", action="store",
                        help="Monomers to be considered ('standard','top25' or'ripps')",
                        choices=['standard', 'top25', 'ripps'], default='standard')

    # parser.add_argument("--preprocess", dest="preprocess", action="store_true",
    #                   help='Filter and Preprocess the spectra (otherwise it\'s assumed it\'s preprocessed\')')
    # parser.add_argument("--reprsentative", dest='represent', help='Median or Average as cluster representatives', default='median')
    parser.add_argument("-k", "--kmer_size", dest='kmerSize', help='k-mer size for building de Bruijn graph', default=int(5),type=int)
    parser.add_argument("-t", "--kmer_threshold", dest='kmerThreshold', help='Threshold value for selecting k-mers for graph', default=int(2),type=int)
    parser.add_argument("-b", "--beta", dest='beta', help='Betta value for calculating multiplicity threshold', default=float(-1),type=float)
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
    import xml.etree.ElementTree
    from distutils import dir_util
    import os
    from os.path import join, isfile, isdir, basename
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
def get_params_and_mapping(opts):
    print opts
    if opts.params is not None:
        verify_file(opts.params, description="params file")
        params, file_mapping = parse_params_xml(opts.params)
    else:
        params = opts.__dict__  # TODO: take a subset of all options
        file_mapping = None
    return params, file_mapping

# def initiations():
#     global representative, peptides_considered, 
    


def get_original_fpath(fpath, file_mapping=None):
    return file_mapping[basename(fpath)] if file_mapping and basename(fpath) in file_mapping else fpath


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
    params, file_mapping = get_params_and_mapping(args)

    building_blocks_main, polymer_repeat_units = getSpecializedBuildingBlocks(params["monomers"], cyclonovoPath)
    thisbuildingblock = building_blocks_main.copy()
    if args.filter:
        Process=Popen(['print_spectrum '+ str(args.mgfFile) +' --scan_num -1 --print_spectrum '],shell=True,stdin=PIPE, stdout=PIPE)
    else:
        Process=Popen(['print_spectrum '+ str(args.mgfFile) +' --scan_num -1 --print_spectrum --no_filter --no_merge'],shell=True,stdin=PIPE, stdout=PIPE)

    peaksfile = Process.communicate()[0]
    peaksnIntensity,pepMasses,charges,retentions,fileLines= readStringMGF(peaksfile)
    outputdir = params['output']
    outputprefix = os.path.join(outputdir, args.pname)
    mkdir_p(outputdir)
    output = os.path.join(outputdir, args.pname)
    verboseprint("# of spectra in the file:\t{}".format(len(peaksnIntensity)))

    representative, peptides_considered, num_raw_spectra, num_spectra_analyzed, allsequenced, e = 'median', 0, 0, 0, 0, float(args.e)
    num_cyclopeptide_spectra_dic, num_polymer_spectra_dic, compound_type_reports_file, benchmark_file, reconstructions_file= {},{},{},{},{}
    beta_values, alpha_values, kValues, kthresholdValues = [params['beta']],[params['alpha']],[int(args.kmerSize)],[int(args.kmerThreshold)]

    for x in alpha_values:
            alpha = round(1/float(x),3)
            for beta in beta_values:
                num_cyclopeptide_spectra_dic[(alpha,beta)] = 0
                num_polymer_spectra_dic[(alpha,beta)] = 0

    for kmerSize in kValues: 
        for kmerThreshold in kthresholdValues:
                benchmark_file[(kmerSize,kmerThreshold)] = 1
                reconstructions_file[(kmerSize,kmerThreshold)] = open(output+"_sequencing_reconstructions.txt","w")
    for peptide in peaksnIntensity:
        peaksnIntensity_peptide,precursorMass = peaksnIntensity[peptide],pepMasses[peptide]
        #initilize the analysis
        realPepMass,total,num_raw_spectra,num_spectra_analyzed = initialize_spectrum(verboseprint,pepMasses[peptide],retentions[peptide],charges[peptide],peaksnIntensity_peptide,total,protonMass,num_raw_spectra,num_spectra_analyzed)
        if realPepMass == -1:
            continue
        #Identify cyclospectra
        final_compound,info,standardAutconvCleaned,building_blocks = identify_cyclospectra(
            thisbuildingblock,peaksnIntensity_peptide,e,pepMasses[peptide],retentions[peptide], charges[peptide],representative,pepMasses[peptide],alpha_values,beta_values,
            polymer_repeat_units,fileLines[peptide],writeOriginalSpectra,num_cyclopeptide_spectra_dic,
            generate_convolutions,find_proteinogenic_clusters,output_cyclopeptide_polymers,find_aa_for_denovo,filter_intensity_kmerScore,verboseprint,output,peptide,generate_All_Kmers,realPepMass)
        #de novo sequence cyclospectra
        if final_compound == "cyclopeptide":
            verboseprint("classified: " + str("cyclospectrum"))
            if not args.denovo:
                continue
            allsequenced = denovo_sequence_spectrum(standardAutconvCleaned,kValues,kthresholdValues,denovo_sequence,peaksnIntensity_peptide,realPepMass,building_blocks,e,verboseprint,output_denovo_results,precursorMass,
                retentions[peptide], charges[peptide],peptide,benchmark_file,reconstructions_file,allsequenced,args.aminoThresh)

    nameofrecontfile = output+"_sequencing_reconstructions.txt"
    nameofcyclospecfile = output+"_cyclospectra.mgf"
    # reconstructions_file[(kmerSize,kmerThreshold)].close()
    verboseprint("calculating P-values and cleaning up ...")
    if args.denovo:
        Popen(["sh " +cyclonovoPath+'/scripts/generate_pvalus_recontfile_mgf.sh '+nameofrecontfile+ " " + nameofcyclospecfile + " " + str(output) + " "+ str(cyclonovoPath) +" " +str(e)]
                    ,shell=True)
    
    #report the number of cyclopeptides identified in output_summary.txt
    finalMGF = args.mgfFile.replace(" ", "_")
    argstoprint = [finalMGF, str(len(peaksnIntensity)), num_spectra_analyzed, num_cyclopeptide_spectra_dic[(alpha,beta)], allsequenced]

    verboseprint("CycloNovo succesfully finished. Find the results in "+str(output))






