def denovo_sequence_spectrum(standardAutconvCleaned,kValues,kthresholdValues,denovo_sequence,peaksnIntensity_peptide,realPepMass,building_blocks,e,verboseprint,
    output_denovo_results,precursorMass,retention,charge,peptide,benchmark_file,reconstructions_file,allsequenced,aminoThresh):         
    for kmerSize in kValues: 
        for kmerThreshold in kthresholdValues: 
            verboseprint("@@@")
            verboseprint("De novo sequence ... \nkmer: size {}, threshold {}".format(kmerSize,kmerThreshold))
            if len(building_blocks)<3:
                verboseprint("Number of predicted cyclopeptidic amino acids is low")
                break
            verboseprint("Predicted Amino Acids:")
            building_blocks = []
            # for p in sorted(standardAutconvCleaned.iteritems(), key=itemgetter(1), reverse=True):
            for p in standardAutconvCleaned.items():
                if p[1]>aminoThresh:
                    verboseprint("{}\t{}".format(round(p[0],2),p[1]))
                    building_blocks.append(p[0])
            candidateSequences = denovo_sequence(
                peaksnIntensity_peptide, realPepMass, building_blocks,kmerSize,e,e,kmerThreshold,verboseprint)
            allsequenced += output_denovo_results(candidateSequences, precursorMass, retention, charge, peptide,
                benchmark_file[(kmerSize,kmerThreshold)], kmerSize, kmerThreshold,reconstructions_file[(kmerSize,kmerThreshold)])
            # nameofrecontfile = output+"_sequencing_reconstructions.txt"
            # nameofcyclospecfile = output+"_cyclospectra.mgf"
            # reconstructions_file[(kmerSize,kmerThreshold)].close()
    return allsequenced