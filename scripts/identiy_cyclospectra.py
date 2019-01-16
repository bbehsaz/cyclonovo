
def identify_cyclospectra(thisbuildingblock,peaksnIntensity_peptide,e,precursorMass,retention,charge,representative,pepMass,
    alpha_values,beta_values,polymer_repeat_units,lines2print,writeOriginalSpectra,num_cyclopeptide_spectra_dic,
    generate_convolutions,find_proteinogenic_clusters,output_cyclopeptide_polymers,find_aa_for_denovo,filter_intensity_kmerScore,verboseprint,output,peptide,generate_All_Kmers,realPepMass):
        finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge, standardAminoMasses, e = generate_convolutions(
                    thisbuildingblock,peaksnIntensity_peptide,e,int(precursorMass/10.0)+10,charge,representative,pepMass)
        
        for x in alpha_values:
                    alpha = round(1/float(x),3)
                    for beta in beta_values:                      
                        suffix = "_"+str(x)+"_"+str(beta)
                        num_protenogenic,N, foundConvolutions, distancesCleaned, thresholdValue, stand_aa_peak_pairs = find_proteinogenic_clusters(
                            finalClusters, ALLcleanedAutoconvolutions, precursorMass, charge,standardAminoMasses,e,alpha,beta)                     
                        compound_type,convolutions, numStnd = output_cyclopeptide_polymers(
                            num_protenogenic,N, foundConvolutions, distancesCleaned,polymer_repeat_units,charge,output+suffix, 
                            pepMass, retention,lines2print,writeOriginalSpectra,len(peaksnIntensity_peptide),peptide, alpha, beta,int(thresholdValue),stand_aa_peak_pairs)
        beta=beta_values[0]
        final_compound = "unclassified"

        if compound_type == 'polymer':
            final_compound = 'notcyclic'
        # if compound_type == "cyclopeptide":
        if True:
            intensities = peaksnIntensity_peptide
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
            num_cyclopeptide_spectra_dic[(alpha,beta)] += 1
            cyclopeptide_spectra_file = open(output+"_cyclospectra.mgf","a")
            writeOriginalSpectra(lines2print,cyclopeptide_spectra_file)
        # if final_compound == 'polymer':
        #     num_polymer_spectra_dic[(alpha,beta)] += 1

        if len(standardAutconvCleaned)>2:
            building_blocks = [key for key, value in standardAutconvCleaned.iteritems() if value>1]
            # sorted_building_blockes = [x for x in standardAutconvCleaned.iteritems() ]
        # info = [precursorMass, retentions[peptide],charges[peptide], peptide, final_compound, percentProteinoPeaks, highetScoredKmer, percentPolymerPeaks]
        if final_compound != "cyclopeptide":
            final_compound = "nc"
        from operator import itemgetter
        building_blocks_multiplicity =  "["+ ",".join([ "("+str(round(key,2))+","+str(value)+")" for key, value in sorted(standardAutconvCleaned.items(), key=itemgetter(1), reverse=True) if value>2]) +"]"
        info = [precursorMass, retention,charge, peptide, final_compound, highetScoredKmer, 
                percentProteinoPeaks, topkmer, building_blocks_multiplicity ] #reports classification measure for each spectrum
        compound_type_reports_file = open(output+"_classification_report.txt","a")
        compound_type_reports_file.write("\t".join(str(z) for z in info)+"\n")

        return final_compound,info,standardAutconvCleaned,building_blocks