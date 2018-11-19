# CycloNovo: Algorithm for de novo cyclopeptide analysis of high-resolution tandem mass spectra

CycloNovo is a new algorithm that identifies spectra generated from cyclopeptides in large mass spectrometry datasets. Moreover, CycloNovo can de novo sequence the cyclopeptides represented by identified cyclospectra.

Developed in University of California San Diego, La Jolla, CA, USA


Dependencies:

- Dereplicator (http://cab.spbu.ru/software/dereplicator/)


System requirements:

- macOS
- Python 2.7


Usage examples: 

     cyclonovo_classify_denovo.py -s data/surugamide_spectrum.mgf -o test_output \
     --denovo -k 5 --verbosity -e 0.015 --kmer_threshold 2 


For the full list of available options please run

     cyclonovo.py -h


Output:
* cyclonovo_classification_report.txt              list spectra analyzed with scores related to identifying cyclospectra\
* cyclonovo_cyclospectra.mgf                       mgf file contaning found cyclospectra\
* cyclonovo_sequencing_reconstructions.txt         reconstructions if de novo sequencing option is used





