# CycloNovo: Algorithm for de novo cyclopeptide analysis of high-resolution tandem mass spectra

CycloNovo is a new algorithm that identifies spectra generated from cyclopeptides in large mass spectrometry datasets. CycloNovo can also de novo sequence the cyclopeptides represented by identified cyclospectra.

Developed in University of California San Diego, La Jolla, CA, USA


Dependencies (add binaries to PATH):

	- NPDtools (https://github.com/ablab/npdtools)

System requirements:

	- Linux or macOS
	- Python 2.7
	- GNU sed 

To configure/test run:

     make Makefile

Usage examples: 

     cyclonovo -s data/surugamide_spectrum.mgf -o test_output \
     --denovo -k 5 --verbosity -e 0.015 --kmer_threshold 2 


For the full list of available options please run

     cyclonovo -h


Output:

* cyclonovo_classification_report.txt:              list spectra analyzed with scores related to identifying cyclospectra
* cyclonovo_cyclospectra.mgf:                       mgf file contaning found cyclospectra
* cyclonovo_sequencing_reconstructions.txt:         reconstructions if de novo sequencing option is used


For the latest version and the web application see (login to GNPS first):
https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=%7B%22workflow%22:%22CYCLONOVO%22%7D
