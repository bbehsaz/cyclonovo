PYTHON_VERSION := $(shell python --version 2>&1)
print_spectrum_help := $(shell print_spectrum -h 2>&1)
print_score_help := $(shell print_score -h 2>&1)

# Run test example for CycNovo.
# surugamide_test/summary.tsv : data/surugamide_spectrum.mgf
# 	echo "Running tests surugamdie";	python cyclonovo.py -s data/surugamide_spectrum.mgf -o surugamide_test --monomers standard --preprocess -e 0.02 --verbosity;


.PHONY: all

all:
ifdef PYTHON_VERSION
	@echo "Found version $(PYTHON_VERSION)"
else
	@echo "python not in PATH.";	exit 1;
endif
ifdef print_spectrum_help
	@echo "print_spectrum found."
else
	@echo $(print_spectrum_help)
	@echo "print_spectrum not in PATH.";	exit 1;
endif
ifdef print_score_help
	@echo "print_score found."
else
	@echo $(print_spectrum_help)
	@echo "print_score not in PATH.";	exit 1;
endif
	@echo "Running CycloNovo tests ..."
	# testsurugsamide:
	cp ./scripts/cyclonovo.py cyclonovo
	chmod a+wrx cyclonovo
	@echo "Running tests: surugamdie";	./cyclonovo -s data/surugamide_spectrum.mgf -o surugamide_test --monomers standard -e 0.015 --verbosity --denovo
