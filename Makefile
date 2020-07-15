.PHONY: run_example
run_example:
	LC_NUMERIC=C progs/genimpute.sh -c toy/config.txt toy/testset_all.bcf.gz
