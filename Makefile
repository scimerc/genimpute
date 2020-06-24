BASEDIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

all: 3rdparty

prepare_offline:
	$(MAKE) -C lib/3rd download

.PHONY: 3rdparty
3rdparty:
	$(MAKE) -C lib/3rd all

.PHONY: distclean
distclean:
	$(MAKE) -C lib/3rd clean

.PHONY: run_example
run_example: 3rdparty
	LC_NUMERIC=C progs/genimpute.sh -c toy/config.txt toy/testset_all.bcf.gz
