# genimpute

integrated tool set to perform quality control and imputation of diploid genotypes

this tool set was written by Florian Krull and me loosely based on the deCODE genetics protocol
adapted by Sudheer Giddaluru.

NOTE: the tool set is currently designed to run in linux.

run

    genimpute.sh -h

to get usage instructions.

**NOTE**: before imputation -- and advisably also before quality control -- reference maps and
haplotypes may have to be pre-processed in order to ensure full compatibility with `genimpute.sh`.
the scripts `genproc.sh` and `refproc.sh` perform the necessary pre-processing steps on genetic map
and reference haplotypes, respectively.

