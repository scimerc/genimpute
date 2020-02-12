# genimpute

integrated tool set to perform quality control and imputation of diploid genotypes

this tool set was written by Florian Krull and me loosely based on the deCODE genetics protocol
adapted by Sudheer Giddaluru.

**NOTE**: the tool set is currently designed to run in linux.

run

    genimpute.sh -h

to get usage instructions.
`genimpute.sh`'s default settings deploy a minimal configuration.
below is a configuration file example to submit `genimpute.sh` jobs to a slurm scheduler using a
set of reference haplotypes (necessary for imputation). the cluster nodes are expected to have
256GB of memory distributed over 64 computing cores.

    # custom configuration file for slurm

    execommand='sbatch --account=myaccount'
    igroupsize=500
    plinkctime=12:00:00
    plinkmem=12000
    refprefix=hrc.r1-1.ega.grch37
    snodecores=64
    snodemem=256000

**NOTE**: before imputation -- and advisably also before quality control -- reference maps and
haplotypes may have to be pre-processed in order to ensure full compatibility with `genimpute.sh`.
the scripts `mapproc.sh` and `refproc.sh` perform the necessary pre-processing steps on genetic map
and reference haplotypes, respectively.

