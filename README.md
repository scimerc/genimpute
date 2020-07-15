# genimpute

Integrated tool set to perform quality control and imputation of diploid genotypes

## Quick start

    git clone https://github.com/scimerc/genimpute.git && cd genimpute

Run

    genimpute.sh -h

to get usage instructions. `genimpute.sh`'s default settings deploy a minimal configuration.

Run

    make run_example
  
to run imputation on a small test data set. 

## Configuration

Below is a configuration file example to submit `genimpute.sh` jobs to a slurm scheduler using a
set of reference haplotypes (necessary for imputation). The cluster nodes are expected to have
256GB of memory distributed over 64 computing cores.

    # custom configuration file for slurm

    execommand='sbatch --account=myaccount'
    igroupsize=500
    plinkctime=12:00:00
    plinkmem=12000
    refprefix=hrc.r1-1.ega.grch37
    snodecores=64
    snodemem=256000

## Requirements 

* Linux
* Bash

**NOTE**: Before imputation -- and advisably also before quality control -- reference maps and
haplotypes may have to be pre-processed in order to ensure full compatibility with `genimpute.sh`.
The scripts `mapproc.sh` and `refproc.sh` perform the necessary pre-processing steps on genetic map
and reference haplotypes, respectively.

## Authors

Written by Francesco Bettella and Florian Krull based on a deCODE genetics protocol
adapted by Sudheer Giddaluru.
