#! /bin/bash
#$ -l mem=10G,time=4:: -cwd -S /bin/bash -N famSKAT

export R_LIBS_USER=~/local/R/hpc
Rscript  --max-ppsize=500000 qsub.R $SGE_TASK_ID

