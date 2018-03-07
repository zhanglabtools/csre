#!/bin/sh
#BSUB -q cpuII
#BSUB -n 24
#BSUB -R "span[ptile=24] rusage[mem=1024]"
#BSUB -W 6:00
#BSUB -o log/sszToPd/%I/%I.out
#BSUB -e log/sszToPd/%I/%I.err

dirDmge=/home/zhanglab/cwang/project/dmge
dirLog=$PWD/log
mkdir -p $dirLog/sszToPd/$LSB_JOBINDEX
cd $dirDmge
Rscript R/up/sszToPd.R roadmap $LSB_JOBINDEX > $dirLog/sszToPd/$LSB_JOBINDEX/$LSB_JOBINDEX.Rout
