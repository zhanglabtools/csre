#!/bin/sh
#BSUB -q cpuII
#BSUB -n 128
#BSUB -R "span[ptile=22] rusage[mem=10240]"
#BSUB -W 6:00
#BSUB -o log/sh/%I/%I.out
#BSUB -e log/sh/%I/%I.err

dirDmge=/home/zhanglab/cwang/project/dmge
dirLog=$PWD/log
mkdir -p $dirLog/sh/$LSB_JOBINDEX
cd $dirDmge
Rscript R/up/sh.R roadmap $LSB_JOBINDEX > $dirLog/sh/$LSB_JOBINDEX/$LSB_JOBINDEX.Rout
