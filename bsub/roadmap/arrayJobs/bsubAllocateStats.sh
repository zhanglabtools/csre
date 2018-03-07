#!/bin/sh
#BSUB -q cpu_dbg
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -W 30
#BSUB -o log/allocateStats/all.out
#BSUB -e log/allocateStats/all.err

dirDmge=/home/zhanglab/cwang/project/dmge
dirLog=$PWD/log
mkdir -p $dirLog/allocateStats
cd $dirDmge
Rscript R/up/allocateStats.R roadmap > $dirLog/allocateStats/all.Rout
