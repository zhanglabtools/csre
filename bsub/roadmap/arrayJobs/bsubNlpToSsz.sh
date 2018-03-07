#!/bin/sh
#BSUB -q cpuII
#BSUB -n 24
#BSUB -R "span[ptile=24] rusage[mem=2048]"
#BSUB -W 6:00
#BSUB -o log/nlpToSsz/%I/%I.out
#BSUB -e log/nlpToSsz/%I/%I.err

dirDmge=/home/zhanglab/cwang/project/dmge
dirLog=$PWD/log
mkdir -p $dirLog/nlpToSsz/$LSB_JOBINDEX
cd $dirDmge
Rscript R/up/nlpToSsz.R roadmap $LSB_JOBINDEX > $dirLog/nlpToSsz/$LSB_JOBINDEX/$LSB_JOBINDEX.Rout
