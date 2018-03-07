#!/bin/sh
#BSUB -q cpu_dbg
#BSUB -n 24
#BSUB -R "span[ptile=24] rusage[mem=9216]"
#BSUB -W 30
#BSUB -o log/paraGetWholeCSRE/all.out
#BSUB -e log/paraGetWholeCSRE/all.err

regionHeight=10
regionWidth=96
dirDmge=/home/zhanglab/cwang/project/dmge
dirLog=$PWD/log
mkdir -p $dirLog/paraGetWholeCSRE
cd $dirDmge
Rscript R/up/paraGetWholeCSRE.R roadmap $regionHeight $regionWidth > $dirLog/paraGetWholeCSRE/${regionHeight}_$regionWidth.Rout
