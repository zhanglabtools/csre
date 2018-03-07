#!/bin/sh
mkdir -p log
date > log/tm
MY_JOB=`bsub -J "nlpToSsz[1-23]" < bsubNlpToSsz.sh | grep -Eo "[0-9]+"`
echo $MY_JOB | tee log/jobs
MY_JOB=`bsub -J "allocateStats" -w "done($MY_JOB)" < bsubAllocateStats.sh | grep -Eo "[0-9]+"`
echo $MY_JOB | tee -a log/jobs
MY_JOB=`bsub -J "sszToPd[1-23]" -w "done($MY_JOB)" < bsubSszToPd.sh | grep -Eo "[0-9]+"`
echo $MY_JOB | tee -a log/jobs
MY_JOB=`bsub -J "sh[1-23]" -w "done($MY_JOB)" < bsubSh.sh | grep -Eo "[0-9]+"`
echo $MY_JOB | tee -a log/jobs
MY_JOB=`bsub -J "paraGetWholeCSRE" -w "done($MY_JOB)" < bsubParaGetWholeCSRE.sh | grep -Eo "[0-9]+"`
echo $MY_JOB | tee -a log/jobs
bsub -w "done($MY_JOB)" -J calTime < calTime.sh
