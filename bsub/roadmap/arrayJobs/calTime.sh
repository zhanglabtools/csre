#!/bin/sh
#BSUB -J calTime
#BSUB -q cpu_dbg
STARTTIME=`date +%s -d "$(cat log/tm)"`
THEEND=`date`
ENDTIME=`date -d "$THEEND" +%s`
TIMESAPN=`expr "(" $ENDTIME - $STARTTIME ")" / 60`
echo $THEEND >> log/tm
echo "Time span: $TIMESAPN mins" >> log/tm
cd log
echo "Sum of time span of each job: `grep -ri total *|grep -Eo \"[0-9]+\.[0-9]+ \"|awk '{sum+=$1}END{print sum}'` mins" >> tm
