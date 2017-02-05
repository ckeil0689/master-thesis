#!/bin/bash
# Substitute the working directory and data directory to run trimmomatic and STAR over all files in the directory
# Important: wd and data directory need to be identical, otherwise it can't find the files

FILES=`ls /work/other/bsteglich/data/`


for FILE in $FILES 
do
  qsub -clear -V -S /bin/bash -q hpc.q -l mem_free=30G -j n -w w -N trim_star -o /work/other/bsteglich/output -e /work/other/bsteglich/no_backup/error_logs -wd /work/other/bsteglich/data ./scripts/160927_Trim_STAR.sh $FILE test
  JOBNAME="trim_star"
  RUNNING=`qstat -u bsteglich | grep $JOBNAME`
  while [ "$RUNNING" != '' ]
  do
    sleep 60
    RUNNING=`qstat -u bsteglich | grep $JOBNAME`
  done
  echo "done with $FILE"
done
