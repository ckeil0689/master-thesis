#!/usr/bin/env bash

echo "Ensuring all required data from NCBI GEO GSE40918 is present."
# Check if in correct environment (setup.sh could have been moved)
if [ ! `git rev-parse --git-dir` ];
then
  echo "Not in git repository."
  exit 1
fi

# Create all sub-dirs, make no noise if they exist.
parentDir=$(pwd)
supplDir="$parentDir/suppl"
dataDir="$supplDir/data"
mkdir -p $dataDir
mkdir -p $dataDir/{chipseq, deseq, inferelator, rnaseq}


echo "Done! All set up and ready to go. To generate a network with the files present in /data/, move into /scripts/ and run 'master.R'"
