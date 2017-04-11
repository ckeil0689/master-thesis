#!/usr/bin/env bash
echo "Ensuring all required data is present."
parentDir=$(pwd)
supplDir="$parentDir/suppl"
if [ -d $supplDir ]
then
  echo "'suppl' subdirectory exists."
else
  echo "$supplDir not found." 
  $(mkdir suppl)
fi


