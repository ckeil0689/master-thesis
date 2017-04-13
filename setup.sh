#!/usr/bin/env bash

echo "Ensuring expected directory structure and that all required data from NCBI GEO GSE40918 is present."
# Create all sub-dirs, make no noise if they exist.
parentDir=$(pwd)
supplDir="$parentDir/suppl"
dataDir="$supplDir/data"
mkdir -p $dataDir
mkdir -p $dataDir/{chipseq,deseq,inferelator,rnaseq}

# Download mmc4 table (experiment library reference). TODO remove mmc5 and create code to extract it from DESeq files
echo "Checking experiment library reference table (mmc4.xlsx) from Cell."
wget -nc --directory-prefix=$supplDir "http://www.cell.com/cms/attachment/2007961119/2030652145/mmc4.xlsx"
if [ $? -ne 0 ];
then
  echo "Problem when attempting to download experiment library reference table (mmc4.xlsx). Stopping."
  exit 1
fi

echo "Checking the presence of GEO GSE40918 data."
if [ -z "$(ls -A $dataDir/chipseq)" ] && [ -z "$(ls -A $dataDir/deseq)" ]; then
  read -p "No ChIP-seq or DESeq data from GEO found. Would you like to download it (~133MB)? [y/n] " -n 1 -r
  echo
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    mkdir geotmp
    #wget --spider -nc --directory-prefix="$parentDir/geotmp" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40918&format=file"
    wget -nc -l1 --directory-prefix="$parentDir/geotmp" "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40918/suppl/*"
    # Extract
    tmpDir="$parentDir/geotmp"
    filepath="$tmpDir/GSE40918_RAW.tar"
    if [ -f $filepath ]; then
       echo "Extracting and moving ChIP-seq and RNA-seq files..."
       tar -xvf $filepath -C $tmpDir
       gunzip -k $tmpDir/*.gz
       # move all files ending in '_genes.txt' to /data/chipseq/
       mv $tmpDir/*_genes.txt $dataDir/chipseq/
       # move all files ending in '_genes.expr.txt' to /data/rnaseq/
       mv $tmpDir/*_genes.expr.txt $dataDir/rnaseq/
       # Move DESeq and Inferelator files
       mv $tmpDir/GSE40918_Th17*.wt.vs.Th17*ko_*.txt $dataDir/deseq/
       mv $tmpDir/GSE40918_Inferelator*.txt $dataDir/inferelator/
       # delete tmp folder and all of its contents
       #rm -rf $parentDir/geotmp/
       echo "Would delete $tmpDir"
    else
      echo "Could not find GSE40918_RAW.tar. Download did not finish as expected. Stopping."
      exit 1
    fi
  else
    echo "Setup not completed, no GEO data available."
    exit 1
  fi
else
  echo "Files present in /data/chipseq/ and /data/deseq/. Not downloading anything."
fi
echo "Done! All set up and ready to go. To generate a network with the files present in /data/, move into /scripts/ and run 'master.R'"

