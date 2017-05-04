#!/usr/bin/env bash
echo "Ensuring expected directory structure and that all required data from NCBI GEO GSE40918 is present."
# Create all sub-dirs, make no noise if they exist.
parentDir=$(pwd)
supplDir="$parentDir/suppl"
dataDir="$supplDir/data"
mkdir -p $dataDir
mkdir -p $dataDir/{chipseq,deseq,inferelator,rnaseq}

load_from_cell () {
   descr="experiment library"
   if [ $1 -eq "mmc5.xls" ]; then descr="z-score" fi
   fullname=$(basename $1)
   filename="${fullname##*.}"
   if [ ! -f $supplDir/$filename.csv ]; then
      if [ ! -f $supplDir/$1 ]; then
         echo "Downloading $1..."   
         curl -o "$supplDir/$1" "http://www.cell.com/cms/attachment/2007961119/2030652145/$1"
         if [ $? -ne 0 ]; then
            echo "Problem when attempting to download $descr reference table ($1). Stopping."
            exit 1
         else
	    echo "Successfully loaded $1."
         fi
      fi
   fi
}

# Download mmc4 table (experiment library reference).
echo "Checking experiment library reference table (mmc4.xlsx/.csv) from Cell."
load_from_cell "mmc4.xlsx"

# Download mmc5 table (z-score reference). TODO remove mmc5 and create code to extract it from DESeq files
echo "Checking z-score reference table (mmc5.xls) from Cell."
load_from_cell "mmc5.xls"

# Convert mmc4 and mmc5 to CSV using LibreOffice if available
echo "Attempting to automatically convert mmc4.xlsx and mmc5.xls to CSV-files using LibreOffice."
if [ type -P soffice 2>/dev/null ]; then
   # if soffice command is set up on OSX with LibreOffice
   soffice --headless --convert-to csv $supplDir/mmc4.xlsx --outdir $supplDir
   soffice --headless --convert-to csv $supplDir/mmc5.xls --outdir $supplDir
elif [ type -P libreoffice 2>/dev/null ]; then
   # linux with libreoffice install works here     
   libreoffice --headless --convert-to csv $supplDir/mmc4.xlsx --outdir $supplDir
   libreoffice --headless --convert-to csv $supplDir/mmc5.xls --outdir $supplDir
else 
   echo "No LibreOffice command found for conversion of XLSX-files to CSV format."
fi   

if [ ! -f $supplDir/mmc4.csv ]; then
   echo "Failed to automatically convert mmc4.xlsx to mmc4.csv. Please convert the file manually (e.g. using spreadsheet software, such as Excel or LibreOffice --> 'Save As'"
   echo "When converted, run the setup script again. Stopping because mmc4 is required in CSV format. File is located at: $supplDir/mmc4.xlsx"
   exit 1
else
   echo "Successfully converted mmc4.xlsx to mmc4.csv"
fi 

if [ ! -f $supplDir/mmc5.csv ]; then
   echo "Failed to automatically convert mmc5.xls to mmc5.csv. Please convert the file manually (e.g. using spreadsheet software, such as Excel or LibreOffice --> 'Save As'"
   echo "Stopping because the z-score table in mmc5 is required to display differential expression based on RNA-seq data (Th17 vs Th0 at 48h). File is located at: $supplDir/mmc5.xls"
   exit 1
else
   echo "Successfully converted mmc5.xls to mmc5.csv"
fi 

# Download and ensure presence of NCBI GEO GSE40918 raw data
echo "Checking the presence of GEO GSE40918 data."
if [ -z "$(ls -A $dataDir/chipseq)" ] && [ -z "$(ls -A $dataDir/deseq)" ]; then
  read -p "No ChIP-seq or DESeq data from GEO found. Would you like to download it (~133MB)? [y/n] " -n 1 -r
  echo
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    mkdir geotmp
    tmpDir="$parentDir/geotmp"
    ftpDir="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40918/suppl"
    if [ hash wget 2>/dev/null ]; then
       #wget --spider -nc --directory-prefix="$parentDir/geotmp" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40918&format=file"
       wget -nc -l1 --directory-prefix="$tmpDir/" "$ftpdir/*"
    else
       filelist=$(curl -i -l "$ftpDir/") 
       for i in $filelist; do 
	   echo $ftpDir/${i}
	   curl -o "$tmpDir/${i}" $ftpDir/${i} 
       done
    fi
    echo "Exit code: $?"
    if [ $? -ne "0" ]; then
	echo "Download of GEO files failed. Stopping."
        exit 1
    fi	
    # Extract
    echo "Extracting and moving ChIP-seq and RNA-seq files..."
    filepath="$tmpDir/GSE40918_RAW.tar"
    if [ -f $filepath ]; then
       tar -xvf $filepath -C $tmpDir
       gunzip -k $tmpDir/*.gz
       # move all files ending in '_genes.txt' to /data/chipseq/
       mv $tmpDir/*_genes.txt $dataDir/chipseq/
       # move all files ending in '_genes.expr.txt' to /data/rnaseq/
       mv $tmpDir/*_genes.expr.txt $dataDir/rnaseq/
       # Move DESeq and Inferelator files
       mv $tmpDir/GSE40918_Th17*.wt.vs.Th17*ko_*.txt $dataDir/deseq/
       mv $tmpDir/GSE40918_Inferelator*.txt $dataDir/inferelator/
       # Move README to /suppl/
       mv $tmpDir/GSE40918_README_supplementary_files.txt $supplDir
       # delete tmp folder and all of its contents
       rm -rf $tmpDir
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

