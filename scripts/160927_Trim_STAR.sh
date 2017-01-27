#! /bin/bash
set -e
set -u
set -o pipefail

# This script can be used to loop over a set of fastq files (using a submit_all_jobs.sh script), running them first through trimmomatic, then STAR. For both, a file with key stats is produced. fastq and SAM files are not stored, only a ReadsPerGene.Out.tab file is produced for each fastq file.

file_path=$1
out_dir=$2
index=/work/other/bsteglich/no_backup/STAR_index/

## check if files and directories are present. The script appends information about the the Trimmomatic and STAR results to the two results files, so they need to be created for the first run, after that the script only adds extra lines to them

if [ ! -d "$out_dir" ]; then mkdir "$out_dir"
fi

if [ ! -d "$out_dir/tmp" ]; then mkdir "$out_dir/tmp"
fi

if [ ! -e "$out_dir/trimmomatic_results.txt" ]; then 
	touch $out_dir/trimmomatic_results.txt
	echo -e "Sample\tSurviving_reads\tDropped_reads" >> $out_dir/trimmomatic_results.txt
fi

if [ ! -e "$out_dir/star_results.txt" ]; then 
	touch $out_dir/star_results.txt
	echo -e "Sample\tUnique_mapped_reads\tMultiple_mapped_reads\tUnmapped_reads" >> $out_dir/star_results.txt
fi

## run trimmomatic, save the resulting fastq file in tmp directory, add info about surviving and dropped reads to results file.

filename=$(basename $file_path)
samplename="${filename%.fastq.gz}"

java -jar /work/gi/software/trimmomatic-0.33/trimmomatic-0.33.jar SE -phred33 -threads 8 $file_path $out_dir/tmp/$filename ILLUMINACLIP:/work/gi/software/trimmomatic-0.33/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> $out_dir/tmp/stderror.txt

surv=`cat $out_dir/tmp/stderror.txt | awk 'NR==5' | awk -F "[()]" '{ for (i=2; i<NF; i+=2) print $i }' | head -1`
drop=`cat $out_dir/tmp/stderror.txt | awk 'NR==5' | awk -F "[()]" '{ for (i=2; i<NF; i+=2) print $i }' | tail -1`
echo -e $samplename"\t"$surv"\t"$drop >> $out_dir/trimmomatic_results.txt

## run STAR on the fastq file produced by trimmomatic, add info about the alignment to results file. Note: STAR will only produce a counts per gene file in this way, if you need a SAM or BAM file, change the --outSAMmode setting.

/work/gi/software/STAR/bin/Linux_x86_64/STAR  --runThreadN 8 --genomeDir $index --readFilesIn $out_dir/tmp/$filename --quantMode GeneCounts --outSAMmode None --readFilesCommand gunzip -c --outFileNamePrefix $out_dir/tmp/$samplename.

unique=`cat $out_dir/tmp/$samplename.Log.final.out | awk 'NR==10' | grep -o -E '[0-9.]+'`
multiple_1=`cat $out_dir/tmp/$samplename.Log.final.out | awk 'NR==25' | grep -o -E '[0-9.]+'`
multiple_2=`cat $out_dir/tmp/$samplename.Log.final.out | awk 'NR==27' | grep -o -E '[0-9.]+'`
multiple_all=$(awk "BEGIN {print $multiple_1 + $multiple_2; exit}")
unmapped_1=`cat $out_dir/tmp/$samplename.Log.final.out | awk 'NR==29' | grep -o -E '[0-9.]+'`
unmapped_2=`cat $out_dir/tmp/$samplename.Log.final.out | awk 'NR==30' | grep -o -E '[0-9.]+'`
unmapped_3=`cat $out_dir/tmp/$samplename.Log.final.out | awk 'NR==31' | grep -o -E '[0-9.]+'`
unmapped_all=$(awk "BEGIN {print $unmapped_1 + $unmapped_2 + $unmapped_3; exit}")
echo -e $samplename"\t"$unique"\t"$multiple_all"\t"$unmapped_all >> $out_dir/star_results.txt

## Move the ReadsPerGene file out of the tmp directory to keep it, then remove the tmp directory.

mv $out_dir/tmp/$samplename.ReadsPerGene.out.tab $out_dir

rm -r $out_dir/tmp/
