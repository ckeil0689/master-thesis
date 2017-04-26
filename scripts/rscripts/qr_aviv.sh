#!/bin/bash

np=12

function usage() {
    echo "usage: $(basename $0) [-np <procs>] [-bigmem] file.R" > /dev/stderr
    echo "submits an R script to the rabbot cluster" > /dev/stderr
    echo "  -np: number of processors to use (default ${np})" > /dev/stderr
    echo "  -bigmem: run on a big memory node only" > /dev/stderr
    exit 1
}

unset bigmem
while [[ "$1" =~ ^- ]]; do
    case $1 in
	-usage)
	    usage
	    ;;
	-np)
	    np=$2
	    shift; shift
	    ;;
	-bigmem)
	    bigmem=":bigmem"
	    shift
	    ;;
	*)
	    echo "unknown argument \"$1\", aborting" > /dev/stderr
	    usage
	    ;;
    esac
done

fn=$1
if [ -z "$fn" ]; then echo "must supply R script name" > /dev/stderr; usage; fi

cat <<. | qsub
#!/bin/sh
#PBS -d $(pwd)
#PBS -l nodes=1:ppn=$np${bigmem}
#PBS -j oe
#PBS -o ${fn/.R/.out}
#PBS -N ${fn%.R}
#PBS -m abe
#PBS -M madaraviv@gmail.com

R CMD BATCH --no-restore --no-save $fn
.
