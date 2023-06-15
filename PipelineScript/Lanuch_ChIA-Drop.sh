#!/bin/bash

print_help()
{
    echo ''
	echo -e "Usage: $0 -t [STR] -n [STR] -p [DIR] -d [STR] -r [STR] -@ [INT] "
    echo -e "\t-t Experiment type: sprite, scsprite, rdsprite, chiadrop, scrna, scatac, scrnawscatac"
	echo -e "\t-n Library name"
    echo -e "\t-p Pipeline processing tool directory"
    echo -e "\t-d DNA alignment file"
    echo -e "\t-r RNA alignment file"
    echo -e "\t-b Identifier to cluster, for example: CELL:CB-COMPLEX:CP"
    echo -e "\t-s Identifier relationship, for example: CELL:COMPLEX"
    echo -e "\t-m Distance to merge in ChIA-Drop"
    echo -e "\t-e Extend from 5' end"
    echo -e "\t-c Configuration file of barcode."
    echo -e "\t-x Genome size file name under BARP directory"
	echo -e "\t-@ Thread to use"
}

if [[ $1 == '' ]]
then
    print_help
    exit 1
fi

LIB_NAME="CHDP"
ChIADrop_FILES="-1 R1.fq -2 R2.fq"
THREAD=10
ChIADrop_CONFIG="-"
BWA_REF_GENOME="REF_OF_BWA"
GENOME_SIZE_FILE="dm3.size.txt"

BARP_DIR="REF_OF_ScSmOP"

while getopts x:c:t:p:n:@:d:r:g:h flag
do
    case "${flag}" in 
        c)  config=${OPTARG};;
        t)  expe_type=${OPTARG};;
        p)  pipe_dir=${OPTARG};;
        n)  name=${OPTARG};;
        d)  dna_alignment=${OPTARG};;
        r)  rna_alignment=${OPTARG};;
        @)  num_thread=${OPTARG};;
        g)  star_ref=${OPTARG};;
        x)  genomesize=${OPTARG};;
        ? | h) print_help
            exit 1;;
    esac
done


