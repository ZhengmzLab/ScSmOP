#!/bin/bash
start=$SECOND
#=============================
print_help()
{
    echo ''
	echo -e "Usage: $0 -t [STR] -n [STR] -p [DIR] -d [STR] -r [STR] -1 [STR] -2 [STR] -@ [INT] -f [DIR] -g [DIR] -k [STR]"
    echo -e "\t-t Experiment type: sprite, scsprite, rdsprite, chiadrop, scrna, scatac, scrnawscatac"
	echo -e "\t-n Library name"
    echo -e "\t-p Pipeline processing tool directory"
    echo -e "\t-d Alignment tool used to align DNA, [bwa]"
    echo -e "\t-r Alignment tool used to align RNA, [STAR]"
	echo -e "\t-1 Read 1 files, separate with ",""
    echo -e "\t-2 Read 2 files, separate with ",""
    echo -e "\t-f BWA reference genome directory"
    echo -e "\t-g STAR reference genome directory"
    echo -e "\t-k Treat NDNR file as DNA or RNA.[DNA]"
	echo -e "\t-@ Thread to use"
    echo -e "\t-s single-cell rna-seq"
}

if [[ $1 == '' ]]
then
    print_help
    exit 1
fi

NDNR_AS="DNA"
run_dir=$(pwd)
name="-"
star_aligner="STAR"
bwa_aligner="bwa"
read_1_list="-"
read_2_list="-"
num_thread=1
dna_aligner=${bwa_aligner}
rna_aligner=${star_aligner}
scrna_flag=0

while getopts t:p:n:f:@:1:2:g:d:r:k:sh flag
do
    case "${flag}" in 
        t)  expe_type=${OPTARG};;
        p)  pipe_dir=${OPTARG};;
        n)  name=${OPTARG};;
        d)  dna_aligner=${OPTARG};;
        1)  read_1_list=${OPTARG};;
        2)  read_2_list=${OPTARG};;
        @)  num_thread=${OPTARG};;
        r)  rna_aligner=${OPTARG};;
        f)  bwa_ref=${OPTARG};;
        k)  NDNR_AS=${OPTARG};;
        g)  star_ref=${OPTARG};;
        s)  scrna_flag=1;;
        ? | h) print_help
            exit 1;;
    esac
done


if [[ $expe_type == "scsprite" ]]
then
    dna_aligner="star"
fi

FinishSuccess=1

PreDir=$(pwd)

if [[ ( ! ${bwa_ref} =~ ^/ ) && ( ! ${bwa_ref} =~ ^~ ) ]]
then
    bwa_ref=$( echo -e "${PreDir}/${bwa_ref}" )
fi

if [[ ( ! ${star_ref} =~ ^/ ) && ( ! ${star_ref} =~ ^~ ) ]]
then
    star_ref=$( echo -e "${PreDir}/${star_ref}" )
fi

if [[ ! -d 02.ReadAlign ]]
then
    mkdir 02.ReadAlign
    cd 02.ReadAlign
else 
    echo -e "\nDirectory 02.ReadAlign exist, will cover the files generated before in 5 seconds...";
    sleep 5
    cd 02.ReadAlign
fi
echo -e "\nSequence alignment..."
echo -e "Working directory: $(pwd)"
echo ""

# dna mapping
function bwa_map
{
    local FileArray
    FileArray=($(echo "$@"))
    if [[ ! -f ${bwa_ref} ]]
    then
        echo "BWA reference did not specified."
        exit 1
    fi
    if [[ ${#FileArray[@]} == 1 ]]
    then
        echo "bwa mapping single end read."
        echo "${pipe_dir}/Tools/${bwa_aligner} mem ${bwa_ref} ${FileArray} -o ${name}_DNA.bam -t ${num_thread} "
        ${pipe_dir}/Tools/${bwa_aligner} mem ${bwa_ref} ${FileArray} -o ${name}_DNA.bam -t ${num_thread} 
        if [[ $? != 0 ]]
        then
            FinishSuccess=0
        fi
    elif [[ ${#FileArray[@]} == 2 ]]
    then
        echo "bwa mapping pair end read."
        echo "${pipe_dir}/Tools/${bwa_aligner} mem -t ${num_thread} -o ${name}_DNA.bam ${bwa_ref} ${FileArray[0]} ${FileArray[1]}"
        ${pipe_dir}/Tools/${bwa_aligner} mem -t ${num_thread} -o ${name}_DNA.bam ${bwa_ref} ${FileArray[0]} ${FileArray[1]}
        if [[ $? != 0 ]]
        then
            FinishSuccess=0
        fi
    else
        echo "Only single end or pair end read is supported."
        exit 2
    fi
}

function star_map
{ 
    local GenomeType
    GenomeType=$@

    if [[ ! -d ${star_ref} ]]
    then
        echo "STAR reference did not specified."
        exit 1
    fi

    if [[ ${GenomeType} == "RNA" ]]
    then
        echo "star mapping single end read."
        ${pipe_dir}/Tools/${rna_aligner} --genomeDir ${star_ref} --readFilesIn ${read1_file} ${read2_file}\
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFileNamePrefix ${name} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode TranscriptomeSAM GeneCounts \
            --outWigType bedGraph \
            --runThreadN ${num_thread}
        if [[ $? != 0 ]]
        then
            FinishSuccess=0
        fi
    elif [[ ${GenomeType} == "DNA" ]]
    then
        echo -e "Mapping ${GenomeType} with STAR.\n"
        ${pipe_dir}/Tools/${rna_aligner} --outFilterMultimapNmax 50 \
        --outFilterScoreMinOverLread 0.30 \
        --outFilterMatchNminOverLread 0.30 \
        --outFilterIntronMotifs None \
        --alignIntronMax 50000 \
        --alignMatesGapMax 1000 \
        --genomeLoad NoSharedMemory \
        --outReadsUnmapped Fastx \
        --alignIntronMin 80 \
        --alignSJDBoverhangMin 5 \
        --sjdbOverhang 100 \
        --genomeDir ${star_ref} \
        --readFilesIn ${read1_file} ${read2_file}\
        --outFileNamePrefix ${name}.Barcoded. \
        --outSAMtype BAM Unsorted \
        --outSAMattributes All \
        --limitOutSJcollapsed 10000000 \
        --limitIObufferSize=300000000 320000000 \
        --runThreadN ${num_thread}
        if [[ $? != 0 ]]
        then
            FinishSuccess=0
        fi
    elif [[ ${GenomeType} == "SCRNA" ]]
    then
        echo -e "Mapping ${GenomeType} with STAR.\n"
        ${pipe_dir}/Tools/${rna_aligner} --readFilesIn ${read1_file} ${read2_file}\
            --genomeDir ${star_ref} \
            --soloType BARP \
            --soloBarpIdentifier CELL \
            --soloCBwhitelist None \
            --runThreadN ${num_thread} \
            --soloUMIfiltering MultiGeneUMI_CR \
            --soloUMIdedup 1MM_CR \
            --outFilterScoreMin 30 \
            --clipAdapterType CellRanger4 \
            --soloBarcodeReadLength 0 \
            --outFileNamePrefix ${name}
        if [[ $? != 0 ]]
        then
            FinishSuccess=0
        fi
    else
        echo "No mapping.\n"
    fi
}


if [[ ${read_1_list} == "-" ]]
then
    echo "Read FASTQ files from 01.BarcodeIden"
    if read_1_list=`ls ${PreDir}/01.BarcodeIden/*_1.fq 2>/dev/null`
    then
        r1_file_list=(${read_1_list})
    fi
else
    r1_file_list=(${read_1_list//,/ })
fi

if [[ ${read_2_list} == "-" ]]
then
    echo "Read FASTQ files from 01.BarcodeIden"
    if read_2_list=`ls ${PreDir}/01.BarcodeIden/*_2.fq 2>/dev/null`
    then
        r2_file_list=(${read_2_list})
    fi
else
    r2_file_list=(${read_2_list//,/ })
fi

pair_end_mode=0

echo "Read 1 List : ${r1_file_list[@]} Read 2 List : ${r2_file_list[@]}"

if [[ (${#r1_file_list[@]} != 0) && (${#r2_file_list[@]} != 0) ]]
then
    echo "Input files are pair end FASTQ"
    pair_end_mode=1
    if [[ ${#r1_file_list[@]} != ${#r2_file_list[@]} ]]
    then
        echo "Input files are pair end but the number of Read 1 and Read 2 don't match."
    fi
else
    echo "Input files are single end FASTQ"
    if [[ (${#r2_file_list[@]} != 0) ]]
    then
        r1_file_list=${r2_file_list}
    fi
fi

for((i=0;i<${#r1_file_list[@]};i++))
do
{
    read1_file=${r1_file_list[i]}
    if [[ ! ${read1_file} =~ ^/ ]]
    then
        read1_file="../${read1_file}"
        if [[ ! -f ${read1_file} ]]
        then
            echo "${read1_file} do not exist, exiting"
            exit 2
        fi
    fi

    if [[ $pair_end_mode == 1 ]]
    then
        read2_file=${r2_file_list[i]}
        if [[ ! ${read2_file} =~ ^/ ]]
        then
            read2_file="../${read2_file}"
            if [[ ! -f ${read2_file} ]]
            then
                echo "${read2_file} do not exist, exiting"
                exit 2
            fi

        fi
    fi

    echo "Mapping: ${read1_file} ${read2_file}"

    if [[ ${scrna_flag} == 1 ]]
    then
        star_map "SCRNA"
    else
        if [[ ( (${read1_file} =~ DNA) || ( (${read1_file} =~ NDNR) && (${NDNR_AS} == "DNA" || ${NDNR_AS} == "dna" ) ) ) && ( (${dna_aligner} == "bwa") || (${dna_aligner} == "BWA") ) ]]
        then

            if [[ ${pair_end_mode} == 1 ]]
            then
                bwa_map ${read1_file} ${read2_file}
            else
                bwa_map ${read1_file}
            fi
        elif [[ ( (${read1_file} =~ DNA) || ( (${read1_file} =~ NDNR) && (${NDNR_AS} == "DNA" || ${NDNR_AS} == "dna" ) ) ) && ( (${dna_aligner} == "star") || (${dna_aligner} == "STAR") ) ]]
        then
            star_map "DNA"
        fi

        if [[ ( (${read1_file} =~ RNA) || ( (${read1_file} =~ NDNR) && (${NDNR_AS} == "RNA" || ${NDNR_AS} == "rna" ) ) ) && ( (${rna_aligner} == "bwa") || (${rna_aligner} == "BWA") ) ]]
        then
            if [[ ${pair_end_mode} == 1 ]]
            then
                bwa_map ${read1_file} ${read2_file}
            else
                bwa_map ${read1_file}
            fi
        elif [[ ( (${read1_file} =~ RNA) || ( (${read1_file} =~ NDNR) && (${NDNR_AS} == "RNA" || ${NDNR_AS} == "rna" ) ) ) && ( (${rna_aligner} == "star") || (${rna_aligner} == "STAR") ) ]]
        then
            star_map "RNA"
        fi
    fi
}
done

if [[ $FinishSuccess == 1 ]]
then
    cd ..
    echo "Finished sequence alignment."
    echo "$0 done" > SequenceAlignment.done
else
    echo "Sequence alignment failed."
fi

