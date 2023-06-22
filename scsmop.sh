#!/bin/bash

print_help()
{
    echo ''
	echo -e "Usage: $0 [options] -t [STR] -n [STR] -1 [FILE] -2 [FILE] -r [DIR]"
    echo -e "\t-t Experiment type: sprite, scsprite, rdsprite, chiadrop, scrna_10x_v3, scrna_10x_v2, dropseq, scatac_10x_v1, scarc_10x_v1."
	echo -e "\t-n Library name."
    echo -e "\t-1 Read 1 FASTQ."
    echo -e "\t-2 Read 2 FASTQ."
    echo -e "\t-3 Read 3 FASTQ."
    echo -e "\t-4 Read 4 FASTQ."
    echo -e "\t-5 Read 5 FASTQ."
    echo -e "\t-r STAR reference genome directory."
    echo -e "\t-b BWA reference genome directory."
    echo -e "\t-c Custom configuration file."
    echo -e "\t-w Custom whitelist."
	echo -e "\t-@ Thread to use."
    echo -e "\t-p Custom ScSmOP directory."
    echo -e "\t-s Chromosome size file, need by chiadrop, scatac_10x_v1, scarc_10x_v1."
    echo -e "\t-s Processing scARC-seq, specify ATAC FASTQ through -1 ATAC R1 -2 ATAC R2 -3 ATAC R3 -4 GEX R1 -5 GEX R2.\n"
    echo -e "Universal pipeline for multi-omics data process: <https://github.com/ZhengmzLab/ScSmOP/wiki>.\n"
}

if [[ $1 == '' ]]
then
    print_help
    exit 1
fi

thread=1
read_1_str="-"
read_2_str="-"
read_3_str="-"
read_4_str="-"
custom_config="-"
custom_whitelist="-"
chrom_size="-"

barp_dir="REF_OF_ScSmOP"

while getopts t:n:1:2:3:4:5:r:b:c:w:@:p:s:h flag
do
    case "${flag}" in 
        t)  expe_type=${OPTARG};;
        n)  name=${OPTARG};;
        1)  read_1_str=${OPTARG};;
        2)  read_2_str=${OPTARG};;
        3)  read_3_str=${OPTARG};;
        4)  read_4_str=${OPTARG};;
        5)  read_5_str=${OPTARG};;
        r)  star_ref=${OPTARG};;
        b)  bwa_ref=${OPTARG};;
        c)  custom_config=${OPTARG};;
        w)  custom_whitelist=${OPTARG};;
        @)  thread=${OPTARG};;
        p)  barp_dir=${OPTARG};;
        s)  chrom_size=${OPTARG};;
        ? | h) print_help
            exit 1;;
    esac
done

if [[ -z ${expe_type} ]]
then
    echo -e "\nCurrently, easy process only support for the following experiment type with/without subsititued whitelist. For more complicated DIY process, please reference to Github: <https://github.com/ZhengmzLab/ScSmOP>."
    print_help
    exit 1
fi

if [[ -z ${name} ]]
then
    name=$( echo "${expe_type}" | tr [:lower:] [:upper:] )
fi

function checkBWARef
{
    if [[ ! -f ${bwa_ref} ]]
    then
        echo -e "\nBWA reference genome is necessary for $@ data process, but do not provided or exist, please check..."
        print_help
        exit 1
    fi
}

function checkSTARRef 
{
    if [[ ! -d ${star_ref} ]]
    then
        echo -e "\nSTAR reference genome is necessary for $@ data process, but do not provided or exist, please check..."
        print_help
        exit 1
    fi
}

function checkChromSize
{
    if [[ ! -f ${chrom_size} ]]
    then
        echo -e "\nChromosome size file is necessary for $@ data process, but do not provided or exist, please check..."
        print_help
        exit 1
    fi
}

if [[ -f ${custom_config} ]]
then 
    custom_config=$( readlink -f ${custom_config} )
fi

if [[ ${expe_type} == "sprite" ]]
then
    if [[ ${custom_whitelist} != "-" ]]; then echo "SPRITE do not support custom whitelist in this version, please refer to DIY."; fi
    checkBWARef ${expe_type}
    ${barp_dir}/PipelineScript/SPRITE.sh ${name} ${read_1_str} ${read_2_str} ${thread} ${custom_config} ${bwa_ref} ${barp_dir}
elif [[ ${expe_type} == "scsprite" ]]
then
    if [[ ${custom_whitelist} != "-" ]]; then echo "scSPRITE do not support custom whitelist in this version, please refer to DIY."; fi
    SPRITE_FILES=$( echo "-1 ${read_1_str} -2 ${read_2_str}" )
    checkSTARRef ${expe_type}
    ${barp_dir}/PipelineScript/scSPRITE.sh ${name} "${SPRITE_FILES}" ${thread} ${custom_config} ${star_ref} ${barp_dir}
elif [[ ${expe_type} == "rdsprite" ]]
then
    if [[ ${custom_whitelist} != "-" ]]; then echo "rdSPRITE do not support custom whitelist in this version, please refer to DIY."; fi
    checkSTARRef ${expe_type}
    checkBWARef ${expe_type}
    ${barp_dir}/PipelineScript/rdSPRITE.sh ${name} ${read_1_str} ${read_2_str} ${thread} ${custom_config} ${star_ref} ${bwa_ref} ${barp_dir}
elif [[ ${expe_type} == "chiadrop" ]]
then
    if [[ ${custom_whitelist} != "-" ]];
    then
        ${barp_dir}/Tools/python3 ${barp_dir}/PythonScript/PipeGenerateConfig.py -i ${barp_dir}/ConfigFiles/.OriginalConfigFile_ChIA-Drop.json -o ${name} -w ${custom_whitelist}
        custom_config=$( readlink -f ${name}_config.json )
    fi
    ChIADrop_FILES=$( echo "-1 ${read_1_str} -2 ${read_2_str}" )
    checkBWARef ${expe_type}
    checkChromSize ${expe_type}
    echo -e "${barp_dir}/PipelineScript/ChIA-Drop.sh ${name} "${ChIADrop_FILES}" ${thread} ${custom_config} ${bwa_ref} ${chrom_size} ${barp_dir}"
    ${barp_dir}/PipelineScript/ChIA-Drop.sh ${name} "${ChIADrop_FILES}" ${thread} ${custom_config} ${bwa_ref} ${chrom_size} ${barp_dir}
elif [[ ${expe_type} == "scrna_10x_v3" ]]
then   
    if [[ ${custom_whitelist} != "-" ]];
    then
        ${barp_dir}/Tools/python3 ${barp_dir}/PythonScript/PipeGenerateConfig.py -i ${barp_dir}/ConfigFiles/.OriginalConfigFile_scRNA_10x.json -o ${name} -w ${custom_whitelist}
        custom_config=$( readlink -f ${name}_config.json )
    fi
    Input_FILES=$( echo "-1 ${read_1_str} -2 ${read_2_str}" )
    checkSTARRef ${expe_type}
    ${barp_dir}/PipelineScript/scRNA.sh ${name} "${Input_FILES}" ${thread} ${custom_config} ${star_ref} ${barp_dir}
elif [[ ${expe_type} == "scrna_10x_v2" ]]
then
    if [[ ${custom_whitelist} != "-" ]];
    then
        ${barp_dir}/Tools/python3 ${barp_dir}/PythonScript/PipeGenerateConfig.py -i ${barp_dir}/ConfigFiles/.OriginalConfigFile_scRNA_10x.json -o ${name} -w ${custom_whitelist}
        custom_config=$( readlink -f ${name}_config.json )
    fi
    Input_FILES=$( echo "-1 ${read_1_str} -2 ${read_2_str}" )
    checkSTARRef ${expe_type}
    ${barp_dir}/PipelineScript/10x_scrna-v2.sh ${name} "${Input_FILES}" ${thread} ${custom_config} ${star_ref} ${barp_dir}
elif [[ ${expe_type} == "dropseq" ]]
then
    if [[ ${custom_whitelist} != "-" ]];
    then
        ${barp_dir}/Tools/python3 ${barp_dir}/PythonScript/PipeGenerateConfig.py -i ${barp_dir}/ConfigFiles/.OriginalConfigFile_scRNA_Dropseq.json -o ${name} -w ${custom_whitelist}
        custom_config=$( readlink -f ${name}_config.json )
    fi
    Input_FILES=$( echo "-1 ${read_1_str} -2 ${read_2_str}" )
    checkSTARRef ${expe_type}
    ${barp_dir}/PipelineScript/dropseq_scrna.sh ${name} "${Input_FILES}" ${thread} ${custom_config} ${star_ref} ${barp_dir}
elif [[ ${expe_type} == "scatac_10x_v1" ]]
then
    if [[ ${custom_whitelist} != "-" ]];
    then
        ${barp_dir}/Tools/python3 ${barp_dir}/PythonScript/PipeGenerateConfig.py -i ${barp_dir}/ConfigFiles/.OriginalConfigFile_scATAC_10x.json -o ${name} -w ${custom_whitelist}
        custom_config=$( readlink -f ${name}_config.json )
    fi
    Input_FILES=$( echo "-1 ${read_1_str} -2 ${read_2_str} -3 ${read_3_str}" )
    checkChromSize ${expe_type}
    ${barp_dir}/PipelineScript/scATAC.sh ${name} "${Input_FILES}" ${thread} ${custom_config} ${bwa_ref} ${chrom_size} ${barp_dir}
elif [[ ${expe_type} == "scarc_10x_v1" ]]
then
    if [[ ${custom_whitelist} != "-" ]]; then echo "scARC-seq do not support custom whitelist in this version, please refer to DIY."; fi
    checkSTARRef ${expe_type}
    checkBWARef ${expe_type}
    checkChromSize ${expe_type}
    ATAC_Input_FILES=$( echo "-1 ${read_1_str} -2 ${read_2_str} -3 ${read_3_str}" )
    GEX_Input_FILES=$( echo "-1 ${read_4_str} -2 ${read_5_str}" )
    atac_config="10x_scarc-atac_config.json"
    gex_config="10x_scarc-rna_config.json"
    checkChromSize ${expe_type}
    ${barp_dir}/PipelineScript/scarc.sh ${name} "${ATAC_Input_FILES}" "${GEX_Input_FILES}" ${thread} ${atac_config} ${gex_config} ${star_ref} ${bwa_ref} ${chrom_size} ${barp_dir}
elif [[ ${expe_type} == "scrna" ]]
then
    echo -e "\nYou need to be more specific about the experiment to \"scrna_10x_v2\", \"scrna_10_v3\" or \"dropseq\"."
    print_help
    exit 1
elif [[ ${expe_type} == "scatac" ]]
then
    if [[ ${custom_whitelist} != "-" ]];
    then
        ${barp_dir}/Tools/python3 ${barp_dir}/PythonScript/PipeGenerateConfig.py -i ${barp_dir}/ConfigFiles/.OriginalConfigFile_scATAC_10x.json -o ${name} -w ${custom_whitelist}
        custom_config=$( readlink -f ${name}_config.json )
    fi
    Input_FILES=$( echo "-1 ${read_1_str} -2 ${read_2_str} -3 ${read_3_str}" )
    checkChromSize ${expe_type}
    ${barp_dir}/PipelineScript/scATAC.sh ${name} "${Input_FILES}" ${thread} ${custom_config} ${bwa_ref} ${chrom_size} ${barp_dir}
    echo -e "\nSetting experiment type to \"scatac_10x_v1\" ...\n"
    expe_type="scatac_10x_v1"
else
    echo -e "\nCurrently, easy process only support for the following experiment type with/without subsititued whitelist. For more complicated DIY process, please reference to Github: <https://github.com/ZhengmzLab/ScSmOP>."
    print_help
    exit 1
fi

