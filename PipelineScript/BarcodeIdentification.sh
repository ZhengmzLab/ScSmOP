#!/bin/bash
start=$SECOND

print_help()
{
    echo ''
	echo -e "Usage: $0 -t [STR] -n [STR] -p [DIR] -c [Config] -1 [STR] -2 [STR] -3 [STR] -4 [STR] -@ [INT] -a -m [STR]"
    echo -e "\t-t Experiment type: sprite, scsprite, rdsprite, chiadrop, scrna, scatac, scrnawscatac"
	echo -e "\t-n Library name"
    echo -e "\t-p Pipeline processing tool directory"
    echo -e "\t-c Config file"
	echo -e "\t-1 List of Read 1 files, separate with ",""
    echo -e "\t-2 List of Read 2 files, separate with ",""
	echo -e "\t-@ Thread to use"
	echo -e "\t-a Retain all reads instead of fully barcoded reads"
	echo -e "\t-m Mode"
}
if [[ $1 == '' ]]
then
    print_help
    exit 1
fi
num_thread=1
config="-"
name="-"
mode="-"
read_1_list="-"
read_2_list="-"
read_3_list="-"
read_4_list="-"

while getopts t:p:n:c:f:@:m:a1:2:3:4:h flag
do
    case "${flag}" in 
        t)  expe_type=${OPTARG};;
        p)  pipe_dir=${OPTARG};;
        n)  name=${OPTARG};;
        c)  config=${OPTARG};;
        1)  read_1_list=${OPTARG};;
        2)  read_2_list=${OPTARG};;
        3)  read_3_list=${OPTARG};;
        4)  read_4_list=${OPTARG};;
        @)  num_thread=${OPTARG};;
        a)  all_read=1;;
        m)  mode=${OPTARG};;
        ? | h) print_help
            exit 1;;
    esac
done

# experiment type can be sprite, chiadrop, scrna, scatac, scrnawscatac

process_date=$(date +%Y%b%d | tr [:lower:] [:upper:])
run_dir=$(pwd)

if [[ (${read_1_list} == "-") && (${read_2_list} == "-") && (${read_3_list} == "-") && (${read_4_list} == "-") ]]
then
    echo "No input file specified."
    exit 1
fi


if [[ ${read_1_list} != "-" ]]
then
    read_1_param="-1 "
    r1_file_list=(${read_1_list//,/ })
    for((i=0;i<${#r1_file_list[@]};i++))
    do
        F=${r1_file_list[i]}
        if [[ ! $F =~ ^/ ]]
        then
            if [[ $i == 0 ]]
            then
                read_1_param=$(echo -e "${read_1_param}../${F}")
            else
                echo "I = $i Read 1 ${read_1_param}"
                read_1_param=$(echo -e "${read_1_param},../${F}")
            fi
        else
            if [[ $i == 0 ]]
            then
                read_1_param=$(echo -e "${read_1_param}${F}")
            else
                read_1_param=$(echo -e "${read_1_param},${F}")
            fi
        fi
    done
fi

if [[ ${read_2_list} != "-" ]]
then
    read_2_param="-2 "
    r2_file_list=(${read_2_list//,/ })
    for((i=0;i<${#r2_file_list[@]};i++))
    do
        F=${r2_file_list[$i]}
        if [[ ! $F =~ ^/ ]]
        then
            if [[ $i == 0 ]]
            then
                read_2_param=$(echo -e "${read_2_param}../${F}")
            else
                read_2_param=$(echo -e "${read_2_param},../${F}")
            fi
        else
            if [[ $i == 0 ]]
            then
                read_2_param=$(echo -e "${read_2_param}${F}")
            else
                read_2_param=$(echo -e "${read_2_param},${F}")
            fi
        fi
    done
fi

if [[ ${read_3_list} != "-" ]]
then
    read_3_param="-3 "
    r3_file_list=(${read_3_list//,/ })
    for((i=0;i<${#r3_file_list[@]};i++))
    do
        F=${r3_file_list[$i]}
        if [[ ! $F =~ ^/ ]]
        then
            if [[ $i == 0 ]]
            then
                read_3_param=$(echo -e "${read_3_param}../${F}")
            else
                read_3_param=$(echo -e "${read_3_param},../${F}")
            fi
        else
            if [[ $i == 0 ]]
            then
                read_3_param=$(echo -e "${read_3_param}${F}")
            else
                read_3_param=$(echo -e "${read_3_param},${F}")
            fi
        fi
    done
fi

if [[ ${read_4_list} != "-" ]]
then
    read_4_param="-4 "
    r4_file_list=(${read_4_list//,/})
    for((i=0;i<${#r4_file_list[@]};i++))
    do
        F=${r4_file_list[$i]}
        if [[ ! $F =~ ^/ ]]
        then
            if [[ $i == 0 ]]
            then
                read_4_param=$(echo -e "${read_4_param}../${F}")
            else
                read_4_param=$(echo -e "${read_4_param},../${F}")
            fi
        else
            if [[ $i == 0 ]]
            then
                read_4_param=$(echo -e "${read_4_param}${F}")
            else
                read_4_param=$(echo -e "${read_4_param},${F}")
            fi
        fi
    done
fi


if [[ ! -d 01.BarcodeIden ]]
then
    mkdir 01.BarcodeIden
    cd 01.BarcodeIden
else 
    echo -e "\nDirectory 01.BarcodeIden exist, will cover the files generated before in 5 seconds...";
    sleep 5
    cd 01.BarcodeIden
    rm * 2>/dev/null
fi
echo -e "\nBarcode identification..."
echo -e "Working directory: $(pwd)"
echo ""
if [[ ${config} == "-" ]]
then
    if [[ ${expe_type} == "sprite" ]]
    then
        ori_config=${pipe_dir}/ConfigFiles/sprite_config.json
    elif [[ ${expe_type} == "chiadrop" ]]
    then
        ori_config=${pipe_dir}/ConfigFiles/chia-drop_config.json
    elif [[ ${expe_type} == "dropseq" ]]
    then
        ori_config=${pipe_dir}/ConfigFiles/dropseq_scrna_config.json
    elif [[ ${expe_type} == "scrna_10x_v3" ]]
    then
        ori_config=${pipe_dir}/ConfigFiles/10x_scrna_v3_config.json
    elif [[ ${expe_type} == "scrna_10x_v2" ]]
    then
        ori_config=${pipe_dir}/ConfigFiles/10x_scrna_v2_config.json
    elif [[ ${expe_type} == "scatac" ]]
    then
        ori_config=${pipe_dir}/ConfigFiles/10x_scatac_v2_config.json
    elif [[ ${expe_type} == "rdsprite" ]]
    then
        ori_config=${pipe_dir}/ConfigFiles/rdsprite_config.json
    elif [[ ${expe_type} == "scsprite" ]]
    then
        ori_config=${pipe_dir}/ConfigFiles/scsprite_config.json
    else
        echo -e "Neither configuration file nor experiment type provided, provide experiment type with -t or configuration file with -c "
        exit 1
    fi
else
    if [[ -f ${config} ]]
    then
        ori_config=${config}
    else
        ori_config=${pipe_dir}/ConfigFiles/${config}
    fi
fi


if [[ ${ori_config} == "-" ]]
then
    echo "No config file specified."
    exit 1
fi


if [[ ${expe_type} == "chiadrop" || ${expe_type} == "scrna" || ${expe_type} == "scatac" ]]
then
    BarcodeMode="atact"
else
    BarcodeMode="default"
fi

if [[ ${mode} != "-" ]]
then
    BarcodeMode=${mode}
fi



if [[ ${BarcodeMode} != "-" ]]
then
    echo -e "${pipe_dir}/Tools/BARP idb ${read_1_param} ${read_2_param} ${read_3_param} ${read_4_param} -c ${ori_config} -@ ${num_thread} -o ${name} --mode ${BarcodeMode}"
    ${pipe_dir}/Tools/BARP idb ${read_1_param} ${read_2_param} ${read_3_param} ${read_4_param} -c ${ori_config} -@ ${num_thread} -o ${name} --mode ${BarcodeMode} 
else
    echo -e "${pipe_dir}/Tools/BARP idb ${read_1_param} ${read_2_param} ${read_3_param} ${read_4_param} -c ${ori_config} -@ ${num_thread} -o ${name} --mode default"
    ${pipe_dir}/Tools/BARP idb ${read_1_param} ${read_2_param} ${read_3_param} ${read_4_param} -c ${ori_config} -@ ${num_thread} -o ${name} 
fi

if [[ $? != 0 ]]
then
    echo "Barcode identification failed."
else
    cd ..
    echo "Finished barcode identification."
    echo "$0 done" > BarcodeIdentification.done
fi


