#!/bin/bash

cnt=1
for i in "$@"
do
    if [[ $cnt == 1 ]]; then LIB_NAME=$i; fi
    if [[ $cnt == 2 ]]; then rdSPRITE_Read_1=$i; fi
    if [[ $cnt == 3 ]]; then rdSPRITE_Read_2=$i; fi
    if [[ $cnt == 4 ]]; then THREAD=$i; fi
    if [[ $cnt == 5 ]]; then CUSTOM_CONFIG_FILE=$i; fi
    if [[ $cnt == 6 ]]; then STAR_REF_GENOME=$i; fi
    if [[ $cnt == 7 ]]; then BWA_REF_GENOME=$i; fi
    if [[ $cnt == 8 ]]; then TBARP_DIR=$i; fi
    (( cnt++ ))
done

if [[ ! -z $TBARP_DIR ]]
then
    BARP_DIR=$TBARP_DIR
else
    BARP_DIR="REF_OF_ScSmOP"
fi

echo "Start: $(date)" > TimeStamp
if [[ ${BARP_DIR} == "-" ]]
then
    echo -e "BARP directory is necessary, but not exist."
    exit 1
fi

${BARP_DIR}/Tools/trim_galore --paired --gzip --cores ${THREAD} --quality 20 --fastqc ${rdSPRITE_Read_1} ${rdSPRITE_Read_2} --basename ${LIB_NAME} --path_to_cutadapt ${BARP_DIR}/Tools/cutadapt

${BARP_DIR}/PipelineScript/BarcodeIdentification.sh -t rdsprite -n ${LIB_NAME} -p ${BARP_DIR} -1 ${LIB_NAME}_val_1.fq.gz -2 ${LIB_NAME}_val_2.fq.gz -@ ${THREAD} -c ${CUSTOM_CONFIG_FILE}
if [[ -f BarcodeIdentification.done ]]
then
    ${BARP_DIR}/PipelineScript/SequenceAlignment.sh -n ${LIB_NAME} -p ${BARP_DIR} -f ${BWA_REF_GENOME} -@ ${THREAD} -g ${STAR_REF_GENOME}
fi

if [[ -f SequenceAlignment.done ]]
then
    ${BARP_DIR}/PipelineScript/GroupAndDataRefine.sh -t rdsprite -n ${LIB_NAME} -p ${BARP_DIR} -d 02.ReadAlign/${LIB_NAME}_DNA.bam -@ ${THREAD} -r 02.ReadAlign/${LIB_NAME}Aligned.sortedByCoord.out.bam
fi

if [[ -f GroupAndDataRefine.done ]]
then
    ${BARP_DIR}/PipelineScript/QualityAssessment.sh -t rdsprite -n ${LIB_NAME} -p ${BARP_DIR}
fi

echo "End: $(date)" >> TimeStamp