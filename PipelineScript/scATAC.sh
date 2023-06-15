#!/bin/bash

cnt=1
for i in "$@"
do
    if [[ $cnt == 1 ]]; then LIB_NAME=$i; fi
    if [[ $cnt == 2 ]]; then ATAC_FILES=$i; fi
    if [[ $cnt == 3 ]]; then THREAD=$i; fi
    if [[ $cnt == 4 ]]; then ATAC_CONFIG=$i; fi
    if [[ $cnt == 5 ]]; then BWA_REF_GENOME=$i; fi
    if [[ $cnt == 6 ]]; then GENOME_SIZE_FILE=$i; fi
    if [[ $cnt == 7 ]]; then TBARP_DIR=$i; fi
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


${BARP_DIR}/PipelineScript/BarcodeIdentification.sh -t scatac -n ${LIB_NAME} -p ${BARP_DIR} ${ATAC_FILES} -@ ${THREAD} -c ${ATAC_CONFIG}

if [[ -f BarcodeIdentification.done ]]
then
    ${BARP_DIR}/PipelineScript/SequenceAlignment.sh -t scatac -n ${LIB_NAME} -p ${BARP_DIR} -f ${BWA_REF_GENOME} -@ ${THREAD} -k "DNA" -1 01.BarcodeIden/${LIB_NAME}_NDNR_1.fq -2 01.BarcodeIden/${LIB_NAME}_NDNR_3.fq
fi

if [[ -f SequenceAlignment.done ]]
then
    ${BARP_DIR}/PipelineScript/GroupAndDataRefine.sh -t scatac -n ${LIB_NAME} -p ${BARP_DIR} -d 02.ReadAlign/${LIB_NAME}_DNA.bam -@ ${THREAD} -x ${GENOME_SIZE_FILE}
fi

if [[ -f GroupAndDataRefine.done ]]
then
    ${BARP_DIR}/PipelineScript/QualityAssessment.sh -t scatac -n ${LIB_NAME} -p ${BARP_DIR}
fi

echo "End: $(date)" >> TimeStamp
