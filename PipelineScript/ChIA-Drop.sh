#!/bin/bash

cnt=1
for i in "$@"
do
    if [[ $cnt == 1 ]]; then LIB_NAME=$i; fi
    if [[ $cnt == 2 ]]; then ChIADrop_FILES=$i; fi
    if [[ $cnt == 3 ]]; then THREAD=$i; fi
    if [[ $cnt == 4 ]]; then ChIADrop_CONFIG=$i; fi
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

echo -e "${BARP_DIR}/PipelineScript/BarcodeIdentification.sh -t chiadrop -n ${LIB_NAME} -p ${BARP_DIR} ${ChIADrop_FILES} -@ ${THREAD} -c ${ChIADrop_CONFIG}"

${BARP_DIR}/PipelineScript/BarcodeIdentification.sh -t chiadrop -n ${LIB_NAME} -p ${BARP_DIR} ${ChIADrop_FILES} -@ ${THREAD} -c ${ChIADrop_CONFIG}

if [[ -f BarcodeIdentification.done ]]
then
    ${BARP_DIR}/PipelineScript/SequenceAlignment.sh -t chiadrop -n ${LIB_NAME} -p ${BARP_DIR} -f ${BWA_REF_GENOME} -@ ${THREAD} -k "DNA"
fi

if [[ -f SequenceAlignment.done ]]
then
    ${BARP_DIR}/PipelineScript/GroupAndDataRefine.sh -t chiadrop -n ${LIB_NAME} -p ${BARP_DIR} -d 02.ReadAlign/${LIB_NAME}_DNA.bam -@ ${THREAD} -x ${GENOME_SIZE_FILE}
fi

if [[ -f GroupAndDataRefine.done ]]
then
    ${BARP_DIR}/PipelineScript/QualityAssessment.sh -t chiadrop -n ${LIB_NAME} -p ${BARP_DIR}
fi

echo "End: $(date)" >> TimeStamp