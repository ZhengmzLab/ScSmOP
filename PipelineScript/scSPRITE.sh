#!/bin/bash

cnt=1
for i in "$@"
do
    if [[ $cnt == 1 ]]; then LIB_NAME=$i; fi
    if [[ $cnt == 2 ]]; then SPRITE_FILES=$i; fi
    if [[ $cnt == 3 ]]; then THREAD=$i; fi
    if [[ $cnt == 4 ]]; then SPRITE_CONFIG=$i; fi
    if [[ $cnt == 5 ]]; then STAR_REF_GENOME=$i; fi
    if [[ $cnt == 6 ]]; then TBARP_DIR=$i; fi
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


${BARP_DIR}/PipelineScript/BarcodeIdentification.sh -t scsprite -n ${LIB_NAME} -p ${BARP_DIR} ${SPRITE_FILES} -@ ${THREAD} -c ${SPRITE_CONFIG} 
if [[ -f BarcodeIdentification.done ]]
then
    ${BARP_DIR}/PipelineScript/SequenceAlignment.sh -t scsprite -n ${LIB_NAME} -p ${BARP_DIR} -g ${STAR_REF_GENOME} -@ ${THREAD} -k "DNA" -d star
fi

if [[ -f SequenceAlignment.done ]]
then
    ${BARP_DIR}/PipelineScript/GroupAndDataRefine.sh -t scsprite -n ${LIB_NAME} -p ${BARP_DIR} -d 02.ReadAlign/${LIB_NAME}.Barcoded.Aligned.out.bam
fi

if [[ -f GroupAndDataRefine.done ]]
then
    ${BARP_DIR}/PipelineScript/QualityAssessment.sh -t scsprite -n ${LIB_NAME} -p ${BARP_DIR} 
fi

echo "End: $(date)" >> TimeStamp
