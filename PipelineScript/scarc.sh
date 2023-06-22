#!/bin/bash

cnt=1
for i in "$@"
do
    if [[ $cnt == 1 ]]; then LIB_NAME=$i; fi
    if [[ $cnt == 2 ]]; then ATAC_FILES=$i; fi
    if [[ $cnt == 3 ]]; then RNA_FILES=$i; fi
    if [[ $cnt == 4 ]]; then THREAD=$i; fi
    if [[ $cnt == 5 ]]; then ATAC_CONFIG=$i; fi
    if [[ $cnt == 6 ]]; then RNA_CONFIG=$i; fi
    if [[ $cnt == 7 ]]; then STAR_REF_GENOME=$i; fi
    if [[ $cnt == 8 ]]; then BWA_REF_GENOME=$i; fi
    if [[ $cnt == 9 ]]; then GENOME_SIZE_FILE=$i; fi
    if [[ $cnt == 10 ]]; then TBARP_DIR=$i; fi
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


${BARP_DIR}/PipelineScript/BarcodeIdentification.sh -t scrna_10x_v3 -n ${LIB_NAME} -p ${BARP_DIR} ${RNA_FILES} -@ ${THREAD} -c ${RNA_CONFIG}

if [[ -f BarcodeIdentification.done ]]
then
    ${BARP_DIR}/PipelineScript/SequenceAlignment.sh -t scrna -n ${LIB_NAME} -p ${BARP_DIR} -g ${STAR_REF_GENOME} -@ ${THREAD} -k "RNA" -s
fi

if [[ -f SequenceAlignment.done ]]
then
    ${BARP_DIR}/PipelineScript/GroupAndDataRefine.sh -t scrna -n ${LIB_NAME} -p ${BARP_DIR} -r 02.ReadAlign/${LIB_NAME}Aligned.sortedByCoord.out.bam -c ${RNA_CONFIG} -g ${STAR_REF_GENOME} -@ ${THREAD}
fi

echo "GroupAndDataRefine done" > GroupAndDataRefine.done

if [[ -f GroupAndDataRefine.done ]]
then
    ${BARP_DIR}/PipelineScript/QualityAssessment.sh -t scrna -n ${LIB_NAME} -p ${BARP_DIR}
fi


if [[ ! -d RNAResult ]]
then
    mkdir RNAResult
fi

mv 01.BarcodeIden 02.ReadAlign 04.QualityAssess *.done RNAResult


${BARP_DIR}/PipelineScript/BarcodeIdentification.sh -t scatac -n ${LIB_NAME} -p ${BARP_DIR} ${ATAC_FILES} -@ ${THREAD} -c ${ATAC_CONFIG}
if [[ -f BarcodeIdentification.done ]]
then
   echo -e "${BARP_DIR}/PipelineScript/SequenceAlignment.sh -t scatac -n ${LIB_NAME} -p ${BARP_DIR} -f ${BWS_REF_GENOME} -@ ${THREAD} -k "DNA" -1 01.BarcodeIden/${LIB_NAME}_NDNR_1.fq -2 01.BarcodeIden/${LIB_NAME}_NDNR_3.fq"
   ${BARP_DIR}/PipelineScript/SequenceAlignment.sh -t scatac -n ${LIB_NAME} -p ${BARP_DIR} -f ${BWA_REF_GENOME} -@ ${THREAD} -k "DNA" -1 01.BarcodeIden/${LIB_NAME}_NDNR_1.fq -2 01.BarcodeIden/${LIB_NAME}_NDNR_3.fq
fi

if [[ -f SequenceAlignment.done ]]
then
   ${BARP_DIR}/PipelineScript/GroupAndDataRefine.sh -t arcatac -n ${LIB_NAME} -p ${BARP_DIR} -d 02.ReadAlign/${LIB_NAME}_DNA.bam -@ ${THREAD} -x ${GENOME_SIZE_FILE}
fi

if [[ -f GroupAndDataRefine.done ]]
then
    ${BARP_DIR}/PipelineScript/QualityAssessment.sh -t scatac -n ${LIB_NAME} -p ${BARP_DIR}
fi

if [[ ! -d ATACResult ]]
then
   mkdir ATACResult
fi

mv 01.BarcodeIden 02.ReadAlign 03.GroupAndRefine 04.QualityAssess *.done ATACResult
echo "End: $(date)" >> TimeStamp