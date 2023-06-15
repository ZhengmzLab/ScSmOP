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

run_dir=$(pwd)
name="-"
num_thread=1
config="-"

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

FinishSuccess=1

if [[ ! -d 03.GroupAndRefine ]]
then
    mkdir 03.GroupAndRefine
    cd 03.GroupAndRefine
else 
    echo -e "\nDirectory 03.GroupAndRefine exist, will cover the files generated before in 5 seconds...";
    sleep 5
    cd 03.GroupAndRefine
fi

dna_alignment="../${dna_alignment}"
rna_alignment="../${rna_alignment}"

echo -e "\nGrouping and data refine..."
echo -e "Working directory: $(pwd)"
echo ""

if [[ ! ${dna_alignment} =~ ^/ ]]
then
    dna_alignment=${dna_alignment}
fi



if [[ ${expe_type} == "scsprite" ]]
then
    echo -e "DNA alignment ${dna_alignment}"
    ${pipe_dir}/Tools/samtools view -@ ${num_thread} -h -b -q 255 ${dna_alignment} -o ${name}.Barcoded.Aligned.out.bam
    ${pipe_dir}/Tools/samtools stats -@ ${num_thread} ${name}.Barcoded.Aligned.out.bam > ${name}.Barcoded.UniqAlign.stat
    ${pipe_dir}/Tools/BARP p2s -b "CELL-COMPLEX" -s "CELL:COMPLEX" ${name}.Barcoded.Aligned.out.bam -o ${name}.DNA
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/RemoveDuplication.py -i ${name}.DNA.cluster -r "CHR-START" -o ${name}.DNA
    ${pipe_dir}/Tools/python3  ${pipe_dir}/PythonScript/GetFrag2FragDistance.py ${name}.DNA.RDP.cluster 3
    # sort -k2 -k6 ${name}.DNA.RDP.cluster.juice > ${name}.DNA.RDP.cluster.sorted.juice
    # awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$9,$7}' ${name}.DNA.RDP.cluster.L2L.stat > ${name}.DNA.loop
    awk '{OFS="\t"; print $1,$2,$3,$5}' ${name}.DNA.RDP.cluster.FragLen.stat > ${name}.DNA.rgn
    # java -jar ${pipe_dir}/Tools/juicer_tools_1.19.01.jar pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 ${name}.DNA.RDP.cluster.sorted.juice ${name}.hic ${pipe_dir}/ChromSize/${genomesize}
    echo -e "SubGEMing..."
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/ChIADropSubGEM.py ${name}.DNA.RDP.cluster 3
    echo -e "Generating region file..."
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GetFrag2FragDistance.py ${name}.DNA.RDP.cluster.SUBGEM 3 
    awk '{OFS="\t"; print $1,$2,$3,$4}' ${name}.DNA.RDP.cluster.SUBGEM.FragLen.stat > ${name}.DNA.SUBGEM.rgn
    if [[ ! -d StatFiles ]]
    then
        mkdir StatFiles
    fi
    mv *stat StatFiles
    rm ${name}.DNA_P2S.bam
    if [[ $? != 0 ]]
    then
        FinishSuccess=0
    fi
elif [[ (${expe_type} == "rdsprite") ]]
then 
    echo -e "DNA alignment ${dna_alignment}"
    ${pipe_dir}/Tools/samtools view -@ ${num_thread} -h -b -q 30 -o ${name}.Barcoded.UniqAlign.DNA.bam -F 2304 ${dna_alignment}
    ${pipe_dir}/Tools/samtools stats -@ ${num_thread} ${name}.Barcoded.UniqAlign.DNA.bam > ${name}.Barcoded.UniqAlign.DNA.stat
    ${pipe_dir}/Tools/BARP p2s -b "COMPLEX" ${name}.Barcoded.UniqAlign.DNA.bam -o ${name}.DNA
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/RemoveDuplication.py -i ${name}.DNA.cluster -r "CHR-START" -o ${name}.DNA
    ${pipe_dir}/Tools/python3  ${pipe_dir}/PythonScript/GetFrag2FragDistance.py ${name}.DNA.RDP.cluster 2
    awk '{OFS="\t"; print $1,$2,$3,$4}' ${name}.DNA.RDP.cluster.FragLen.stat > ${name}.DNA.rgn
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/ChIADropSubGEM.py ${name}.DNA.RDP.cluster 2
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GetFrag2FragDistance.py ${name}.DNA.RDP.cluster.SUBGEM 2 
    awk '{OFS="\t"; print $1,$2,$3,$4}' ${name}.DNA.RDP.cluster.SUBGEM.FragLen.stat > ${name}.DNA.SUBGEM.rgn
    
    # awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$9,$7}' L2L.stat > ${name}.DNA.loop
    # sort -k2 -k6 ${name}.DNA.RDP.cluster.SUBGEM.juice > ${name}.DNA.RDP.cluster.SUBGEM.sorted.juice
    # java -jar ${pipe_dir}/Tools/juicer_tools_1.19.01.jar pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 ${name}.RDP.DNA.cluster.SUBGEM.sorted.juice ${name}.hic ${pipe_dir}/ChromSize/${genomesize}

    echo -e "RNA alignment ${rna_alignment}"
    ${pipe_dir}/Tools/samtools view -@ ${num_thread} -h -b -q 255 ${rna_alignment} -o ${name}.Barcoded.UniqAlign.RNA.bam
    ${pipe_dir}/Tools/samtools stats -@ ${num_thread} ${name}.Barcoded.UniqAlign.RNA.bam > ${name}.Barcoded.UniqAlign.RNA.stat
    ${pipe_dir}/Tools/BARP p2s -b "COMPLEX" ${name}.Barcoded.UniqAlign.RNA.bam -o ${name}.RNA
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/RemoveDuplication.py -i ${name}.RNA.cluster -r "CHR-START" -o ${name}.RNA
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GetFrag2FragDistance.py ${name}.RNA.RDP.cluster 2
    awk '{OFS="\t"; print $1,$2,$3,$4}' ${name}.RNA.RDP.cluster.FragLen.stat > ${name}.RNA.rgn
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/ChIADropSubGEM.py ${name}.RNA.RDP.cluster 2
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GetFrag2FragDistance.py ${name}.RNA.RDP.cluster.SUBGEM 2 
    awk '{OFS="\t"; print $1,$2,$3,$4}' ${name}.RNA.RDP.cluster.SUBGEM.FragLen.stat > ${name}.RNA.SUBGEM.rgn
    awk '{OFS="\t"; print $0,"DNA"}' ${name}.DNA.rgn > ${name}.Mixed.rgn
    awk '{OFS="\t"; print $0,"RNA"}' ${name}.RNA.rgn >> ${name}.Mixed.rgn
    if [[ ! -d StatFiles ]]
    then
        mkdir StatFiles
    fi
    mv *stat StatFiles
    rm ${name}.DNA_P2S.bam
    rm ${name}.RNA_P2S.bam
    if [[ $? != 0 ]]
    then
        FinishSuccess=0
    fi
elif [[ (${expe_type} == "sprite") ]]
then 
    echo -e "DNA alignment ${dna_alignment}"
    ${pipe_dir}/Tools/samtools view -@ ${num_thread} -h -b -q 30 -o ${name}.Barcoded.UniqAlign.DNA.bam -F 2304 ${dna_alignment}
    ${pipe_dir}/Tools/samtools stats -@ ${num_thread} ${name}.Barcoded.UniqAlign.DNA.bam > ${name}.Barcoded.UniqAlign.DNA.stat
    ${pipe_dir}/Tools/bedtools intersect -v -a ${name}.Barcoded.UniqAlign.DNA.bam -b ${pipe_dir}/RepeatMask/hg38_blacklist_rmsk.milliDivLessThan140.bed.gz > ${name}.Barcoded.UniqAlign.Masked.DNA.bam
    ${pipe_dir}/Tools/BARP p2s -b "COMPLEX" ${name}.Barcoded.UniqAlign.Masked.DNA.bam -o ${name}.DNA
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/RemoveDuplication.py -i ${name}.DNA.cluster -r "CHR-START" -o ${name}.DNA
    ${pipe_dir}/Tools/python3  ${pipe_dir}/PythonScript/GetFrag2FragDistance.py ${name}.DNA.RDP.cluster 2
    awk '{OFS="\t"; print $1,$2,$3,$4}' ${name}.DNA.RDP.cluster.FragLen.stat > ${name}.DNA.rgn
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/ChIADropSubGEM.py ${name}.DNA.RDP.cluster 2
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GetFrag2FragDistance.py ${name}.DNA.RDP.cluster.SUBGEM 2 
    awk '{OFS="\t"; print $1,$2,$3,$4}' ${name}.DNA.RDP.cluster.SUBGEM.FragLen.stat > ${name}.DNA.SUBGEM.rgn
    # awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$9,$7}' L2L.stat > ${name}.DNA.loop
    if [[ ! -d StatFiles ]]
    then
        mkdir StatFiles
    fi
    mv *stat StatFiles
    rm ${name}.DNA_P2S.bam
    # sort -k2 -k6 ${name}.DNA.RDP.cluster.SUBGEM.juice > ${name}.DNA.RDP.cluster.SUBGEM.sorted.juice
    # java -jar ${pipe_dir}/Tools/juicer_tools_1.19.01.jar pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 ${name}.RDP.DNA.cluster.SUBGEM.sorted.juice ${name}.hic ${pipe_dir}/ChromSize/${genomesize}
elif [[ (${expe_type} == "chiadrop") ]]
then 
    echo -e "${pipe_dir}/Tools/samtools view -@ ${num_thread} -h -b -F 2304 -q 30 ${dna_alignment} -o ${name}.F2304.q30.bam"
    ${pipe_dir}/Tools/samtools view -@ ${num_thread} -h -b -F 2304 -q 30 ${dna_alignment} -o ${name}.F2304.q30.bam
    echo -e "Statistic for Uniquely mappable reads."
    ${pipe_dir}/Tools/samtools stats -@ ${num_thread} ${name}.F2304.q30.bam > ${name}.UniqAlign.stat
    echo "Grouping reads to complex."
    ${pipe_dir}/Tools/BARP p2s -b "COMPLEX" --read1 ${name}.F2304.q30.bam -o ${name}
    echo "Merging fragments with complex."
    ${pipe_dir}/Tools/BARP merge -d 3000 -e 100 ${name}.cluster
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/ChIADropSubGEM.py Merged_Cluster.txt 2
    mv Merged_Cluster.txt.SUBGEM ${name}.ChIADrop.SubGEM
    grep -v Het ${name}.ChIADrop.SubGEM | grep -v chrU | grep -v chrM > ${name}.ChIADrop.SubGEM.Clean
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GetFrag2FragDistanceChIADrop.py ${name}.ChIADrop.SubGEM.Clean 2
    awk '{OFS="\t"; print $1,$2,$3,$4}' ${name}.ChIADrop.SubGEM.Clean.FragLen.stat > ${name}.rgn
    awk '{OFS="\t"; print $1,$2,$3,$4,$5,$6,$9,$7}' ${name}.ChIADrop.SubGEM.Clean.L2L.stat > ${name}.loop
    sort -k2 -k6 ${name}.ChIADrop.SubGEM.Clean.juice > ${name}.ChIADrop.SubGEM.Clean.sorted.juice
    if [[ ${genomesize} =~ ${pipe_dir} ]] 
    then
        tgenome=${genomesize}
    else
        tgenome=${pipe_dir}/ChromSize/${genomesize}
    fi
    java -jar ${pipe_dir}/Tools/juicer_tools_1.19.01.jar pre -r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,1000 ${name}.ChIADrop.SubGEM.Clean.sorted.juice ${name}.hic ${tgenome}
    if [[ ! -d StatFiles ]]
    then
        mkdir StatFiles
    fi
    mv *stat StatFiles
    rm ${name}_P2S.bam
    if [[ $? != 0 ]]
    then
        FinishSuccess=0
    fi
elif [[ ${expe_type} == "arcatac" ]]
then
    echo -e "Filtering qualified reads."
    echo -e "samtools view -F 2304 ${dna_alignment} -h -q 30 -o ${name}_DNA.F2304.q30.bam -@ ${num_thread}"
    ${pipe_dir}/Tools/samtools view -F 2304 ${dna_alignment} -h -q 30 -o ${name}_DNA.F2304.q30.bam -@ ${num_thread}
    echo -e "Clustering fragments"
    echo -e "${pipe_dir}/Tools/BARP p2s -b CELL --strict --atac ${name}_DNA.F2304.q30.bam"
    ${pipe_dir}/Tools/BARP p2s -b CELL --strict --atac ${name}_DNA.F2304.q30.bam -o ${name} 2>/dev/null
    echo -e "Removing duplicates"
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/RemoveDuplication.py -i ${name}.cluster -r "CHR-START-END" -o ${name}
    echo -e "Converting inline bed file to bed file"
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/ConvertInlineBedToBed.py -i ${name}.RDP.cluster -a 1 -o ${name}.Qualified
    echo -e "Sorting bed."
    sort -k1,1 -k2,2n ${name}.Qualified.bed > ${name}.Qualified.Sorted.bed
    echo -e "Generating coverage"
    if [[ ${genomesize} =~ ${pipe_dir} ]] 
    then
        tgenome=${genomesize}
    else
        tgenome=${pipe_dir}/ChromSize/${genomesize}
    fi
    ${pipe_dir}/Tools/bedtools genomecov -bg -i ${name}.Qualified.Sorted.bed -g ${tgenome} > ${name}.QualifiedFragments.bedgraph
    echo -e "Call peak"
    ${pipe_dir}/Tools/macs2 bdgpeakcall -i ${name}.QualifiedFragments.bedgraph -o ${name}.narrowPeaks
    echo -e "Intersect single cell fragment with peak"
    ${pipe_dir}/Tools/bedtools intersect -a ${name}.Qualified.Sorted.bed -b ${name}.narrowPeaks -wa -wb | ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GeneratePeakCellCount.py -o ${name}.PeakCellCount
    echo -e "Generate matrix"
    if [[ ! -d raw_matrix ]]
    then
        mkdir raw_matrix
    fi
    cd raw_matrix
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GenerateMatrixForATAC.py ../${name}.PeakCellCount.bed
    Rscript ${pipe_dir}/Rscript/RemoveEmptyDropletForATAC.R 
    cd ..
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/SynchronizeATACAndRNA.py ../RNAResult/01.BarcodeIden/${name}.CELL.idb ../01.BarcodeIden/${name}.CELL.idb ${name}.PeakCellCount.bed 4
    grep -v ATAC ${name}.PeakCellCount.bedSyn > ${name}.Overlap.PeakCellCount.bed
    if [[ ! -d filtered_matrix ]]
    then
        mkdir filtered_matrix
    fi
    cd filtered_matrix
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/RemovePredominantlyDoublets.py ../${name}.RDP.cluster ../raw_matrix/filtered_cells.tsv
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/FilterCellForATAC.py ../${name}.PeakCellCount.bed Cell
    cd ..

    if [[ ! -d translated_raw_matrix ]]
    then
        mkdir translated_raw_matrix
    fi
    cd translated_raw_matrix
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GenerateMatrixForATAC.py ../${name}.Overlap.PeakCellCount.bed
    Rscript ${pipe_dir}/Rscript/RemoveEmptyDropletForATAC.R 
    cd ..

    if [[ ! -d translated_filtered_matrix ]]
    then
        mkdir translated_filtered_matrix
    fi
    cd translated_filtered_matrix
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/RemovePredominantlyDoublets.py ../${name}.RDP.cluster ../translated_raw_matrix/filtered_cells.tsv
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/FilterCellForATAC.py ../${name}.Overlap.PeakCellCount.bed Cell
    if [[ $? != 0 ]]
    then
        FinishSuccess=0
    fi
    cd ..
elif [[ ${expe_type} == "scatac" ]]
then
    echo -e "Filtering qualified reads."
    echo -e "samtools view -F 2304 ${dna_alignment} -h -q 30 -o ${name}_DNA.F2304.q30.bam -@ ${num_thread}"
    ${pipe_dir}/Tools/samtools view -F 2304 ${dna_alignment} -h -q 30 -o ${name}_DNA.F2304.q30.bam -@ ${num_thread}
    echo -e "Clustering fragments"
    echo -e "${pipe_dir}/Tools/BARP p2s -b CELL --strict --atac ${name}_DNA.F2304.q30.bam"
    ${pipe_dir}/Tools/BARP p2s -b CELL --strict --atac ${name}_DNA.F2304.q30.bam -o ${name} 2>/dev/null
    echo -e "Removing duplicates"
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/RemoveDuplication.py -i ${name}.cluster -r "CHR-START-END" -o ${name}
    echo -e "Converting inline bed file to bed file"
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/ConvertInlineBedToBed.py -i ${name}.RDP.cluster -a 1 -o ${name}.Qualified
    echo -e "Sorting bed."
    sort -k1,1 -k2,2n ${name}.Qualified.bed > ${name}.Qualified.Sorted.bed
    echo -e "Generating coverage"
    if [[ ${genomesize} =~ ${pipe_dir} ]] 
    then
        tgenome=${genomesize}
    else
        tgenome=${pipe_dir}/ChromSize/${genomesize}
    fi
    ${pipe_dir}/Tools/bedtools genomecov -bg -i ${name}.Qualified.Sorted.bed -g ${tgenome} > ${name}.QualifiedFragments.bedgraph
    echo -e "Call peak"
    ${pipe_dir}/Tools/macs2 bdgpeakcall -i ${name}.QualifiedFragments.bedgraph -o ${name}.narrowPeaks
    echo -e "Intersect single cell fragment with peak"
    ${pipe_dir}/Tools/bedtools intersect -a ${name}.Qualified.Sorted.bed -b ${name}.narrowPeaks -wa -wb | ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GeneratePeakCellCount.py -o ${name}.PeakCellCount
    echo -e "Generate matrix"
    if [[ ! -d raw_matrix ]]
    then
        mkdir raw_matrix
    fi
    cd raw_matrix
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/GenerateMatrixForATAC.py ../${name}.PeakCellCount.bed
    Rscript ${pipe_dir}/Rscript/RemoveEmptyDropletForATAC.R 
    cd ..
    if [[ ! -d filtered_matrix ]]
    then
        mkdir filtered_matrix
    fi
    cd filtered_matrix
    echo -e "Removing barcode of predominatlyDoublets..."
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/RemovePredominantlyDoublets.py ../${name}.RDP.cluster ../raw_matrix/filtered_cells.tsv
    ${pipe_dir}/Tools/python3 ${pipe_dir}/PythonScript/FilterCellForATAC.py ../${name}.PeakCellCount.bed Cell
    cd ..
    if [[ $? != 0 ]]
    then
        FinishSuccess=0
    fi
fi

if [[ $FinishSuccess != 1 ]]
then
    echo "Group and data refinement failed."
else
    cd ..
    echo "$0 done" > GroupAndDataRefine.done
fi

