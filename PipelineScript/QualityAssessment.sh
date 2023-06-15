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
	echo -e "\t-@ Thread to use"
}
name="-"
while getopts t:p:n:f:@:1:2:g:d:r:k:h flag
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
        ? | h) print_help
            exit 1;;
    esac
done



if [[ ! -d 04.QualityAssess ]]
then
    mkdir 04.QualityAssess
    cd 04.QualityAssess
else 
    echo -e "\nDirectory 04.QualityAssess exist, will cover the files generated before in 5 seconds...";
    sleep 5
    cd 04.QualityAssess
    rm * 2>/dev/null
fi
process_date=$(date +%Y%b%d | tr [:lower:] [:upper:])

if [[ $name == "-" ]]
then
    name=${process_date}
fi
echo -e "\nQuality assessment..."
echo -e "Working directory: $(pwd)"
echo ""

if [[ ${expe_type} == "chiadrop" ]]
then
    Total_read_pairs=$(grep "Total reads" ../01.BarcodeIden/${name}.stat | awk -F ":" '{print $2}')
    Read_pairs_with_full_barcodes=$(grep -F 'Fully barcoded' ../01.BarcodeIden/${name}.stat | awk -F ":" '{sum+=$2} END {print sum}')
    Fully_barcode_rate=$(echo "scale=6; ${Read_pairs_with_full_barcodes}/${Total_read_pairs}" | bc)
    UniqMapReads=$(grep "^SN" ../03.GroupAndRefine/StatFiles/${name}.UniqAlign.stat | head -n 6 | tail -n 1 | awk '{print $4}')
    ComplexCountAtFastq=$(grep COMPLEX ../01.BarcodeIden/BarcodeIdentification_statistic.txt | awk -F ":" '{print $2}')
    TotalFragment=$(awk '{sum+=$2} END {print sum}' ../03.GroupAndRefine/${name}.ChIADrop.SubGEM.Clean)
    DuplicatedFragment="-"
    DuplicationRate="-"
    RefinedComplex=$(wc -l ../03.GroupAndRefine/${name}.ChIADrop.SubGEM.Clean | awk '{print $1}')
    echo -e "Total_read_pairs\t${Total_read_pairs}" > ${name}_final_stat.tsv
    echo -e "Read_pairs_with_full_barcodes\t${Read_pairs_with_full_barcodes}" >> ${name}_final_stat.tsv
    echo -e "Fully_barcode_rate\t${Fully_barcode_rate}" >> ${name}_final_stat.tsv
    echo -e "Complex_count_at_fastq\t${ComplexCountAtFastq}" >> ${name}_final_stat.tsv
    echo -e "Total_fragments\t${TotalFragment}" >> ${name}_final_stat.tsv
    echo -e "Duplicated_fragments\t${DuplicatedFragment}" >> ${name}_final_stat.tsv
    echo -e "Duplication_rate\t${DuplicationRate}" >> ${name}_final_stat.tsv
    echo -e "Uniquely_mapped_reads\t${UniqMapReads}" >> ${name}_final_stat.tsv
    echo -e "Refined_complex\t${RefinedComplex}" >> ${name}_final_stat.tsv
    cat ../03.GroupAndRefine/StatFiles/${name}.ChIADrop.SubGEM.Clean.FragCount.stat >> ${name}_final_stat.tsv
    ${pipe_dir}/Tools/Rscript ${pipe_dir}/Rscript/PlotComplexDist.R ../03.GroupAndRefine/StatFiles/${name}.ChIADrop.SubGEM.Clean.ComplexSpan.stat
    ${pipe_dir}/Tools/Rscript ${pipe_dir}/Rscript/PlotF2FDist.R ../03.GroupAndRefine/StatFiles/${name}.ChIADrop.SubGEM.Clean.F2F.stat
    ${pipe_dir}/Tools/Rscript ${pipe_dir}/Rscript/PlotFragLen.R ../03.GroupAndRefine/StatFiles/${name}.ChIADrop.SubGEM.Clean.FragLen.stat
    ${pipe_dir}/Tools/Rscript ${pipe_dir}/Rscript/PlotLoopSpans.R ../03.GroupAndRefine/StatFiles/${name}.ChIADrop.SubGEM.Clean.L2L.stat
elif [[ ${expe_type} == "sprite" ]]
then
    Total_read_pairs=$(grep "Total reads" ../01.BarcodeIden/${name}.stat | awk -F ":" '{print $2}')
    Read_pairs_with_full_barcodes=$(grep -F 'Fully barcoded' ../01.BarcodeIden/${name}.stat | awk -F ":" '{sum+=$2} END {print sum}')
    Fully_barcode_rate=$(echo "scale=6; ${Read_pairs_with_full_barcodes}/${Total_read_pairs}" | bc)
    UniqMapReads=$(grep "^SN" ../03.GroupAndRefine/StatFiles/${name}.Barcoded.UniqAlign.DNA.stat | head -n 7 | tail -n 1 | awk '{print $4}')
    ComplexCountAtFastq=$(grep COMPLEX ../01.BarcodeIden/BarcodeIdentification_statistic.txt | awk -F ":" '{print $2}')
    TotalFragment=$(awk '{print $2}' ../03.GroupAndRefine/StatFiles/${name}.DNA.Dup.stat | awk -F ":" '{print $2}' )
    DuplicatedFragment=$(awk '{print $1}' ../03.GroupAndRefine/StatFiles/${name}.DNA.Dup.stat | awk -F ":" '{print $2}' )
    DuplicationRate=$(echo "scale=6; ${DuplicatedFragment}/${TotalFragment}" | bc)
    RefinedComplex=$(echo "`wc -l ../03.GroupAndRefine/${name}.DNA.RDP.cluster | awk '{print $1}'` -2" | bc)
    echo -e "Total_read_pairs\t${Total_read_pairs}" > ${name}_final_stat.tsv
    echo -e "Read_pairs_with_full_barcodes\t${Read_pairs_with_full_barcodes}" >> ${name}_final_stat.tsv
    echo -e "Fully_barcode_rate\t${Fully_barcode_rate}" >> ${name}_final_stat.tsv
    echo -e "Complex_count_at_fastq\t${ComplexCountAtFastq}" >> ${name}_final_stat.tsv
    echo -e "Uniquely_mapped_reads\t${UniqMapReads}" >> ${name}_final_stat.tsv
    echo -e "Total_fragments\t${TotalFragment}" >> ${name}_final_stat.tsv
    echo -e "Duplicated_fragments\t${DuplicatedFragment}" >> ${name}_final_stat.tsv
    echo -e "Duplication_rate\t${DuplicationRate}" >> ${name}_final_stat.tsv
    echo -e "Total_qualified_fragment\t${RefinedComplex}" >> ${name}_final_stat.tsv
    cat ../03.GroupAndRefine/StatFiles/${name}.DNA.RDP.cluster.FragCount.stat >> ${name}_final_stat.tsv
elif [[ ${expe_type} == "scrna" ]]
then
    Total_read_pairs=$(grep "Total reads" ../01.BarcodeIden/${name}.stat | awk -F ":" '{print $2}')
    Read_pairs_with_full_barcodes=$(grep -F 'Fully barcoded' ../01.BarcodeIden/${name}.stat | awk -F ":" '{sum+=$2} END {print sum}')
    Fully_barcode_rate=$(echo "scale=6; ${Read_pairs_with_full_barcodes}/${Total_read_pairs}" | bc)
    UniqMapReads=$(grep "Unique Reads in Cells Mapped to Gene" ../02.ReadAlign/${name}Solo.out/Gene/Summary.csv | awk -F "," '{print $2}')
    CellCountAtFastq=$(grep CELL ../01.BarcodeIden/BarcodeIdentification_statistic.txt | awk -F ":" '{print $2}')
    CellCountEstimate=$(grep "Estimated Number of Cells" ../02.ReadAlign/${name}Solo.out/Gene/Summary.csv | awk -F "," '{print $2}')
    Total_Gene_Detected=$(grep "Total Gene Detected" ../02.ReadAlign/${name}Solo.out/Gene/Summary.csv | awk -F "," '{print $2}')
    Mean_Gene_per_Cell=$(grep "Mean Gene per Cell" ../02.ReadAlign/${name}Solo.out/Gene/Summary.csv | awk -F "," '{print $2}')
    Median_UMI_per_Cell=$(grep "Median UMI per Cell" ../02.ReadAlign/${name}Solo.out/Gene/Summary.csv | awk -F "," '{print $2}')
    echo -e "Total_read_pairs\t${Total_read_pairs}" > ${name}_final_stat.tsv
    echo -e "Read_pairs_with_full_barcodes\t${Read_pairs_with_full_barcodes}" >> ${name}_final_stat.tsv
    echo -e "Fully_barcode_rate\t${Fully_barcode_rate}" >> ${name}_final_stat.tsv
    echo -e "Cell_count_at_fastq\t${CellCountAtFastq}" >> ${name}_final_stat.tsv
    echo -e "Cell_count_estimated\t${CellCountEstimate}" >> ${name}_final_stat.tsv
    echo -e "Total_gene_detected\t${Total_Gene_Detected}" >> ${name}_final_stat.tsv
    echo -e "Mean_gene_per_cell\t${Mean_Gene_per_Cell}" >> ${name}_final_stat.tsv
    echo -e "Median_UMI_per_cell\t${Median_UMI_per_Cell}" >> ${name}_final_stat.tsv
elif [[ ${expe_type} == "scatac" ]]
then
    Total_read_pairs=$(grep "Total reads" ../01.BarcodeIden/${name}.stat | awk -F ":" '{print $2}')
    Read_pairs_with_full_barcodes=$(grep -F 'Fully barcoded' ../01.BarcodeIden/${name}.stat | awk -F ":" '{sum+=$2} END {print sum}')
    Fully_barcode_rate=$(echo "scale=6; ${Read_pairs_with_full_barcodes}/${Total_read_pairs}" | bc)
    CellCountAtFastq=$(grep CELL ../01.BarcodeIden/BarcodeIdentification_statistic.txt | awk -F ":" '{print $2}')
    TotalFragment=$(awk '{print $2}' ../03.GroupAndRefine/${name}.Dup.stat | awk -F ":" '{print $2}' )
    DuplicatedFragment=$(awk '{print $1}' ../03.GroupAndRefine/${name}.Dup.stat | awk -F ":" '{print $2}' )
    DuplicationRate=$(echo "scale=6; ${DuplicatedFragment}/${TotalFragment}" | bc)
    Fragments_overlaps_peak=$(awk '{sum+=$5} END {print sum}' ../03.GroupAndRefine/${name}.PeakCellCount.bed)
    PeakCount=$(wc -l ../03.GroupAndRefine/filtered_matrix/Peak | awk '{print $1}')
    CellCount=$(wc -l ../03.GroupAndRefine/filtered_matrix/Cell | awk '{print $1}')

    echo -e "Total_read_pairs\t${Total_read_pairs}" > ${name}_final_stat.tsv
    echo -e "Read_pairs_with_full_barcodes\t${Read_pairs_with_full_barcodes}" >> ${name}_final_stat.tsv
    echo -e "Fully_barcode_rate\t${Fully_barcode_rate}" >> ${name}_final_stat.tsv
    echo -e "Cell_count_at_fastq\t${CellCountAtFastq}" >> ${name}_final_stat.tsv
    echo -e "Total_fragments\t${TotalFragment}" >> ${name}_final_stat.tsv
    echo -e "Duplicated_fragments\t${DuplicatedFragment}" >> ${name}_final_stat.tsv
    echo -e "Duplication_rate\t${DuplicationRate}" >> ${name}_final_stat.tsv
    echo -e "Fragments_overlap_peak\t${Fragments_overlaps_peak}" >> ${name}_final_stat.tsv
    echo -e "Peak_count\t${PeakCount}" >> ${name}_final_stat.tsv
    echo -e "Cell_count\t${CellCount}" >> ${name}_final_stat.tsv
elif [[ ${expe_type} == "scsprite" ]]
then
    Total_read_pairs=$(grep "Total reads" ../01.BarcodeIden/${name}.stat | awk -F ":" '{print $2}')
    Read_pairs_with_full_barcodes=$(grep -F 'Fully barcoded' ../01.BarcodeIden/${name}.stat | awk -F ":" '{sum+=$2} END {print sum}')
    Fully_barcode_rate=$(echo "scale=6; ${Read_pairs_with_full_barcodes}/${Total_read_pairs}" | bc)
    UniqMapReads=$(grep "^SN" ../03.GroupAndRefine/StatFiles/${name}.Barcoded.UniqAlign.stat | head -n 7 | tail -n 1 | awk '{print $4}')
    ComplexCountAtFastq=$(grep COMPLEX ../01.BarcodeIden/BarcodeIdentification_statistic.txt | awk -F ":" '{print $2}')
    CellCountAtFastq=$(grep CELL ../01.BarcodeIden/BarcodeIdentification_statistic.txt | awk -F ":" '{print $2}')
    TotalFragment=$(awk '{print $2}' ../03.GroupAndRefine/StatFiles/${name}.DNA.Dup.stat | awk -F ":" '{print $2}' )
    DuplicatedFragment=$(awk '{print $1}' ../03.GroupAndRefine/StatFiles/${name}.DNA.Dup.stat | awk -F ":" '{print $2}' )
    DuplicationRate=$(echo "scale=6; ${DuplicatedFragment}/${TotalFragment}" | bc)
    RefinedComplex=$(echo "`wc -l ../03.GroupAndRefine/${name}.DNA.RDP.cluster | awk '{print $1}'` -2" | bc)
    echo -e "Total_read_pairs\t${Total_read_pairs}" > ${name}_final_stat.tsv
    echo -e "Read_pairs_with_full_barcodes\t${Read_pairs_with_full_barcodes}" >> ${name}_final_stat.tsv
    echo -e "Fully_barcode_rate\t${Fully_barcode_rate}" >> ${name}_final_stat.tsv
    echo -e "Cell_count_at_fastq\t${CellCountAtFastq}" >> ${name}_final_stat.tsv
    echo -e "Complex_count_at_fastq\t${ComplexCountAtFastq}" >> ${name}_final_stat.tsv
    echo -e "Uniquely_mapped_reads\t${UniqMapReads}" >> ${name}_final_stat.tsv
    echo -e "Total_fragments\t${TotalFragment}" >> ${name}_final_stat.tsv
    echo -e "Duplicated_fragments\t${DuplicatedFragment}" >> ${name}_final_stat.tsv
    echo -e "Duplication_rate\t${DuplicationRate}" >> ${name}_final_stat.tsv
    echo -e "Total_qualified_fragment\t${RefinedComplex}" >> ${name}_final_stat.tsv
    cat ../03.GroupAndRefine/StatFiles/${name}.DNA.RDP.cluster.FragCount.stat >> ${name}_final_stat.tsv
elif [[ ${expe_type} == "rdsprite" ]]
then
    Total_read_pairs=$(grep "Total reads" ../01.BarcodeIden/${name}.stat | awk -F ":" '{print $2}')
    Read_pairs_with_full_barcodes=$(grep -F 'Fully barcoded' ../01.BarcodeIden/${name}.stat | awk -F ":" '{sum+=$2} END {print sum}')
    DNA_read_with_full_barcodes=$(grep -F 'Fully barcoded DNA reads' ../01.BarcodeIden/${name}.stat | awk -F ":" '{print $2}')
    RNA_read_with_full_barcodes=$(grep -F 'Fully barcoded RNA reads' ../01.BarcodeIden/${name}.stat | awk -F ":" '{print $2}')
    UniqMapDNAReads=$(grep "^SN" ../03.GroupAndRefine/StatFiles/${name}.Barcoded.UniqAlign.DNA.stat | head -n 7 | tail -n 1 | awk '{print $4}')
    UniqMapRNAReads=$(grep "^SN" ../03.GroupAndRefine/StatFiles/${name}.Barcoded.UniqAlign.RNA.stat | head -n 7 | tail -n 1 | awk '{print $4}')
    Fully_barcode_rate=$(echo "scale=6; ${Read_pairs_with_full_barcodes}/${Total_read_pairs}" | bc)
    DNA_cluster=$(wc -l ../03.GroupAndRefine/${name}.DNA.RDP.cluster | awk '{print $1}')
    RNA_cluster=$(wc -l ../03.GroupAndRefine/${name}.RNA.RDP.cluster | awk '{print $1}')
    TotalDNAFragment=$(awk '{print $2}' ../03.GroupAndRefine/StatFiles/${name}.DNA.Dup.stat | awk -F ":" '{print $2}' )
    TotalRNAFragment=$(awk '{print $2}' ../03.GroupAndRefine/StatFiles/${name}.RNA.Dup.stat | awk -F ":" '{print $2}' )
    DuplicatedDNAFragment=$(awk '{print $1}' ../03.GroupAndRefine/StatFiles/${name}.DNA.Dup.stat | awk -F ":" '{print $2}' )
    DuplicatedRNAFragment=$(awk '{print $1}' ../03.GroupAndRefine/StatFiles/${name}.RNA.Dup.stat | awk -F ":" '{print $2}' )
    DuplicationDNARate=$(echo "scale=6; ${DuplicatedDNAFragment}/${TotalDNAFragment}" | bc)
    DuplicationRNARate=$(echo "scale=6; ${DuplicatedRNAFragment}/${TotalRNAFragment}" | bc)
    echo -e "Total_read_pairs\t${Total_read_pairs}" > ${name}_final_stat.tsv
    echo -e "Read_pairs_with_full_barcode\t${Read_pairs_with_full_barcodes}" >> ${name}_final_stat.tsv
    echo -e "DNA_read_pairs_with_full_barcode\t${DNA_read_with_full_barcodes}" >> ${name}_final_stat.tsv
    echo -e "RNA_read_pairs_with_full_barcode\t${RNA_read_with_full_barcodes}" >> ${name}_final_stat.tsv
    echo -e "Fully_barcode_rate\t${Fully_barcode_rate}" >> ${name}_final_stat.tsv
    echo -e "Uniquely_mapped_DNA_reads\t${UniqMapDNAReads}" >> ${name}_final_stat.tsv
    echo -e "Uniquely_mapped_RNA_reads\t${UniqMapRNAReads}" >> ${name}_final_stat.tsv
    echo -e "DNA_fragments\t${TotalDNAFragment}" >> ${name}_final_stat.tsv
    echo -e "RNA_fragments\t${TotalRNAFragment}" >> ${name}_final_stat.tsv
    echo -e "DNA_duplicated_fragments\t${DuplicatedDNAFragment}" >> ${name}_final_stat.tsv
    echo -e "RNA_duplicated_fragments\t${DuplicatedRNAFragment}" >> ${name}_final_stat.tsv
    echo -e "DNA_duplicate_rate\t${DuplicationDNARate}" >> ${name}_final_stat.tsv
    echo -e "RNA_duplicate_rate\t${DuplicationRNARate}" >> ${name}_final_stat.tsv
    echo -e "DNA_complex\t${DNA_cluster}" >> ${name}_final_stat.tsv
    cat ../03.GroupAndRefine/StatFiles/${name}.DNA.RDP.cluster.FragCount.stat >> ${name}_final_stat.tsv
    echo -e "RNA_complex\t${RNA_cluster}" >> ${name}_final_stat.tsv
    cat ../03.GroupAndRefine/StatFiles/${name}.RNA.RDP.cluster.FragCount.stat >> ${name}_final_stat.tsv
fi

cd ..
echo "Finished quality assessment."
echo "$0 done" > QualityAssessment.done

