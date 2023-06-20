# ScSmOP

Single cell Single Molecule Multiple Omics Pipeline.

-----------

## Summary

ScSmOP is a universal pipeline capable of performing data processing for a wide range of single-cell single-molecule omics like scRNA-seq, scATAC-seq, ChIA-Drop, SPRITE and its derivates, Drop-seq and some other techniques which based on barcode and UMI. It have pre-prepared several pipelines for current popular which can be easily launched with some simple editions.

## OS env

ScSmOP is running in Ubuntu 16.4 and above.

## Installation

First, download ScSmOP to your computer:

```
    # Get ScSmOP source using git
    :~$ git clone https://github.com/kasenjing/ScSmOP.git
    :~$ cd ScSmOP

    # Alternatively, get latest ScSmOP source from releases
    :~$ wget https://github.com/kasenjing/ScSmOP/archive/refs/tags/v0.1.3.tar.gz
    :~$ tar -zxvf v0.1.3.tar.gz
    :~$ cd ScSmOP-0.1.3
```

### Conda installed user

If you are familiar with linux operate system and conda, after activate conda, you can install ScSmOP using the one-click `local_install_ScSmOP_dependencies_require_conda_exist.sh`

```
    (base) :~$ bash local_install_ScSmOP_dependencies_require_conda_exist.sh
```

It will automatically create a conda environment `ScSmOP` and install all the softwares ScSmOP required.

### Do not have conda in computer

For users who are not familiar with linux operate system or conda, we prepared one-click `local_install_ScSmOP_dependencies_install_conda_for_you.sh` to install ScSmOP:

```
    :~$ bash local_install_ScSmOP_dependencies_install_conda_for_you.sh
```

It will create conda environment and install all the softwares ScSmOP required.

***Warings during installation do not affect ScSmOP usage, ignore them.***

#### Activate conda env

```
    :~$ ./Tools/miniconda3/bin/conda init bash
    :~$ exit
    # Reopen a terminal
    (base) :~$ conda activate ScSmOP
    (ScSmOP) :~$ 
```

#### Suppose ScSmOP has been installed @

```
    ~/ScSmOP/
```

## Make index for BWA

*If you have BWA index already, please ignore this step.*
Reference [BWA manual](https://bio-bwa.sourceforge.net/bwa.shtml)

```
    (ScSmOP) :~$ mkdir RefGenome
    (ScSmOP) :~$ cd RefGenome
    (ScSmOP) :~/RefGenome$ mkdir bwa_hg38_index
    (ScSmOP) :~/RefGenome$ cd bwa_hg38_index
    (ScSmOP) :~/RefGenome/bwa_hg38_index$ wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    (ScSmOP) :~/RefGenome/bwa_hg38_index$ gunzip hg38.fa.gz
    (ScSmOP) :~/RefGenome/bwa_hg38_index$ bwa index hg38.fa
```


#### Suppose BWA index has been generated @

```
    ~/RefGenome/bwa_hg38_index/
```

## Make index for STAR

*If you have STAR index already, please ignore this step.*
Reference [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

Download annotation file (gencode.v29.primary_assembly.annotation_UCSC_names.gtf) from UCSC Genome Browser <https://genome.ucsc.edu/cgi-bin/hgTables>
Or make your single-cell RNA-seq result more like 10× Genomics, download the genome fasta file and gene annotation at <https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest>. The FASTA and GTF files:

```
    refdata-gex-GRCh38-2020-A/genes/genes.gtf
    refdata-gex-GRCh38-2020-A/fasta/genome.fa
```

Then build STAR reference genome:

```
    (ScSmOP) :~/RefGenome$ mkdir refdata-gex-GRCh38-2020-A-STAR
    (ScSmOP) :~/RefGenome$ cd refdata-gex-GRCh38-2020-A-STAR
    (ScSmOP) :~/RefGenome/refdata-gex-GRCh38-2020-A-STAR$ ~/ScSmOP/Tools/STAR --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir . \
    --genomeFastaFiles genome.fa \
    --sjdbGTFfile genes.gtf 
```

*Note: to build hg38 STAR reference genome require at least 32 GB memory.*

#### Suppose STAR index has been generated @

```
    ~/RefGenome/refdata-gex-GRCh38-2020-A-STAR/
```

## Run ScSmOP on different experiment data

We've prepared several pipelines for easy data process, if your experiment meet the listed experiment, you can process your data as in [Easy Process](#easy-process).

<html>
 <head>
  <meta charset="utf-8" />
  <title>SheetJS Table Export</title>
 </head>
 <body>
  <table>
   <tr>
    <td t="s" id="sjs-A1">Omics</td>
    <td t="s" id="sjs-B1">Experiment</td>
    <td t="s" id="sjs-C1">Barcode(s) length</td>
    <td t="s" id="sjs-D1">UMI length</td>
    <td t="s" id="sjs-E1" white-space: nowrap>Read structure</td>
    <td t="s" id="sjs-F1" white-space: nowrap>Parameter (-t)</td>
    <td t="s" id="sjs-G1">Source</td>
   </tr>
   <tr>
    <td rowspan="3" t="s" id="sjs-A2" white-space: nowrap>scRNA-seq</td>
    <td t="s" id="sjs-B2">10× Genomics (v2) Single Cell Gene Expression</td>
    <td t="n" id="sjs-C2">16</td>
    <td t="n" id="sjs-D2">10</td>
    <td t="s" id="sjs-E2">R1: Barcode-UMI&#x000d;
     <br />R2: Transcript</td>
    <td t="s" id="sjs-F2">scrna_10x_v2</td>
    <td t="s" id="sjs-G2">https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt</td>
   </tr>
   <tr>
    <td t="s" id="sjs-B3">10× Genomics (v3) Single Cell Gene Expression</td>
    <td t="n" id="sjs-C3">16</td>
    <td t="n" id="sjs-D3">12</td>
    <td t="s" id="sjs-E3">R1: Barcode-UMI&#x000d;
     <br />R2: Transcript</td>
    <td t="s" id="sjs-F3">scrna_10x_v3</td>
    <td t="s" id="sjs-G3">https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz</td>
   </tr>
   <tr>
    <td t="s" id="sjs-B4">Drop-seq</td>
    <td t="n" id="sjs-C4">12</td>
    <td t="n" id="sjs-D4">8</td>
    <td t="s" id="sjs-E4" white-space: nowrap>R1: Barcode-UMI&#x000d;
     <br />R2: Transcript</td>
    <td t="s" id="sjs-F4">dropseq</td>
    <td t="z" id="sjs-G4"></td>
   </tr>
   <tr>
    <td t="s" id="sjs-A5">ChIA-Drop</td>
    <td t="s" id="sjs-B5">ChIA-Drop</td>
    <td t="n" id="sjs-C5">16</td>
    <td t="s" id="sjs-D5">N/A</td>
    <td t="s" id="sjs-E5" white-space: nowrap>R1: Barcode-DNA&#x000d;
     <br />R2: DNA</td>
    <td t="s" id="sjs-F5">chiadrop</td>
    <td t="s" id="sjs-G5">longranger-2.2.2/longranger-cs/2.2.2/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt-</td>
   </tr>
   <tr>
    <td t="s" id="sjs-A6" white-space: nowrap>scATAC-seq</td>
    <td t="s" id="sjs-B6">10× Genomics (v2) Single Cell ATAC</td>
    <td t="n" id="sjs-C6">16</td>
    <td t="s" id="sjs-D6">N/A</td>
    <td t="s" id="sjs-E6"white-space: nowrap>R1: DNA&#x000d;
     <br />R2: Barcode&#x000d;
     <br />R3: DNA</td>
    <td t="s" id="sjs-F6">scatac_10x_v1</td>
    <td t="s" id="sjs-G6">cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz</td>
   </tr>
   <tr>
    <td rowspan="3" t="s" id="sjs-A7" white-space: nowrap>SPRITE and derivates</td>
    <td t="s" id="sjs-B7">SPRITE</td>
    <td t="s" id="sjs-C7">4 kinds with 5 rounds</td>
    <td t="s" id="sjs-D7">N/A</td>
    <td t="s" id="sjs-E7" white-space: nowrap>R1: DPM-DNA&#x000d;
     <br />R2: Y-ODD-EVEN-ODD</td>
    <td t="s" id="sjs-F7">sprite</td>
    <td t="s" id="sjs-G7">GSE114242: GSE114242_human_config.txt.gz</td>
   </tr>
   <tr>
    <td t="s" id="sjs-B8">RD-SPRITE</td>
    <td t="s" id="sjs-C8">5 kinds with 5 rounds (label DNA with DPM, label RNA with RPM)</td>
    <td t="s" id="sjs-D8">N/A</td>
    <td t="s" id="sjs-E8" white-space: nowrap>R1: DNA/RNA&#x000d;
     <br />R2: Y-ODD-EVEN-ODD-DPM/RPM <br />
     (DPM: DNA, RPM: RNA)</td>
    <td t="s" id="sjs-F8">rdsprite</td>
    <td t="s" id="sjs-G8">https://github.com/GuttmanLab/sprite2.0-pipeline/blob/master/config.txt</td>
   </tr>
   <tr>
    <td t="s" id="sjs-B9">Single-cell SPRITE</td>
    <td t="s" id="sjs-C9">4 kinds with 6 rounds (label cells with first last 3 barcodes in Read 2, label complex with all six barcodes)</td>
    <td t="s" id="sjs-D9">N/A</td>
    <td t="s" id="sjs-E9" white-space: nowrap>R1: DNA&#x000d;
     <br />R2: Y-EVEN-ODD-EVEN-ODD-DPM</td>
    <td t="s" id="sjs-F9">scsprite</td>
    <td t="s" id="sjs-G9">https://github.com/caltech-bioinformatics-resource-center/Guttman_Ismagilov_Labs/blob/master/scSPRITE/misc/config_dpm6_y-stag_scSPRITE.txt</td>
   </tr>
   <tr>
    <td t="s" id="sjs-A10">scARC-seq</td>
    <td t="s" id="sjs-B10">10× Genomics Single Cell Multiome ATAC + Gene Expression</td>
    <td t="s" id="sjs-C10">Gene Expression library: 16, ATAC library: 16</td>
    <td t="s" id="sjs-D10">Gene Expression library: 12; ATAC library: N/A</td>
    <td t="s" id="sjs-E10">R1: DNA (R1 in ATAC)&#x000d;
     <br />R2: Barcode (R2 in ATAC)&#x000d;
     <br />R3: DNA (R3 in ATAC)&#x000d;
     <br />R4: Barcode-UMI (R1 in GEX)&#x000d;
     <br />R5: Transcript (R2 in GEX)</td>
    <td t="s" id="sjs-F10">scarc_10x_v1</td>
    <td t="s" id="sjs-G10">737K-arc-v1-scrna.txt:cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz; 737K-arc-v1-scatac.txt:cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz</td>
   </tr>
  </table>
 </body>
</html>

### Easy process

Copy the corresponding process pipeline to the directory you want to store the result, for example: ChIA-Drop.
```
    (ScSmOP) ~:$ ~/ScSmOP/scsmop.sh  -t chiadrop -n CHDP -1 CHDP_R1.fastq.gz -2 CHDP_R2.fastq.gz -b ~/RefGenome/bwa_dm3_ref/dm3.fa -s dm3.size.txt
```
Different experiment require different parameters:
```
    (ScSmOP) ~:$ ~/ScSmOP/scsmop.sh 

    Usage: scsmop.sh [options] -t [STR] -n [STR] -1 [FILE] -2 [FILE] -r [DIR]
        -t Experiment type: sprite, scsprite, rdsprite, chiadrop, scrna_10x_v3, scrna_10x_v2, dropseq, scatac_10x_v1, scarc_10x_v1.
        -n Library name.
        -1 Read 1 FASTQ.
        -2 Read 2 FASTQ.
        -3 Read 3 FASTQ.
        -4 Read 4 FASTQ.
        -5 Read 5 FASTQ.
        -r STAR reference genome directory.
        -b BWA reference genome directory.
        -c Custom configuration file.
        -w Custom whitelist.
        -@ Thread to use.
        -p Custom ScSmOP directory.
        -s Chromosome size file, need by chiadrop, scatac_10x_v1, scarc_10x_v1.

    Processing scARC-seq, specify FASTQs through -1 ATAC R1 -2 ATAC R2 -3 ATAC R3 -4 GEX R1 -5 GEX R2, others refer to the table at Github.
    Processing chiadrop, scatac_10x_v1, scarc_10x_v1 requires genome chromosome sizes file showing chromosome length, specify with -s [FILE].

    Universal pipeline for multi-omics data process: <https://github.com/ZhengmzLab/ScSmOP/wiki>.
```

## Customize for your own barcoding method

Please refer to the paper or the [wiki](https://github.com/ZhengmzLab/ScSmOP/wiki).
