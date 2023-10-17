# ScSmOP

Single cell Single Molecule Multiple Omics Pipeline.

Kai Jing, Yewen Xu, Yang Yang, Pengfei Yin, Duo Ning, Guangyu Huang, Yuqing Deng, Gengzhan Chen, Guoliang Li, Simon Zhongyuan Tian, Meizhen Zheng, **ScSmOP: a universal computational pipeline for single-cell single-molecule multiomics data analysis**, Briefings in Bioinformatics, Volume 24, Issue 6, November 2023, bbad343, https://doi.org/10.1093/bib/bbad343

-----------

ScSmOP is a general-purpose pipeline that works for processing the “one barcode one group” type of state-of-art multiomics data, such as the (i) “single-cell single-molecule multiple chromatin interaction” type of scSPRITE method, (ii) “single-molecule multi-way chromatin (DNA/RNA) interaction” type of ChIA-Drop, SPRITE, and RD-SPRITE methods, (iii) “single-cell multiomics methods” type of scRNA-seq, Drop-seq, scATAC-seq, sc(ATAC+RNA)-seq methods, (iv) and even “single-position single-cell expression” type of spatial transcriptomics. In summary, these “one barcode one group” barcoding methods are based on two types of strategies: the ligation-based barcoding strategy and synthesis-based barcoding strategy.  In this study, the developed ScSmOP can not only be applied to a variety of “barcode” types of single-cell and/or single-molecule multiplex omics techniques but also supplies an interface for users to “DIY” a pipeline for their customized barcodes. Detail description and cases of study please see [ScSmOP Wiki](https://github.com/ZhengmzLab/ScSmOP/wiki).

![image](https://github.com/tianzhongyuan/material/blob/main/ScSmOP_operation.png)

## 1. OS env

ScSmOP is running in Ubuntu 16.4 and above.

## 2. Installation

First, download ScSmOP to your computer:

```
# Get ScSmOP source using git
$ git clone https://github.com/ZhengmzLab/ScSmOP.git
$ cd ScSmOP

# Alternatively, get latest ScSmOP source from releases
$ wget https://github.com/ZhengmzLab/ScSmOP/archive/refs/tags/v1.0.tar.gz
$ tar -zxvf v1.0.tar.gz
$ cd ScSmOP-1.0/
```

### 2.1 Linux has Conda:
To install ScSmOP on a Linux system with a Conda virtual environment, you need to activate Conda first, and then install ScSmOP by executing the following command:
```
$ conda activate 
(base)$ bash local_install_ScSmOP_dependencies_require_conda_exist.sh
```
This will automatically create a Conda environment "ScSmOP" and install all the software required for ScSmOP.


### 2.2 Linux has no Conda:

To install ScSmOP on a Linux system without a Conda virtual environment, you only need to install ScSmOP by executing the following command:
```
$ bash local_install_ScSmOP_dependencies_install_conda_for_you.sh
```
This will automatically install Conda tools, create a Conda environment "ScSmOP" and install all the software required for ScSmOP.


***Warings during installation do not affect ScSmOP usage, ignore them.***

**Activate Conda env**

```
$ ./Tools/miniconda3/bin/conda init bash
$ exit

# Reopen a terminal
(base)$ conda activate ScSmOP
(ScSmOP)$ 
```

**Suppose ScSmOP has been installed**

```
~/ScSmOP/
```

## 3. Make index for BWA

*If you have BWA index already, please ignore this step.*
Reference [BWA manual](https://bio-bwa.sourceforge.net/bwa.shtml)

```
(ScSmOP)$ mkdir RefGenome
(ScSmOP)$ cd RefGenome
(ScSmOP)/RefGenome$ mkdir bwa_hg38_index
(ScSmOP)/RefGenome$ cd bwa_hg38_index
(ScSmOP)/RefGenome/bwa_hg38_index$ wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
(ScSmOP)/RefGenome/bwa_hg38_index$ gunzip hg38.fa.gz
(ScSmOP)/RefGenome/bwa_hg38_index$ bwa index hg38.fa
```


**Suppose BWA index has been generated**

```
~/RefGenome/bwa_hg38_index/
```

## 4. Make index for STAR

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
(ScSmOP)/RefGenome$ mkdir refdata-gex-GRCh38-2020-A-STAR
(ScSmOP)/RefGenome$ cd refdata-gex-GRCh38-2020-A-STAR
(ScSmOP)/RefGenome/refdata-gex-GRCh38-2020-A-STAR$ ~/ScSmOP/Tools/STAR --runThreadN 16 \
                                                                       --runMode genomeGenerate \
                                                                       --genomeDir . \
                                                                       --genomeFastaFiles genome.fa \
                                                                       --sjdbGTFfile genes.gtf 
```

*Note: to build hg38 STAR reference genome require at least 32 GB memory.*

**Suppose STAR index has been generated**

```
~/RefGenome/refdata-gex-GRCh38-2020-A-STAR/
```

## 5. Run ScSmOP on various types of data

We've prepared several pipelines for easy data process, if your experiment meet the listed experiment, you can process your data as in Easy Process.

### 5.1 [ScSmOP Supported Techology](#ScSmOPSuppTech)

|Omics|Technology|Barcode(s) length|UMI length|Read structure|Method ID (-t)|
|:----|:----|:----|:----|:----|:----|
|scRNA-seq|10× Genomics (v2) Single Cell Gene Expression|16|10|R1: Barcode-UMI; R2: Transcript|scrna_10x_v2|
|scRNA-seq|10× Genomics (v3) Single Cell Gene Expression|16|12|R1: Barcode-UMI; R2: Transcript|scrna_10x_v3|
|scRNA-seq|Drop-seq|12|8|R1: Barcode-UMI; R2: Transcript|dropseq|
|ChIA-Drop|ChIA-Drop|16|N/A|R1: Barcode-DNA; R2: DNA|chiadrop|
|scATAC-seq|10× Genomics (v2) Single Cell ATAC|16|N/A|R1: DNA; R2: Barcode; R3: DNA|scatac_10x_v1|
|SPRITE and derivates|SPRITE|4 kinds with 5 rounds|N/A|R1: DPM-DNA; R2: Y-ODD-EVEN-ODD|sprite|
|SPRITE and derivates|RD-SPRITE|5 kinds with 5 rounds (label DNA with DPM, label RNA with RPM)|N/A|R1: DNA/RNA; R2: Y-ODD-EVEN-ODD-DPM/RPM (DPM: DNA, RPM: RNA)|rdsprite|
|SPRITE and derivates|Single-cell SPRITE|4 kinds with 6 rounds (label cells with first last 3 barcodes in Read 2, label complex with all six barcodes)|N/A|R1: DNA; R2: Y-EVEN-ODD-EVEN-ODD-DPM|scsprite|
|scARC-seq|10× Genomics Single Cell Multiome ATAC + Gene Expression|Gene Expression library: 16, ATAC library: 16|Gene Expression library: 12; ATAC library: N/A|R1: DNA (R1 in ATAC); R2: Barcode (R2 in ATAC); R3: DNA (R3 in ATAC); R4: Barcode-UMI (R1 in GEX); R5: Transcript (R2 in GEX)|scarc_10x_v1|



### 5.2 Easy process

Copy the corresponding process pipeline to the directory you want to store the result, for example: ChIA-Drop.
```
(ScSmOP) $ ~/ScSmOP/scsmop.sh  -t chiadrop -n CHDP -1 CHDP_R1.fastq.gz -2 CHDP_R2.fastq.gz -b ~/RefGenome/bwa_dm3_ref/dm3.fa -s dm3.size.txt
```
Different experiment require different parameters:
```
(ScSmOP) $ ~/ScSmOP/scsmop.sh 

Usage: ./scsmop.sh [options] -t [STR] -n [STR] -1 [FILE] -2 [FILE] -r [DIR]
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

   
```

## 6 Other Information in ScSmOP Wiki
* [Case Study](https://github.com/ZhengmzLab/ScSmOP/wiki/Case-Study)
* [ScSmOP Module Configuration](https://github.com/ZhengmzLab/ScSmOP/wiki/ScSmOP-Module-Configuration)
* [ScSmOP Standard Output](https://github.com/ZhengmzLab/ScSmOP/wiki/ScSmOP-Standard-Output)
* [Case Study of DIY](https://github.com/ZhengmzLab/ScSmOP/wiki/Case-Study#11-diy)


