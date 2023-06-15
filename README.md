# ScSmOP

Single cell Single Molecule Multiple Omics Pipeline.

-----------

## Summary

ScSmOP is a universal pipeline capable of performing data processing for a wide range of single-cell single-molecule omics like scRNA-seq, scATAC-seq, ChIA-Drop, SPRITE and its derivates, Drop-seq and some other techniques which based on barcode and UMI. It have pre-prepared several pipelines for current popular which can be easily launched with some simple editions.

## Installation

First, download ScSmOP to your computer:

```bash
# Get ScSmOP source using git
git clone https://github.com/kasenjing/ScSmOP.git
cd ScSmOP

# Alternatively, get latest ScSmOP source from releases
wget https://github.com/kasenjing/ScSmOP/archive/refs/tags/v0.1.3.tar.gz
tar -zxvf v0.1.3.tar.gz
cd ScSmOP-0.1.3
```

### Advanced user

If you are familiar with linux operate system and conda, after activate conda, you can install ScSmOP using the one-click `conda_local_install_ScSmOP_dependencies.sh`

    bash conda_local_install_ScSmOP_dependencies.sh

It will automatically create a conda environment `ScSmOP` and install all the softwares ScSmOP required.

### New to linux and conda

For users who are not familiar with linux operate system or conda, we prepared one-click `local_install_ScSmOP_dependencies.sh` to install ScSmOP:

    bash local_install_ScSmOP_dependencies.sh

It will create conda environment as well as install all the softwares ScSmOP required.

## Run ScSmOP on different experiment data

We've prepared several pipelines for easy data process, if your experiment meet the listed experiment, you can process your data as in [Easy Process](#easy-process).

<html>
 <body>
  <table>
   <tr>
    <td t="s" id="sjs-A1" white-space: nowrap><b>Omics<b></td>
    <td t="s" id="sjs-B1"><b>Experiment<b></td>
    <td t="s" id="sjs-C1"><b>Barcode(s) length<b></td>
    <td t="s" id="sjs-D1"><b>UMI length<b></td>
    <td t="s" id="sjs-E1"><b>Pipeline script<b></td>
    <td t="s" id="sjs-F1"><b>Configuration file<b></td>
    <td t="s" id="sjs-G1"><b>Whitelist source<b></td>
   </tr>
   <tr>
    <td rowspan="3" t="s" id="sjs-A2" white-space: nowrap>scRNA-seq</td>
    <td t="s" id="sjs-B2">10× Genomics (v2) Single Cell Gene Expression</td>
    <td t="n" id="sjs-C2">16</td>
    <td t="n" id="sjs-D2">10</td>
    <td t="s" id="sjs-E2">10x_scrna-v2.sh</td>
    <td t="s" id="sjs-F2">10x_rna_v2_config.json</td>
    <td t="s" id="sjs-G2">https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt</td>
   </tr>
   <tr>
    <td t="s" id="sjs-B3">10× Genomics (v3) Single Cell Gene Expression</td>
    <td t="n" id="sjs-C3">16</td>
    <td t="n" id="sjs-D3">12</td>
    <td t="s" id="sjs-E3">scRNA.sh</td>
    <td t="s" id="sjs-F3">10x_rna_v3_config.json</td>
    <td t="s" id="sjs-G3">https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz</td>
   </tr>
   <tr>
    <td t="s" id="sjs-B4" >Drop-seq</td>
    <td t="n" id="sjs-C4">12</td>
    <td t="n" id="sjs-D4">8</td>
    <td t="s" id="sjs-E4">dropseq_scrna.sh</td>
    <td t="s" id="sjs-F4">dropseq_rna_config.json</td>
    <td t="z" id="sjs-G4"></td>
   </tr>
   <tr>
    <td t="s" id="sjs-A5" white-space: nowrap>ChIA-Drop</td>
    <td t="s" id="sjs-B5">ChIA-Drop</td>
    <td t="n" id="sjs-C5">16</td>
    <td t="s" id="sjs-D5">N/A</td>
    <td t="s" id="sjs-E5">ChIA-Drop.sh</td>
    <td t="s" id="sjs-F5">chiadrop_config.json</td>
    <td t="s" id="sjs-G5">longranger-2.2.2/longranger-cs/2.2.2/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt-</td>
   </tr>
   <tr>
    <td t="s" id="sjs-A6" white-space: nowrap>scATAC-seq</td>
    <td t="s" id="sjs-B6">10× Genomics (v2) Single Cell ATAC</td>
    <td t="n" id="sjs-C6">16</td>
    <td t="s" id="sjs-D6">N/A</td>
    <td t="s" id="sjs-E6">scATAC.sh</td>
    <td t="s" id="sjs-F6">10x_atac_v2_config.json</td>
    <td t="s" id="sjs-G6">cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz</td>
   </tr>
   <tr>
    <td rowspan="3" t="s" id="sjs-A7" white-space: nowrap>SPRITE and derivates</td>
    <td t="s" id="sjs-B7">SPRITE</td>
    <td t="s" id="sjs-C7">4 kinds with 6 rounds</td>
    <td t="s" id="sjs-D7">N/A</td>
    <td t="s" id="sjs-E7">SPRITE.sh</td>
    <td t="s" id="sjs-F7">sprite_config.json</td>
    <td t="s" id="sjs-G7">GSE114242: GSE114242_human_config.txt.gz</td>
   </tr>
   <tr>
    <td t="s" id="sjs-B8">RD-SPRITE</td>
    <td t="s" id="sjs-C8">5 kinds with 5 rounds (label DNA with DPM, label RNA with RPM)</td>
    <td t="s" id="sjs-D8">N/A</td>
    <td t="s" id="sjs-E8">rdSPRITE.sh</td>
    <td t="s" id="sjs-F8">rdsprite_config.json</td>
    <td t="s" id="sjs-G8">https://github.com/GuttmanLab/sprite2.0-pipeline/blob/master/config.txt</td>
   </tr>
   <tr>
    <td t="s" id="sjs-B9">Single-cell SPRITE</td>
    <td t="s" id="sjs-C9">4 kinds with 6 rounds (label cells with first 3 barcodes, label complex with all six barcodes)</td>
    <td t="s" id="sjs-D9">N/A</td>
    <td t="s" id="sjs-E9">scSPRITE.sh</td>
    <td t="s" id="sjs-F9">scsprite_config.json</td>
    <td t="s" id="sjs-G9">https://github.com/caltech-bioinformatics-resource-center/Guttman_Ismagilov_Labs/blob/master/scSPRITE/misc/config_dpm6_y-stag_scSPRITE.txt</td>
   </tr>
   <tr>
    <td t="s" id="sjs-A10" white-space: nowrap>scARC-seq</td>
    <td t="s" id="sjs-B10">10× Genomics Single Cell Multiome ATAC + Gene Expression</td>
    <td t="s" id="sjs-C10">Gene Expression library: 16, ATAC library: 16</td>
    <td t="s" id="sjs-D10">Gene Expression library: 12; ATAC library: N/A</td>
    <td t="s" id="sjs-E10">scarc.sh</td>
    <td t="s" id="sjs-F10">scarc_rna_config.json; scarc_atac_config.json</td>
    <td t="s" id="sjs-G10">737K-arc-v1-scrna.txt:cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz; 737K-arc-v1-scatac.txt:cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz</td>
   </tr>
  </table>
 </body>
</html>

### Easy process

Copy the corresponding process pipeline to the directory you want to store the result, for example: ChIA-Drop.

    cp /path/to/ScSmOP/PipelineScript/ChIA-Drop.sh .

Edit the parameters in the pipeline script including library name, input fastq, reference genome, directory of BARP.

    vi ChIA-Drop.sh

Things you need to edit:
Category|Description
--------|-----------
LIB_NAME | Library name
ChIADrop_FILES | "-1 Sample_L001_R1,Sample_L002_R1 -2 Sample_L001_R2,Sample_L002_R2" If you have multiple read, separate them with "," like the example. There can be at most four kinds of FASTQs from Read 1 (R1) to Read 4 (R4) specified with -1 to -4 respectively.
THREAD | Thread to use.
ChIADrop_CONFIG | If you want to use your own configuration file instead of default, specifies the directory to the config file.
STAR_REF_GENOME | STAR genome index directory.
BWA_REF_GENOME | bwa genome index directory (same as bwa).

Then you can run the pipeline:

    ./ChIA-Drop.sh

## Customize for your own barcoding method

Please refer to the paper or the wiki.
