# ScSmOP
ScSmOP, a universal pipeline for barcode-indexed single-cell single-molecule multiple omics data analysis. 

## Installation 

    git clone https://github.com/Zhengmzlab/ScSmOP.git
    unzip ScSmOP
    cd ScSmOP
    bash local_install_ScSmOP_dependencies.sh

## Run BARP on different experiment data
Copy the corresponding process pipeline to the directory you want to store the result, for example: ChIA-Drop.

    cp /path/to/ScSmOP/PipelineScript/ChIA-Drop.sh .
    
Edit the parameters in the pipeline script including library name, input fastq, reference genome, directory of BARP.
    
    vi ChIA-Drop.sh

Things you need to edit:
Category|Description
--------|-----------
LIB_NAME | Library name
ChIADrop_FILES | "-1 Read_1.1,Read_1.2 -2 Read_2.1,Read_2.2" If you have multiple read, separate them with "," like the example. There can be at most four FASTQs from Read 1 to Read 4 specified with -1 to -4 respectively.
THREAD | Thread to use.
ChIADrop_CONFIG | If you want to use your own configuration file instead of default, specifies the directory to the config file.
STAR_REF_GENOME | STAR genome index.
BWA_REF_GENOME | bwa genome index.

Then you can run the pipeline:

    ./ChIA-Drop.sh

## Customize for your own barcoding method
Please refer to the paper or the [wiki](https://github.com/Zhengmzlab/ScSmOP/wiki).
## Understand output.
Please refer to the [wiki](https://github.com/Zhengmzlab/ScSmOP/wiki).
