#!/bin/bash
## perform local installation of dependencies for ScSmOP


echo "Unzip the configuration files..."
cd BarcodeBucket
gunzip *
cd ../ConfigFiles
gunzip *
cd ../PipelineScript
chmod +x *

cd ../Tools
echo "installing trimgalore..."
chmod +x TrimGalore-0.6.6/trim_galore 
ln -s TrimGalore-0.6.6/trim_galore trim_galore 

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ./miniconda3
ln -s miniconda3/bin/conda conda
which conda
CondaDir="$(pwd)/miniconda3/"
source miniconda3/bin/activate
conda config --add channels conda-forge
conda create -n ScSmOP python=3.10 -y -c conda-forge
ScSmOPCondaDir="$(pwd)/miniconda3/envs/ScSmOP"

source ${CondaDir}/bin/activate ScSmOP

echo "installing gcc"

conda config --add channels conda-forge
conda config --add channels bioconda


which conda 
if [[ ! -f $GCCDir ]]
then
	conda install -c "conda-forge/label/cf202003" cxx-compiler -y -n ScSmOP
	ln -s ${ScSmOPCondaDir}/bin/x86_64-conda_cos6-linux-gnu-gcc ${ScSmOPCondaDir}/bin/gcc
	ln -s ${ScSmOPCondaDir}/bin/x86_64-conda_cos6-linux-gnu-g++ ${ScSmOPCondaDir}/bin/g++
	ln -s ${ScSmOPCondaDir}/bin/x86_64-conda_cos6-linux-gnu-ar ${ScSmOPCondaDir}/bin/ar
fi


echo "installing make"
if [[ ! -f $MakeDir ]]
then
	conda install make -n ScSmOP -y -c conda-forge
fi

cd ../BARP-prog
chmod +x third_party/htslib-1.16/configure
chmod +x third_party/zlib-1.2.11/configure
make
cd ../Tools
ln -s ../BARP-prog/BARP BARP
ln -s ${ScSmOPCondaDir}/bin/python python3
# install modified STAR 
tar -zxvf STARM.tar.gz
cd STARM/source
chmod +x htslib/configure
chmod +x zlib/configure
make STAR
cd ../..
mv STARM/source/STAR .
rm -r STARM


# install samtools if no samtools found in PATH
conda install samtools=1.18 -n ScSmOP -y 
ln -s ${ScSmOPCondaDir}/bin/samtools samtools 

# install bedtools if no bedtools found in PATH
conda install bedtools=2.31.0 -n ScSmOP -y 
ln -s ${ScSmOPCondaDir}/bin/bedtools bedtools

conda install macs2=2.2.9.1 -n ScSmOP -y 
ln -s ${ScSmOPCondaDir}/bin/macs2 macs2

${ScSmOPCondaDir}/bin/pip install numpy
${ScSmOPCondaDir}/bin/pip install pandas

conda install cutadapt=4.5 -n ScSmOP -y 
ln -s ${ScSmOPCondaDir}/bin/cutadapt cutadapt

conda install pairtools=1.0.2 -n ScSmOP -y 
${ScSmOPCondaDir}/bin/pip install pairtools

conda install pigz=2.4 -n ScSmOP -y
ln -s ${ScSmOPCondaDir}/bin/pigz pigz

# install bwa if no bwa found in PATH
conda install bwa=0.7.17 -n ScSmOP -y 
ln -s ${ScSmOPCondaDir}/bin/bwa bwa

# install Rscript
conda install r-base=4.3.1 -n ScSmOP -y -c conda-forge
ln -s ${ScSmOPCondaDir}/bin/Rscript Rscript 

# install r packages used by scsmop
conda install -c bioconda -n ScSmOP -y bioconductor-dropletutils=1.20.0
conda install -c r -n ScSmOP -y r-ggplot2=3.4.4
conda install -c r -n ScSmOP -y r-dplyr=1.1.4


cd ..
ScSmOP_dir=$(pwd | sed 's/\//\\\//g')

sed -i "s/REF_OF_ScSmOP/${ScSmOP_dir}/g" PipelineScript/*
sed -i "s/REF_OF_ScSmOP/${ScSmOP_dir}/g" scsmop.sh
chmod +x scsmop.sh
# check installation 
InstallSuccesss=1
for tool in samtools macs2 bedtools Rscript bwa
do
    if [[ ! -f Tools/${tool} ]]
    then
        InstallSuccesss=0
        echo -e "Installation of ${tool} failed, you may need to reinstall them with \n\n\tconda install -n ScSmOP ${tool}\n\nor you can link one of your previously installed ${tool} to the directory of ScSmOP/Tools"
    fi
done

for tool in BARP STAR
do
    if [[ ! -f Tools/${tool} ]]
    then
        InstallSuccesss=0
        echo -e "Installation of ${tool} failed, it may caused by the failed installation of \"gcc\", \"g++\", reinstall these with \n\tconda install -c "conda-forge/label/cf202003" cxx-compiler -y -n ScSmOP\nand then recompile it."
        if [[ ${tool} == "BARP" ]]
        then
            echo -e "\n\tcd BARP-prog \n\t make \n\tcd .. \n\t ln -s BARP-prog/BARP Tools/BARP\n"
        elif [[ ${tool} == "STAR" ]]
        then
            echo -e "\n\tcd Tools \n\ttar -zxvf STARM.tar.gz \n\tcd STARM/source \n\tmake STAR\n\tcd ../..\n\tmv STARM/source/STAR . \n\trm -r STARM\n"
        else
            echo ""
        fi
    fi
done

if [[ $InstallSuccesss == 0 ]]
then
    echo -e "Installation failed. Check the information above to complete the installation or re-install."
else
    echo -e "Installation succeeded."
fi

echo "Setup finished."
