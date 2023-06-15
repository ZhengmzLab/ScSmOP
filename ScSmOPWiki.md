Welcome to the ScSmOP wiki!

# What is ScSmOP?

**S**ingle-**c**ell **S**ingle-**m**olecule **m**ultiple **O**mics **P**ipeline is a general purpose pipeline, which works for processing the "one barcode one group" type of state-of-art multi-omics data.<br>Currently, ScSmOP support for data processing of
 - ChIA-Drop, scRNA-seq, scATAC-seq, scARC-seq (10x genomics)
 - SPRITE, scSPRITE, rdSPRITE (split-pool) data.<br>
 - DIY for your data.

## Quick start
ScSmOP have prepared configuration files for techniques mentioned above using corresponding barcode whitelist:
|Technique|Configuration file name|Barcode whitelist used|Source|
|---------|--------|--------|---|
|ChIA-Drop|chia-drop_config.json|4M-with-alts-february-2016.txt|longranger-2.2.2/longranger-cs/2.2.2/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt|
|SPRITE|sprite_config.json|SPRITE.EVEN.txt;SPRITE.Y.txt;SPRITE.DPM.txt;SPRITE.ODD.txt|GSE114242: GSE114242_human_config.txt.gz|
|scSPRITE|scsprite_config.json|scSPRITE_Y_Barcode.txt;scSPRITE_ODD_Barcode.txt;scSPRITE_EVEN_Barcode.txt;scSPRITE_DPM_Barcode.txt;scSPRITE_DPMPRE_Barcode.txt|github.com/caltech-bioinformatics-resource-center/Guttman_Ismagilov_Labs/blob/master/scSPRITE/misc/config_dpm6_y-stag_scSPRITE.txt|
|rdSPRITE|rdsprite_config.json|rdSPRITE.EVEN.txt;rdSPRITE.Y.txt;rdSPRITE.DPM.txt;rdSPRITE.ODD.txt;rdSPRITE.RPM.txt|github.com/GuttmanLab/sprite2.0-pipeline/blob/master/config.txt|
|scATAC-seq|scatac_config.json|737K-cratac-v1.txt|cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt.gz|
|scRNA-seq|scrna_config.json|3M-february-2018.txt.gz|cellranger-7.0.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz|
|scARC-seq|scarc-rna_config.json; scarc-atac_config.json|737K-arc-v1-scrna.txt;737K-arc-v1-scatac.txt|737K-arc-v1-scrna.txt:cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz; 737K-arc-v1-scatac.txt:cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt.gz|
|Spatial|scRNA-spatial_config.json|visium-v2.txt;visium-v2_coordinates.txt|spaceranger-2.0.0/lib/python/cellranger/barcodes/visium-v2.txt; spaceranger-2.0.0/lib/python/cellranger/barcodes/visium-v2_coordinates.txt|

*Notice that different assay kit will incorporate different sets of barcodes whitelist to the libary, you can generate your own configuration file using the sets of whitelist corresponding to your experiment as in [Generate configuration file for your own data](#generate-configuration-file-for-your-own-data)*

### Run ScSmOP
To run ScSmOP on your ChIA-Drop data:

```
    # First copy ChIA-Drop.sh to your processing directory.
    cp ScSmOP/PipelineScript/ChIA-Drop.sh .
```

Then edit the `ChIA-Drop.sh` for your data:
```
    LIB_NAME="Example"  # Name of your library
    ChIADrop_FILES="-1 R1.fq -2 R2.fq" # FASTQ files of your ChIA-Drop library, if you have multiple Read 1 and Read 2s, separate them with ,
    THREAD=10   # Thread to use
    ChIADrop_CONFIG="NewChIA-Drop_config.json" # Leave it as "-" if your data using the barcode whitelist "4M-with-alts-february-2016.txt", else generate new config file like above.
    BWA_REF_GENOME="/path/to/bwa/reference" # bwa reference genome directory.

    # The rest of the file do not need to edit.
    BARP_DIR="/mnt/e/FinalTest/ScSmOP-main"

    echo "Start: $(date)" > TimeStamp
    if [ ${BARP_DIR} == "-" ]
    then
        echo -e "BARP directory is necessary, but not exist."
        exit 1
    fi

    if [ ${ChIADrop_CONFIG} == "-" ]
    then
        ${BARP_DIR}/PipelineScript/BarcodeIdentification.sh -t chiadrop -n ${LIB_NAME} -p ${BARP_DIR} ${ChIADrop_FILES} -@ ${THREAD} -c ${ChIADrop_CONFIG}
    else
        ${BARP_DIR}/PipelineScript/BarcodeIdentification.sh -t chiadrop -n ${LIB_NAME} -p ${BARP_DIR} ${ChIADrop_FILES} -@ ${THREAD} -c ${ChIADrop_CONFIG}
    fi

    if [ -f BarcodeIdentification.done ]
    then
        ${BARP_DIR}/PipelineScript/SequenceAlignment.sh -t chiadrop -n ${LIB_NAME} -p ${BARP_DIR} -f ${BWA_REF_GENOME} -@ ${THREAD} -k "DNA"
    fi

    if [ -f SequenceAlignment.done ]
    then
        ${BARP_DIR}/PipelineScript/GroupAndDataRefine.sh -t chiadrop -n ${LIB_NAME} -p ${BARP_DIR} -d 02.ReadAlign/CHDP_DNA.bam -@ ${THREAD} -x dm3.size.txt
    fi

    if [ -f GroupAndDataRefine.done ]
    then
        ${BARP_DIR}/PipelineScript/QualityAssessment.sh -t chiadrop -n ${LIB_NAME} -p ${BARP_DIR}
    fi

    echo "End: $(date)" >> TimeStamp
```
*Note that the NewChIA-Drop_config.json need to exist under ScSmOP/ConfigFiles or specifies it with absolute directory.*
Run `ChIA-Drop.sh`

```
    ./ChIA-Drop.sh
```
### Generate configuration file for your own data
First edit the OriginalConfigFile.json corresponding to you experiment: e.g. OriginalConfigFile_ChIA-Drop.json
```
    {
        "barcode chain" : [ {"BC-GENOMEA|1": "R1:1"}, {"GENOMEA|2":"R2:1"}],
        "identifier" : [{"COMPLEX":"BC"}],
        "barcode type" : 
        {
            "GENOMEA|1":
            {
                "SPACE": 8
            },
            "BC": 
            {
                "SPACE": 0,
                "LAXITY": 0,
                "LENGTH": "16",
                "DENSE":1,
                "MISMATCH": 1,  
                "WHITE LIST":"../BarcodeBucket/4M-with-alts-february-2016.txt" # barcode whitelist file directory.
            }
        }
    }
```    
Then generate the actual configuration file used by BARP
```
    python3 ScSmOP/PythonScript/GenerateConfigFile.py -i OriginalConfigFile_ChIA-Drop.json -o NewChIA-Drop
```    
This will generate configuration file: *NewChIA-Drop_config.json*

### Output of ScSmOP

There are four directories storing the processed files by ScSmOP:

- 01.BarcodeIden
  - Results of barcode identification
    - [Result FASTQs](#resulting-fastqs) with idenfified barcodes and identifier information 
    - [Barcoding statistics](#barcoding-statistics)
    - [Genomic material statistics](#genomic-material-statistics)
    - [Library idb files](#idb-files)
- 02.ReadAlign
  - [Results of read alignment](#result-of-read-alignment)
- 03.GroupAndRefine
  - Results of barcode grouping and data refinement:
    - [Alignment filteration](#filtered-alignment)
    - [Barcode grouping](#cluster-and-rdpcluster)
    - [Removing duplication/Merging fragment](#cluster-and-rdpcluster)
    - [Further refinement](#raw_matrix)
    - [Final results](#filtered_matrix)
- 04.QualityAssess
  - Result of library statistic
    - Statistic table
    - Statistic plot if have

### File format explanation

#### Resulting FASTQs
There were three possible kinds of output FASTQs: DNA, RNA, NDNR, each kinds of FASTQs have *n* read corresponding to the input FASTQs from Read 1 to Read 4.
- DNA FASTQs: Reads sorted to DNA based on identified barcode.
- RNA FASTQs: Reads sorted to RNA based on identified barcode.
- NDNR FASTQs: Reads can not sorted to DNA or RNA based on identifed barcode.

The format of FASTQ is just like the original FASTQ format with identified barcode and identifiers store in the name each record.
```
    @SRR7216005.man.87|||DPM6D6|NYBot17_Stg|Odd2Bo30|Even2Bo22|Odd2Bo50|||COMPLEX_1
    TTTGAGGGCTGACTCTTTACTTGCACTGTCCTAGGTAGGAAGTGGAAGATAAAAAA
    +
    EEEEEEEEEEEAEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
```
```
    @SRR7722051.1.1|||BC:AGTACCACAAGAAGAG|||COMPLEX_1
    ATTGGGGTATAAAATGAGAAAATTTTAAGTCGATTTACAAGTGTAAAAAAAATGTCAAAAAATATCACCTTATTTTTCGGAAGTGTGGGCGTGACAGTTTTGTGCGGCACGGAAGAGCACACGNCTG
    +
    <-7<J<AJAFAFFJJJ-A7--<-<AFJA--<FAFAFFFJ<F<AFFF<F-7F<FAJFAAFFFA<JJAJJJJFJ<F<F7FF<<-)7-7-))7)--FFJ<--A7FF-A)-))--)7-<-<F<)A77#F<7
```
```
    @SRR12212044.sra.3|||NYbotLigEven_D12_Stg|Even2Bo70|Odd2Bo1|Even2Bo75|Odd2Bo88|DPMPRE|DPM6bot91|||CELL_1|COMPLEX_1
    NTGCCAGGGGGTTGAATGTCTTTTTCCTTTTCTTACTAAGAATATAGTACTTGACAACACGCTGCCATTAGGAAGAAGAAAATAATCTTACGAGAAGAAA
    +
    #FFFFFFFF,,FFFFFFFFF::FFFFFFFFFFFF:FF:FFF::FFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
The name is separated to 3 fields with `|||`:
|Field|Content|Description|
|----|----|----|
|1<sup>st<sup>|Sequence name|Only retain the name before the first space appear the original FASTQ|
|2<sup>nd<sup>|Identified barcode|Barcodes are separated with `|`, each barcode can be name or "Barcode type":"Corresponding barcode in whitelist"; if no corresponding White barcode, "NOT_FOUND" is write.|
|3<sup>rd<sup>|Identifier|Identifiers are separated with `|`, each identifier is specified with "Type"_"ID", the "ID" of each identifier is independent.

#### Barcoding statistics
```
    *** $ head BarcodeIdentification_statistic.txt ***
    COMPLEX : 7132895
    [DPM] [GENOMEA|1] | Read Count : 42888991
    [Not_Found] [GENOMEA|1] | Read Count : 1528760
    [DPM] [Not_Found] | Read Count : 0
    [Not_Found] [Not_Found] | Read Count : 0
    [Y] [ODD] [EVEN] [ODD] | Read Count : 28032533
    [Y] [ODD] [Not_Found] [Not_Found] | Read Count : 7792381
    [Y] [ODD] [EVEN] [Not_Found] | Read Count : 3383577
    [Y] [Not_Found] [Not_Found] [Not_Found] | Read Count : 2404195
    [Not_Found] [ODD] [EVEN] [ODD] | Read Count : 929331
```
The first several lines specifies the unique group count in the processed library.

The next lines specifies the read count of each barcode permutation in the processed library, the read count should be modified to recover real situation according to [Barcode identification algothrim](#barcode-identification).

#### Genomic material statistics
```
    *** $ head stat ***
    Total reads:44417751
    Total DNA reads:0
    Fully barcoded DNA reads:0
    Not fully barcoded DNA reads:0
    Total RNA reads:0
    Fully barcoded RNA reads:0
    Not fully barcoded RNA reads:0
    Total non-defined reads:44417751
    Fully barcoded non-defined reads:27219050
    Not fully barcoded non-defined reads:17198701
```
Reads sorted to DNA, RNA or NDNR and fully barcoded count in each categories.
#### IDB files
Please refer to the supplementary note of the paper for detailed IDB value.
```
    *** $ head rdSPRITE.COMPLEX.idb ***
    --- IDB value - IDB ID - IDB count in the library ---
    83886080        3639261 2
    138414081       1243089 1
    25169922        134200  1
    125833218       5075147 1
```
#### Result of Read alignment
|File name|Description|Appear in|
|---|---|---|
|LIB_DNA.bam|DNA alignment|ChIA-Drop, SPRITE, rdSPRITE, scATAC-seq|
|LIB.Barcoded.Aligned.out.bam|Fully barcoded DNA alignment|scSPRITE|
|LIBAligned.sortedByCoord.out.bam|Fully barcode RNA alignment|scRNA-seq, Spatial-RNA-seq, rdSPRITE|
|LIBSignal.Unique.str1.out.bg; LIBSignal.Unique.str2.out.bg|RNA Coverage|scRNA-seq, Spatial-RNA-seq, rdSPRITE|
#### Filtered alignment
Alignments are filtered to retain uniquely mapped reads.
|File name|Filter with samtools|Description|Appear in|
|---|---|---|---|
|LIB.F2304.q30.bam| -F 2304 -q 30| Uniquely aligned DNA with valid barcode|ChIA-Drop|
|LIB.PrimaryAlign.out.bam| -F 0x100|One alignment for one read|scRNA-seq, scRNA-seq in scARC-seq|
|LIB_DNA.F2304.q30.bam| -F 2304 -q 30| Uniquely aligned DNA with valid barcode|scATAC-seq, scATAC-seq in scARC-seq|
|LIB.Barcoded.Aligned.out.bam| -q 255 | Uniquely aligned DNA with full set of barcode|scSPRITE|
|LIB.Barcoded.UniqAlign.DNA.bam|-F 2304 -q 30|Uniquely aligned DNA with full set of barcode|rdSPRITE,SPRITE|
|LIB.Barcoded.UniqAlign.RNA.bam|-q 255|Uniquely aligned RNA with full set of barcode|rdSPRITE|
|LIB_P2S.bam|-|One alignment for one read with barcode and UMI information stored in bam field: CB and UR respectively|scATAC-seq, scATAC-seq in scARC-seq|
#### cluster and RDP.cluster
```
    *** $ head rdSPRITE.DNA.cluster ***
    #This file is produced by barp p2s.
    @CN     COMPLEX Count
    COMPLEX_4196353 2   chr2    96439353    96439478    chr2    96445065    96445120
    COMPLEX_4200451 2   chr1    85968959    85969031    chr1    85968959    85969031
    COMPLEX_8196    5   chr2    157832187   157832245   chr2    157883432   157883582   chr2    157832187   157832245   chr2    157832187    157832245        chr2    157883432 157883582
    COMPLEX_4204549 1   chr8    14266475    14266539
    COMPLEX_4208647 1   chr7    27303335    27303459
```
```
    *** $ head rdSPRITE.DNA.RDP.cluster ***
    # This file is produced by barp rdp
    @CN     complex COUNT
    COMPLEX_4196353 2   chr2    96439353    96439478    chr2    96445065    96445120    1   1
    COMPLEX_4200451 1   chr1    85968959    85969031    2
    COMPLEX_8196    2   chr2    157832187   157832245   chr2    157883432   157883582   3   2
    COMPLEX_4204549 1   chr8    14266475    14266539    1
    COMPLEX_4208647 1   chr7    27303335    27303459    1
```
```
    *** $ head scSPRITE.DNA.RDP.cluster ***
    # This file is produced by barp rdp
    @CN     cell    complex COUNT
    CELL_926    COMPLEX_2049    1   chr16   38350928    38350960    2
    CELL_778    COMPLEX_4098    1   chr8    9041233     9041287     1
    CELL_1192   COMPLEX_6147    1   chr10   17566046    17566141    1
    CELL_157    COMPLEX_8196    8   chr2    98667125    98667213    chr1    3037492 3037582 chr2    98667125    98667215    chr2    98666550    98666640    chr3    7058394 7058484 chrUn_GL456383  28713   28802   chr2    98666874    98666964    chr2    98662870    98662960    1   1    7   1   1   1   1   1
```
**cluster** stores the deconvolved biological groups with *one-line-one-group* fanshion as we called [***inline BED format***](#cluster-and-rdpcluster). The format split each line as 4 fields that the field are not splited by certain symbols, instead these fields are organized based the count field. 
**RDP.cluster** stores the deduplicated groups, the duplicates times are stored in the forth field for each fragment.
|Field|Content|Description|
|---|-----|---|
|1<sup>st<sup>|Group field|This field contains columns before "Count" field, specifying the biological group of all the fragments within the line.|
|2<sup>nd<sup>|Count field|This field has only one column specifies the fragment count within the group.|
|3<sup>rd<sup>|Fragments field| Fragments within the group that each fragment is specified with chromosome start end.|
|4<sup>th<sup>|Additional field| This field stores additional information for each fragment, if there are multiple additional information for each fragment, the information should be recorded separately.|

 An example of [***inline BED format***](#cluster-and-rdpcluster) with two addition information for each fragments: InfoA and InfoB.
```
    COMPLEX_8196    2   chr2    157832187   157832245   chr2    157883432   157883582   InfoA   InfoA   InfoB   InfoB
```
#### Merged cluster
**Merged_cluster.txt** stores the information of merged fragments in ChIA-Drop with [***inline BED format***](#cluster-and-rdpcluster).
#### Qualified fragments in scATAC-seq
[**LIB.Qualified.bed**](#qualified-fragments-in-scatac-seq) converted the [***inline BED format***](#cluster-and-rdpcluster) RDP.cluster to bed format that the columns corresponding with: Fragment chromosome | start | end | Cell ID | Duplicate time. 
```
    *** $ head PBMC.Qualified.bed ***
    chr10   116264702   116264748   CELL_2049   1
    chr18   23631039    23631090    CELL_4098   2
    chr5    175075551   175075637   CELL_6147   1
```
[**LIB.Qualified.Sorted.bed**](#qualified-fragments-in-scatac-seq) sorting the fragments from **LIB.Qualified.bed** according to genomic coordinates.
```
    *** $ head PBMC.Qualified.Sorted.bed ***
    chr1    10074   10322   CELL_63 1
    chr1    10098   10334   CELL_220        1
    chr1    10152   10438   CELL_118        1
    chr1    10169   10341   CELL_185        1
    chr1    10229   10304   CELL_550        1
```
[**LIB.QualifiedFragments.bedgraph**](#qualified-fragments-in-scatac-seq) generated coverage from **LIB.Qualified.Sorted.bed** using *bedtools genomecov*.
```
    *** $ head PBMC.QualifiedFragments.bedgraph ***
    chr1    10074   10098   1
    chr1    10098   10152   2
    chr1    10152   10169   3
    chr1    10169   10229   4
    chr1    10229   10248   5
```
[**LIB.narrowPeaks**](#qualified-fragments-in-scatac-seq) called peaks from **LIB.QualifiedFragments.bedgraph** using *macs2 bdgpeakcall*.
```
    *** $ head PBMC.narrowPeaks ***
    track type=narrowPeak name="PBMC.narrowPeaks" description="PBMC.narrowPeaks" nextItemButton=on
    chr1    778240  779288  PBMC.narrowPeaks_narrowPeak1    1550    .       0       0       0       466
    chr1    817218  817486  PBMC.narrowPeaks_narrowPeak2    320     .       0       0       0       149
    chr1    826703  826928  PBMC.narrowPeaks_narrowPeak3    70      .       0       0       0       28
    chr1    827068  827852  PBMC.narrowPeaks_narrowPeak4    1240    .       0       0       0       523
```
[**LIB.PeakCellCount.bed**](#qualified-fragments-in-scatac-seq) Fragment counts fall into each called peak for each cell.
```
    *** $ head PBMC.PeakCellCount.bed ***
    -- Peak coordinate --  -- Cell ID --    -- Fragment count --
    chr1    778240  779288  CELL_187        2
    chr1    778240  779288  CELL_154        2
    chr1    778240  779288  CELL_275        2
    chr1    778240  779288  CELL_38 4
    chr1    778240  779288  CELL_341        2
```
#### SubGEM
**SubGEM** is generated for ChIA-Drop type data with [***inline BED format***](#cluster-and-rdpcluster), complexes are splited into subGEMs by chromosome. 
```
    *** $ head CHDP.ChIADrop.SubGEM ***
    --- The complex are named as "COMPLEX_ID-nth SubGEM in the complex" ---
    COMPLEX_2049-0  1   chr2    12713870    12714100
    COMPLEX_2049-1  1   chr2R   7863021     7863695
    COMPLEX_2049-2  1   chr3L   11185707    11185889
    COMPLEX_2049-3  1   chr3R   3018595     3018777
    COMPLEX_4098-0  2   chr2R   6408379     6408606     chr2R   12027465    12027974
    COMPLEX_4098-1  4   chr3L   1490949     1491382     chr3L   4460139     4460366     chr3L   21232289    21232501    chr3L   22314225    22314879
    COMPLEX_4098-2  5   chr3R   3671537     3671764     chr3R   5063085     5063312     chr3R   9839373     9839600     chr3R   19349478    19349626    chr3R   25367853    25368080
    COMPLEX_4098-3  2   chrX    9123872     9124094     chrX    17155881    17156108
```
#### rgn file
**rgn file** is generated for visualized by ChIA-View with one line one fragment, the forth line in the **.rgn** files is the unique ID of the fragment.
```
    *** $ head rdSPRITE.DNA.rgn ***
    --- For data do not require SubGEM, the complexes are named as "COMPLEX_ID-Total fragment in the complex-nth fragment" ---
    chr2    96439353        96439478        COMPLEX_4196353-2-0
    chr2    96445065        96445120        COMPLEX_4196353-2-1
    chr1    85968959        85969031        COMPLEX_4200451-1-0
    chr2    157832187       157832245       COMPLEX_8196-2-0
```
```
    $ head CHDP.rgn
    --- For data requires SubGEM, the complexes are names as "COMPLEX_ID-nth SubGEM-Total fragment in the SubGEM-nth fragment" ---
    chr2L   12713870        12714100        COMPLEX_2049-0-1-0
    chr2R   7863021 7863695 COMPLEX_2049-1-1-0
    chr3L   11185707        11185889        COMPLEX_2049-2-1-0
    chr3R   3018595 3018777 COMPLEX_2049-3-1-0
    chr2R   6408379 6408606 COMPLEX_4098-0-2-0
    chr2R   12027465        12027974        COMPLEX_4098-0-2-1
```
#### raw_matrix
There were at least three files under the directory. Showing the result for all the barcodes. The name of the directory do not necessarily named as raw_matrix, it can be raw for scRNA-seq related pipeline under Solo.out/Gene/.
|File name|Description|Experiment|
|---|---|---|
|features.tsv|Gene information as row name of expression matrix|scRNA-seq, scRNA-seq in scARC-seq|
|barcodes.tsv|Barcode information as column name of expression matrix|scRNA-seq, scRNA-seq in scARC-seq|
|matrix.mtx|Expression matrix/peak cell fragment count matrix as matrix market format|scRNA-seq, scATAC-seq, scRNA-seq in scARC-seq, scATAC-seq in scARC-seq|
|barcodes.tsv_BARPID.txt|Converted barcode sequence to BARP identifier ID|scRNA-seq in scARC-seq|
|BarcodeToBARPIDConvertTable.tsv|Pairwise conversion from barcode to BARP ID| scRNA-seq in scARC-seq|
|Peak|Peak information as row name of peak cell fragment count matrix|scATAC-seq, scATAC-seq in scARC-seq|
|Cell|Cell information as column name of peak cell fragment count matrix| scATAC-seq, scATAC-seq in scARC-seq|
|filtered_cells.tsv| Filtered barcodes with emptyDrops|scATAC-seq, scATAC-seq in scARC-seq|

```
    *** $ head features.tsv ***
    --- Gene ID - Gene symbol - Gene Expression ---
    ENSG00000223972.5       DDX11L1 Gene Expression
    ENSG00000227232.5       WASH7P  Gene Expression
    ENSG00000278267.1       MIR6859-1       Gene Expression
    ENSG00000243485.5       MIR1302-2HG     Gene Expression
    ENSG00000284332.1       MIR1302-2       Gene Expression
```
```
    *** $ head barcodes.tsv ***
    AAACAGCCAAACCTAT
    AAACAGCCAAACCTTG
    AAACAGCCAAACGGGC
    AAACAGCCAAACTAAG
    AAACAGCCAAACTGCC
```
```
    *** $ head matrix.mtx ***
    --- Header line (necessary) tells the attributes of MatrixMarket ---
    --- %%MatrixMarket - object - format - field - symmetry ---
    %%MatrixMarket matrix coordinate integer general
    --- Annotation line ---
    %
    --- row number - column number - total value ---
    58780 303576 3981045
    --- row index column index - value ---
    665 1 1
    2910 1 1
    3236 1 1
    6355 1 1
    6475 1 1
    17632 1 1
```
Detailed information about MatrixMarket([English](https://math.nist.gov/MatrixMarket/formats.html), [中文](https://zhuanlan.zhihu.com/p/147640256)).
```
    $ head barcodes.tsv_BARPID.txt
    CELL_16833
    CELL_60050
    CELL_149669
    CELL_174442
    CELL_120091
    CELL_119651
```

```
    $ head BarcodeToBARPIDConvertTable.tsv
    AAACAGCCAAACCTAT        CELL_16833
    AAACAGCCAAACCTTG        CELL_60050
    AAACAGCCAAACGGGC        CELL_149669
    AAACAGCCAAACTAAG        CELL_174442
    AAACAGCCAAACTGCC        CELL_120091
```
```
    head filtered_cells.tsv
    CELL_2585
    CELL_53
    CELL_1174
    CELL_2272
    CELL_664
```
#### filtered_matrix
There were at least three files under the directory. Showing the result filtered for real cells. The name of the directory do not necessarily named as raw_matrix, it can be raw for scRNA-seq related pipeline under Solo.out/Gene/.
|File name|Description|Experiment|
|---|---|---|
|features.tsv|Gene information as row name of expression matrix|scRNA-seq, scRNA-seq in scARC-seq|
|barcodes.tsv|Barcode information as column name of expression matrix|scRNA-seq, scRNA-seq in scARC-seq|
|matrix.mtx|Expression matrix/peak cell fragment count matrix as matrix market format|scRNA-seq, scATAC-seq, scRNA-seq in scARC-seq, scATAC-seq in scARC-seq|
|barcodes.tsv_BARPID.txt|Converted barcode sequence to BARP identifier ID|scRNA-seq in scARC-seq|
|BarcodeToBARPIDConvertTable.tsv|Pairwise conversion from barcode to BARP ID| scRNA-seq in scARC-seq|
|Peak|Peak information as row name of peak cell fragment count matrix|scATAC-seq, scATAC-seq in scARC-seq|
|Cell|Cell information as column name of peak cell fragment count matrix| scATAC-seq, scATAC-seq in scARC-seq|
|results.tsv|Pairwise cell to cell same transposition event count matrix [scATAC-seq algorithm](#scatac-seq-algorithm)| scATAC-seq, scATAC-seq in scARC-seq|

|-|CELL_10|CELL_11|CELL_524|CELL_13|CELL_39|CELL_101|CELL_105|CELL_106|CELL_108|CELL_958|
|---|---|---|---|---|---|---|---|---|---|---|
|**CELL_10**|31758|200|142|331|118|294|284|344|387|32|
|**CELL_11**|200|2540|45|93|23|77|110|98|91|6|
|**CELL_524**|142|45|2656|68|15|64|75|96|64|527|
|**CELL_13**|331|93|68|5106|65|100|177|203|191|11|
|**CELL_39**|118|23|15|65|786|41|54|50|69|3|
|**CELL_101**|294|77|64|100|41|3214|148|175|128|18|
|**CELL_105**|284|110|75|177|54|148|5204|225|148|13
|**CELL_106**|344|98|96|203|50|175|225|4784|159|24|
|**CELL_108**|387|91|64|191|69|128|148|159|6388|17|
|**CELL_958**|32|6|527|11|3|18|13|24|17|78|

## Algorithm
### Barcode identification
Please refer to the paper for detailed barcode identification algorithm.
### Identifier 
Please refer to the supplementary note for detailed identifier algorithm.
### scATAC-seq algorithm
Please refer to the supplementary note for detailed scATAC-seq algorithm.
