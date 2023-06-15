import pandas as pd
import os, sys
import re




def main():

    GenomeType = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                  'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

    LibInfo = "rdSPRITE.xlsx"
    ReadStructureInfo = pd.DataFrame
    if os.path.exists(LibInfo):
        FileType = os.path.splitext(LibInfo)[-1]
        if FileType in ["", ".txt", ".tsv"]:
            ReadStructureInfo = pd.read_csv(LibInfo, sep="\t", header=None, encoding='utf-8', comment="#")
        elif FileType in [".xlsx", ".xls"]:
            ReadStructureInfo = pd.read_excel(LibInfo, header=None, comment="#")
        elif FileType in [".csv"]:
            ReadStructureInfo = pd.read_csv(LibInfo, header=None, encoding='utf-8', comment="#")
        else:
            sys.stderr.write("Unknow type of file, exiting...\n")
    else:
        sys.stderr.write("Config file do not exist, exiting...\n")

    ReadStructureInfo = ReadStructureInfo.dropna(axis=0, how="all")

    print(ReadStructureInfo)

    GenomeStrList = []
    ReadStructurePattern = re.compile(r'(readstructure\s*)(\S*)', re.I)
    for index, row in ReadStructureInfo.iterrows():
        if len(row) < 2:
            continue
        Option = "".join(row[0].lower().split())

        # if the option read structure 
        ReadStructureMatch = ReadStructurePattern.search(Option)
        StructureExist = False
        if ReadStructureMatch:
            GenomeStr = str(ReadStructureMatch.group(2))
            for i in range(len(GenomeStrList)):
                if  GenomeStr == GenomeStrList[i]:
                    CurrentGenome = GenomeType[i]
                    StructureExist = True
                    break
            if not StructureExist:
                GenomeStrList.append(GenomeStr)
                CurrentGenome = GenomeType[len(GenomeStrList) - 1]
            
            # read structure to barcode chain:
            BarcodeStrList = row[1].split('-')
            if len(BarcodeStrList) < 2:
                sys.stderr.write("Read structure shoud have at least 1 barcode or GENOME")
                exit(1)
            # barcode chain read position
            ReadPosPattern = re.compile(r'^R(ead|d)?\s*([1234]{1})', re.I)
            ReadPosMatch = ReadPosPattern.search(BarcodeStrList[0])
            if ReadPosMatch:
                ReadPosStr = "R" + str(ReadPosMatch.group(2)) + ":1"
            else:
                sys.stderr.write("Read structure shoud start with Read 1 or Read 2")
                exit(1)
            
            # parse barcode chain
            BracketPattern = re.compile(r'(\(\))', re.I)
            NucleotidePattern = re.compile(r'([RYMKSWBVDHN])', re.I)
            # space accumulate till a barcode find
            Space = 0
            for i in range(1, len(BarcodeStrList)):
                BracketMatch = BracketPattern.search(BarcodeStrList[i])
                print(BarcodeStrList[i])
                # if the item is barcode or space
                if not BracketMatch:
                    # nucleotide sequence 
                    NucleotideMatch = NucleotidePattern.search(BarcodeStrList[i])
                    if NucleotideMatch:
                        print("Nucleotide: " + str(NucleotideMatch.group(1))) 


                # # if the barcode is ligation based:
                # BarcodeStruc = BarcodeStrList[i].split('~')
                # if len(BarcodeStruc > 1):
                #     # barcode is ligation based, the first part is space.
                #     Space = BarcodeStruc[0].strip("bp")
                #     Laxity = 6
                #     BarcodeBracket = BarcodeStruc[1]
                # else:
                #     # barcode is not ligation based.
                #     BarcodeStruc = 
                # print(BarcodeStrList[i])
            


                
                
    
    print(GenomeStrList)






if __name__ == "__main__":
    main()


 