import os, sys
import getopt
import time


class ArrstrWithFlag:
    StrArr = []
    FlagArr = []
    LoopFlag = 0

def Usage():
    print("Usage: barp rdp [options...] [in.cluster]\n")
    print("Options:\n")
    print("  -i FILE    Inline bed file.\n")
    print("  -r STR     Remove duplicates by field.\n")
    print("  -f STR     Re-specify the representation of field.\n")
    print("  -o FILE    Output prefix.\n")
    print("  -d CHAR    Split symbol used in inline bed file.\n")

def parseSparseStr(SparseStr):
#     flag shows what this position represent: 0:2 field and order of barcode; 1:2 field and name of barcode; 2:3 field
#     and order of barcode;3:3 field and name of barcode; 4: tagname
    Ret = ArrstrWithFlag()
    AddFieldArr = SparseStr.split("-")
    for i in range(len(AddFieldArr)):
        CurADF = AddFieldArr[i]
        PosArr = CurADF.split(":")
        if len(PosArr) == 2:
            if len(PosArr[1]) != 2:
                print("Field specified error, should be like: 3:UMI:UM-DPM:DP-5:CP")
                exit(1)
            else:
                if PosArr[0].isdigit():
                    Ret.FlagArr.append(0)
                    Ret.FlagArr.append(4)
                else:
                    Ret.FlagArr.append(1)
                    Ret.FlagArr.append(4)
                    Ret.LoopFlag = 1
                Ret.StrArr.append(PosArr[0])
                Ret.StrArr.append(PosArr[1])
        elif len(PosArr) == 3:
            if len(PosArr[2]) != 2:
                print("Field specified error, should be like: 3:UMI:UM-DPM:DP-5:CP")
                exit(1)
            if PosArr[0].isdigit():
                Ret.FlagArr.append(2)
                Ret.FlagArr.append(3)
                Ret.FlagArr.append(4)
                Ret.StrArr.append(PosArr[0])
                Ret.StrArr.append(PosArr[1])
                Ret.StrArr.append(PosArr[2])
        else:
            print("Field specified error, should be like: 3:UMI:UM-DPM:DP-5:CP")
            exit(1)
    return Ret


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:s:o:", ["input=", "sparse-str", "output="])
    except getopt.GetoptError:
        Usage()
        sys.exit(2)

    InputFile = ''
    OutputPrefix = time.strftime("%b%d", time.localtime()).upper()
    UserPrefix = ''
    FieldName = ''
    RemoveStr = ''
    FileSplit = '\t'
    CommentSplit = '\t'
    CommandLineFlag = 1

    for opt, arg in opts:
        if opt == 'h':
            Usage()
        elif opt in ("-i", "--input"):
            InputFile = arg
        elif opt in ("-r", "--sparse-str"):
            SparseStr = arg
        elif opt in "-o":
            UserPrefix = arg
        else:
            print("Unrecognized parameter " + opt + ", ignoring...\n")

    FIN = open(InputFile, "r", encoding='utf-8')

    AddFieldStr = parseSparseStr(SparseStr)
    LineCount = 0
    while True:
        line = FIN.readline().strip()
        LineCount = LineCount + 1
        if not line:
            break
        if LineCount%4 == 1:
            FieldToBam = []
            BarChain = line.split("|||")[1].split("|")
            if AddFieldStr.LoopFlag == 1:
                for i in range(len(BarChain)):
                    TagContent = BarChain[i].split(":")
                    CodeName = TagContent[0]
                    CodeContent = TagContent[1]
                    for j in range(len(AddFieldStr.StrArr)):
                        if AddFieldStr.FlagArr[j] == 1 or AddFieldStr.FlagArr[j] == 3:
                            if AddFieldStr.StrArr[j] == CodeName:
                                FieldToBam.append(AddFieldStr.StrArr[j+1][:2] + ":Z:" + CodeContent)
                                j = j + 1





