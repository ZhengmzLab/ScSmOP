import os, sys
import getopt
import time
from typing import List, Any


class Fragment:
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end
        self.target = ""
        self.adds = []
        self.dup = 1

    def getFragLen(self):
        return int(self.end) - int(self.start)


class RPST:
    def __init__(self):
        # field record name
        self.Frn = 0
        self.Fathers = []
        self.AddAr = []
        self.frags = []


def Usage():
    print("Usage: barp rdp [options...] [in.cluster]\n")
    print("Options:\n")
    print("  -i FILE    Inline bed file.\n")
    print("  -r STR     Remove duplicates by field.\n")
    print("  -f STR     Re-specify the representation of field.\n")
    print("  -o FILE    Output prefix.\n")
    print("  -d CHAR    Split symbol used in inline bed file.\n")


def readHeader(fi):
    with open(fi, "r", encoding="utf-8") as file:
        while True:
            line = file.readline().strip()
            if not line:
                break
            if not (line[0] == '#' or line[0] == '@'):
                return "-"
            if line.startswith("@CN"):
                return line


def parseField(fdn, split, CommandLineFlag):
    rpst = RPST()
    Tmp = fdn.split(split)
    CountFlag = 1
    for i in range(CommandLineFlag, len(Tmp)):
        if CountFlag:
            if Tmp[i].lower() == "count":
                CountFlag = 0
                continue
            rpst.Fathers.append(Tmp[i].lower())
        else:
            rpst.AddAr.append(Tmp[i].lower())
    return rpst


class ColParam:
    def __init__(self):
        self.TargetCol = []
        self.RestCol = []
        self.Skip = int


def parseCol(rmstr, rpst):
    Final = ColParam()
    FieldBefore = len(rpst.Fathers)
    Tmp = rmstr.split('-')

    for i in range(len(Tmp)):
        for j in range(len(rpst.AddAr)):
            if Tmp[i].lower() == rpst.AddAr[i]:
                Final.TargetCol.append(j + 3)
        if Tmp[i].lower() == "chr":
            Final.TargetCol.append(0)
        if Tmp[i].lower() == "start":
            Final.TargetCol.append(0)
            Final.TargetCol.append(1)
        elif Tmp[i].lower() == "end":
            Final.TargetCol.append(0)
            Final.TargetCol.append(2)
    Final.TargetCol = list(set(Final.TargetCol))
    Final.TargetCol.sort()
    for i in range(len(rpst.AddAr)):
        if rpst.AddAr[i] not in Final.TargetCol:
            Final.RestCol.append(i)
    Final.RestCol = list(set(Final.RestCol))
    Final.RestCol.sort()
    Final.Skip = FieldBefore

    return Final


def removeCore(fi, fo, col, fsplit, dupfile):
    fn = open(fi, "r", encoding='utf-8')
    fp = open(fo, "a", encoding='utf-8')

    skip = col.Skip
    LineCount = 0
    Duplicates = 0
    TotalReads = 0

    while True:
        line = fn.readline().strip()
        if not line:
            break
        if line[0] == '#' or line[0] == '@':
            continue
        else:
            LineCount = LineCount + 1
            InfoArray = line.split(fsplit)
            NewCount = 0
            QuaFragCount = 0
            if len(InfoArray) > skip + 1:
                for i in range(skip):
                    fp.write(InfoArray[i] + "\t")
                FragCount = int(InfoArray[skip])
                TotalReads = TotalReads + FragCount
                FragCompList = {}
                FragSetComp = set()
                for i in range(FragCount):
                    Frag = Fragment(InfoArray[i * 3 + skip + 1], InfoArray[i * 3 + skip + 2],
                                    InfoArray[i * 3 + skip + 3])
                    for j in range(len(col.TargetCol)):
                        Field = col.TargetCol[j]
                        if Field < 3:
                            Frag.target = Frag.target + "-" + InfoArray[i * 3 + skip + 1 + Field]
                        else:
                            Frag.target = Frag.target + "-" + InfoArray + InfoArray[3 * FragCount + skip + 1 + Field - 3]
                    for j in range(len(col.RestCol)):
                        Field = col.TargetCol[j]
                        Frag.adds.append(InfoArray[3 * FragCount + skip + 1 + Field - 3])
                    if Frag.target not in FragSetComp:
                        FragSetComp.add(Frag.target)
                        FragCompList.update({Frag.target: Frag})
                    else:
                        Duplicates = Duplicates + 1
                        FragCompList[Frag.target].dup = FragCompList[Frag.target].dup + 1
                        # print(FragCompList.index(Frag.target))
                        # FragList[FragCompList.index(Frag.target)].dup = 2
                        continue
                fp.write(str(len(FragCompList)) + "\t")
                for ReFrag in FragCompList:
                    Frag = FragCompList[ReFrag]
                    # print("One" + ReFrag.chr + "\t" + ReFrag.start + "\t" + ReFrag.end + "\t")
                    fp.write(Frag.chr + "\t" + Frag.start + "\t" + Frag.end + "\t")
                for ReFrag in FragCompList:
                    Frag = FragCompList[ReFrag]
                    fp.write(str(Frag.dup) + "\t")
                # for ReFrag in FragSet:
                #     print(ReFrag.target)
                #     for adds in ReFrag.target:
                #         fp.write(adds + "\t")
                #     for adds in ReFrag.adds:
                #         print(adds + "\t")
                #         fp.write(adds + "\t")
                fp.write("\n")
        LineCount = LineCount + 1
        if LineCount % 10000 == 0:
            sys.stderr.write("[processed] " + str(LineCount) + "\n")
    with open(dupfile, 'w') as F:
        F.write("Duplicates:" + str(Duplicates) + "\tTotalFragments:" + str(TotalReads) + "\n")


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:r:f:o:d:", ["input=", "remove-by=", "field-order=", "output="])
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
        elif opt in ("-r", "--remove-by"):
            RemoveStr = arg
        elif opt in ("-f", "--field-order"):
            FieldName = arg
        elif opt in "-o":
            UserPrefix = arg
        elif opt in "-d":
            FileSplit = arg
        else:
            print("Unrecognized parameter " + opt + ", ignoring...\n")

    if FieldName == '':
        FieldName = readHeader(InputFile)
    else:
        CommandLineFlag = 0
        CommentSplit = '-'

    rpst = parseField(FieldName, CommentSplit, CommandLineFlag)
    ColumnArray = parseCol(RemoveStr, rpst)

    if UserPrefix == '':
        FOUT = OutputPrefix + "RDP.cluster"
        DupFOUT = OutputPrefix + ".Dup.stat"
    else:
        FOUT = UserPrefix + ".RDP.cluster"
        DupFOUT = UserPrefix + ".Dup.stat"

    with open(FOUT, "w", encoding='utf-8') as file:
        file.write("# This file is produced by barp rdp\n@CN")
        for father in rpst.Fathers:
            file.write("\t" + father)
        file.write("\tCOUNT")
        for adds in ColumnArray.TargetCol:
            if adds > 2:
                file.write("\t" + rpst.AddAr[adds - 3])
        for adds in ColumnArray.RestCol:
            file.write("\t" + rpst.AddAr[adds - 3])
        file.write("\n")

    if len(ColumnArray.TargetCol) == 0:
        print("No column specified to remove duplicates\n")
        exit(2)
    removeCore(InputFile, FOUT, ColumnArray, FileSplit, DupFOUT)


if __name__ == '__main__':
    main()
