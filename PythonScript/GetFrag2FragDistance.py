# This scripts only applies for BARP generated COMPLEX file, because it is sorted.

import sys

FI = open(sys.argv[1], "r")
# PairLoopOUT = open("L2L.stat", "w")
FLOUT = open(sys.argv[1] + ".FragLen.stat", "w")
F2FOUT = open(sys.argv[1] + ".F2F.stat", "w")
CLOUT = open(sys.argv[1] + ".ComplexSpan.stat", "w")

# JuiceFile = open(sys.argv[1] + ".juice", 'w')

CountColumn = int(sys.argv[2])


class Fragment:
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end

    def getFragLen(self):
        return int(self.end) - int(self.start)


Count = 0
FragCountList = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
while True:
    line = FI.readline()
    if not line:
        break
    if line.startswith("@") or line.startswith("#"):
        continue
    word = line.strip().split("\t")
    Complex = ""
    for i in range(CountColumn - 1):
        Complex = Complex + word[i] + "\t"
    Complex = Complex[0:-1]
    FragCount = int(word[CountColumn - 1])
    if FragCount > 15:
        FragCountList[15] = FragCountList[15] + 1
    else:
        FragCountList[FragCount - 1] = FragCountList[FragCount - 1] + 1
    CoordList = []
    i = 0
    while i < FragCount:
        FragMent = Fragment(word[i * 3 + CountColumn], word[i * 3 + CountColumn + 1], word[i * 3 + CountColumn + 2])
        FragMentLen = FragMent.getFragLen()
        CoordList.append(int(word[i * 3 + CountColumn + 1]))
        CoordList.append(int(word[i * 3 + CountColumn + 2]))
        FLenLine = FragMent.chr + "\t" + FragMent.start + "\t" + FragMent.end + "\t" + Complex + "-" + str(
            FragCount) + "-" + str(i) + "\t" + str(FragMentLen) + "\n"
        FLOUT.writelines(FLenLine)
        F2FFlag = 1
        j = i + 1
        # while j < int((len(word) - 2) / 3):
        #     NeighBorFrag = Fragment(word[j * 3 + 2], word[j * 3 + 3], word[j * 3 + 4])
        #     PairLoop = int(NeighBorFrag.start) - int(FragMent.end)
        #     FDistLine = FragMent.chr + "\t" + FragMent.start + "\t" + FragMent.end + "\t" + NeighBorFrag.chr + "\t" + NeighBorFrag.start + "\t" + NeighBorFrag.end + "\t" + Complex + "-" + str(
        #         i) + "-" + str(j) + "\t" + str(PairLoop) + "\t" + str(FragCount) + "\n"
        #     JuiceLine = "0" + " " + FragMent.chr.strip("chr") + " " + str(int((int(FragMent.start) + int(FragMent.end))/2)) + " " + "0" + " " + \
        #                 "0" + " " + NeighBorFrag.chr.strip("chr") + " " + str(int((int(NeighBorFrag.start) + int(NeighBorFrag.end))/2)) + " " + "1" + "\n"
        #     JuiceFile.writelines(JuiceLine)
        #     if F2FFlag == 1:
        #         F2FOUT.writelines(FDistLine)
        #         F2FFlag = 0
        #     PairLoopOUT.writelines(FDistLine)
        #     j = j + 1
        i = i + 1
    FragLen = max(CoordList) - min(CoordList)
    CLLine = Complex + "\t" + str(FragLen) + "\t" + str(FragCount) + "\n"
    CLOUT.writelines(CLLine)

    Count = Count + 1
    if Count % 10000 == 0:
        sys.stderr.write("[processed] " + str(Count) + "\n")

with open(sys.argv[1] + ".FragCount.stat", "w") as file:
    for i in range(len(FragCountList)):
        if i < 15:
            file.write("F = " + str(i + 1) + "\t" + str(FragCountList[i]) + "\n")
        else:
            file.write("F > 15" + "\t" + str(FragCountList[i]) + "\n")
