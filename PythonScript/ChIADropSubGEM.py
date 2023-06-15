import os,sys

IN = sys.argv[1]
FIN = open(IN, "r")
FOUT = open(IN+".SUBGEM", "w")
CountColumn = int(sys.argv[2])


class Fragment():
    def __init__(self, chr, start, end):
        self.chr = chr
        self.s = start
        self.e = end
        self.add = []

while True:
    line = FIN.readline()
    if not line:
        break
    if line.startswith("#") or line.startswith("@"):
        continue
    ComplexInfo = line.strip().split("\t")
    FragInfo = ComplexInfo[0:CountColumn]
    FragCount = ComplexInfo[CountColumn - 1]
    Frags = ComplexInfo[CountColumn:]
    SubGEMList = []
    CurGEM = []
    CurChrList = []
    for i in range(int(FragCount)):
        CurFrag = Fragment(Frags[i * 3], Frags[i * 3 + 1], Frags[i * 3 + 2])
        Pos = CurChrList.index(CurFrag.chr) if CurFrag.chr in CurChrList else -1
        if Pos < 0:
            CurChrList.append(CurFrag.chr)
            SubGEMList.append([CurFrag])
        else:
            SubGEMList[Pos].append(CurFrag)
    for i in range(len(SubGEMList)):
        SubGEM = SubGEMList[i]
        LineWrite = ""
        for j in range(len(FragInfo) - 2):
            LineWrite = LineWrite + FragInfo[j] + "\t"
        LineWrite = LineWrite + FragInfo[CountColumn - 2] + "-" + str(i) + "\t" + str(len(SubGEM))
        FOUT.write(LineWrite)
        for frag in SubGEM:
            FOUT.write("\t" + frag.chr + "\t" + frag.s + "\t" + frag.e)
        FOUT.write("\n")
