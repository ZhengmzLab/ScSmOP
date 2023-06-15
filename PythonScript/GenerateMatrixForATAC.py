import sys

FI = open(sys.argv[1], "r")
FPeak = open("Peak", "w")
FCell = open("Cell", "w")
FMatrix = open("matrix.mtx", "w")


class Fragment:
    def __init__(self, chr, start, end):
        self.chr = chr
        self.start = start
        self.end = end


class CellCount:
    def __init__(self, name, count):
        self.name = name
        self.count = count


PeakList = {}
CellList = {}
Result = {}

MatrixMarketComment = "%%MatrixMarket matrix coordinate integer general\n"
FMatrix.write(MatrixMarketComment)

PeakCount = 1
CellCount = 1
LineCount = 1
MatrixList = []


while True:
    line = FI.readline().strip()
    if not line:
        break
    line = line.split()
    CurPeak = line[0] + "-" + line[1] + "-" + line[2]
    CurCell = line[3]

    LineCount = LineCount + 1
    if not PeakList.get(CurPeak):
        PeakList.update({CurPeak: PeakCount})
        PeakLine = line[0] + "\t" + line[1] + "\t" + line[2] + "\n"
        FPeak.write(PeakLine)
        PeakIndex = PeakCount
        PeakCount = PeakCount + 1
    else:
        PeakIndex = PeakList.get(CurPeak)
    if not CellList.get(CurCell):
        CellList.update({CurCell: CellCount})
        CellLine = CurCell + "\n"
        FCell.write(CellLine)
        CellIndex = CellCount
        CellCount = CellCount + 1
    else:
        CellIndex = CellList.get(CurCell)

    FragCount = line[4]
    MatrixLine = str(PeakIndex) + "\t" + str(CellIndex) + "\t" + FragCount + "\n"
    MatrixList.append(MatrixLine)

    if LineCount % 10000 == 0:
        sys.stderr.write("[processed] " + str(LineCount) + "\n")

DimensionLine = str(PeakCount - 1) + "\t" + str(CellCount - 1) + "\t" + str(LineCount - 1) + "\n"
FMatrix.write(DimensionLine)

for PeakCellCount in MatrixList:
    FMatrix.write(PeakCellCount)

FCell.close()
FPeak.close()
FMatrix.close()
FI.close()
print("Generate Matrix Finished. Processed: " + str(LineCount) + ".")
