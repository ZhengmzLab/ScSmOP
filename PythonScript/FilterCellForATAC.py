import sys

FI = open(sys.argv[1], "r")
FFilter = open(sys.argv[2], "r")
FPeak = open("Peak", "w")
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
Result = {}


MatrixMarketComment = "%%MatrixMarket matrix coordinate integer general\n"
FMatrix.write(MatrixMarketComment)

PeakCount = 1
CellCount = 1
LineCount = 1
MatrixList = []

# read filtered cell
FilteredCellDict = {}
while True:
    line = FFilter.readline().strip()
    if not line:
        break
    line = line.split()
    CellName = line[0]
    FilteredCellDict.update({CellName: CellCount})
    CellCount = CellCount + 1

# filter peak
while True:
    line = FI.readline().strip()
    if not line:
        break
    line = line.split()
    CurPeak = line[0] + "-" + line[1] + "-" + line[2]
    CurCell = line[3]

    if not FilteredCellDict.get(CurCell):
        continue
    else:
        CellIndex = FilteredCellDict.get(CurCell)

        if not PeakList.get(CurPeak):
            PeakList.update({CurPeak: PeakCount})
            PeakLine = line[0] + "\t" + line[1] + "\t" + line[2] + "\n"
            FPeak.write(PeakLine)
            PeakIndex = PeakCount
            PeakCount = PeakCount + 1
        else:
            PeakIndex = PeakList.get(CurPeak)
        FragCount = line[4]
        MatrixLine = str(PeakIndex) + "\t" + str(CellIndex) + "\t" + FragCount + "\n"
        MatrixList.append(MatrixLine)
        LineCount = LineCount + 1


DimensionLine = str(PeakCount - 1) + "\t" + str(CellCount - 1) + "\t" + str(LineCount - 1) + "\n"
FMatrix.write(DimensionLine)

for PeakCellCount in MatrixList:
    FMatrix.write(PeakCellCount)

FPeak.close()
FMatrix.close()
FFilter.close()
FI.close()
print("Generating Matrix Finished. Processed: " + str(LineCount) + ".")