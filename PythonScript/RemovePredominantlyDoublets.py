import sys
from itertools import combinations_with_replacement
import numpy as np
import pandas as pd

FIN = open(sys.argv[1], "r")
FOUT = open("Cell", "w")
SelectedCell = open(sys.argv[2], "r")


class TranspositionEvent:
    def __init__(self):
        self.ID = 0
        self.FragCount = 0
        self.Start = set()
        self.End = set()
        self.Coord = ""


def compareTranspositionEvent(A):
    # print(list(A[0].Start.intersection(A[1].End)), list(set(A[1].Start).intersection(A[0].End)))
    return len(A[0].Start.intersection(A[1].End)) + len(A[1].Start.intersection(A[0].End))


BarcodeFragmentList = []
fsplit = "\t"
skip = 1
LineCount = 0
IDList = []
GenerateWholeTable = False
print("Reading fragments...")
SelectedCellSet = set()


#
while True:
    line = SelectedCell.readline().strip()
    if not line:
        break
    SelectedCellSet.add(line)

FragCountList = []
while True:
    line = FIN.readline().strip()
    if not line:
        break
    if line[0] == '#' or line[0] == '@':
        continue
    else:
        InfoArray = line.split(fsplit)
        Cell = InfoArray[0]
        if Cell in SelectedCellSet:
            if len(InfoArray) > skip + 1:
                FragCount = int(InfoArray[skip])
                BarcodeEvent = TranspositionEvent()
                BarcodeEvent.FragCount = FragCount
                FragCountList.append(FragCount)
                # BarcodeEvent.Coord = line
                BarcodeEvent.ID = LineCount
                IDList.append("|".join(InfoArray[0:skip]))
                for i in range(FragCount):
                    BarcodeEvent.Start.add(InfoArray[i * 3 + skip + 1] + "-" + InfoArray[i * 3 + skip + 2])
                    BarcodeEvent.End.add(InfoArray[i * 3 + skip + 1] + "-" + InfoArray[i * 3 + skip + 3])
                BarcodeFragmentList.append(BarcodeEvent)
            LineCount = LineCount + 1

print("Calculating similarity")
# if GenerateWholeTable:
ComparePairList = combinations_with_replacement(BarcodeFragmentList, 2)
CompareResult = pd.DataFrame(columns=range(len(BarcodeFragmentList)), index=range(len(BarcodeFragmentList)))
#
Count = 0
CellCount = 0
for ComparePair in ComparePairList:
    Count = Count + 1
    if Count % 768 == 0:
        CellCount = CellCount + 1
        sys.stderr.write("[processed] " + str(CellCount) + "\n")
    # print(ComparePair[0].ID, ComparePair[1].ID)
    EventCount = compareTranspositionEvent(ComparePair)
    # print(EventCount)
    # print(ComparePair[0].ID, ComparePair[1].ID)
    CompareResult.iloc[ComparePair[0].ID, ComparePair[1].ID] = EventCount
    CompareResult.iloc[ComparePair[1].ID, ComparePair[0].ID] = EventCount

CompareResult.columns = IDList
CompareResult.index = IDList
CompareResult = CompareResult.astype(int)
for index, row in CompareResult.iterrows():
    if index == row.idxmax():
        FOUT.write(index + "\n")
CompareResult.to_csv('result.csv')


