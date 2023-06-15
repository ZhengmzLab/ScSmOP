import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate
import math
import sys
from itertools import combinations_with_replacement
import getopt
import time


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


# filter cells with the number of fragments/ adapted from dropletutils barcoderanks
# m: a numeric matrix-like object containing UMI counts, column represent droplet and rows represent genes/peaks; dataframe, matrix
# lower: a numeric scalar specifying
# fit bounds: a numeric vector of length 2, specifying the lower and upper bounds on total UMI count from which to
# obtain a of the curve for spline fitting.
# exclude from: an integer scalar specifying the number of highest ranking barcodes to exclude from spline fitting
# df : integer scalar specifying the number of degrees of freedom
def findCurveBound(x, y, ExcludeFrom):
    #     remove the first ExcludeFrom beads for finding of inflection
    DLN = np.diff(y) / np.diff(x)
    # print(DLN)
    if ExcludeFrom != 0:
        Skip = min(len(DLN) - 1, sum(x <= math.log10(ExcludeFrom)))
    else:
        Skip = 0

    DLN = DLN[Skip:]
    RightEdge = np.argmin(DLN)
    LeftEdge = np.argmax(DLN)

    return LeftEdge + Skip, RightEdge + Skip


def barcodeRanks(M, Lower=100, FitBounds=None, ExcludeFrom=50):
    if isinstance(M, pd.DataFrame):
        BarcodeFragCounts = M.apply(lambda x: x.sum())
    elif isinstance(M, np.ndarray):
        if M.ndim == 1:
            BarcodeFragCounts = M
        elif M.ndim == 2:
            BarcodeFragCounts = np.sum(M, axis=0)
        else:
            raise Exception("Matrix dimension is not 1 or 2")
    elif isinstance(M, list):
        M = np.array(M)
        if M.ndim == 1:
            BarcodeFragCounts = M
        elif M.ndim == 2:
            BarcodeFragCounts = np.sum(M, axis=0)
        else:
            raise Exception("Matrix dimension is not 1 or 2")
    else:
        raise Exception("Matrix is not nd array or data frame")

    SortedBarcodeCount = abs(np.sort(-BarcodeFragCounts))
    RankArray = np.arange(1, len(SortedBarcodeCount) + 1, 1)
    DuplicatedBarcodeCount = np.array(np.unique(SortedBarcodeCount, return_counts=True))
    BarcodeRankFlip = np.flip(DuplicatedBarcodeCount[1])
    #     take the middle position of each same counted barcode as the rank
    BarcodeRank = np.cumsum(BarcodeRankFlip) - (BarcodeRankFlip - 1) / 2
    BarcodeCount = np.flip(DuplicatedBarcodeCount[0])

    # print("BarcodeRank", BarcodeRank)
    # print("BarcodeCount", BarcodeCount)
    Keep = BarcodeCount > Lower
    if sum(Keep) < 3:
        raise Exception("Insufficient unique position for computing knee/inflection points")

    x = np.log10(BarcodeRank[Keep])
    y = np.log10(BarcodeCount[Keep])

    EdgeOut = findCurveBound(x, y, ExcludeFrom)
    LeftEdge = EdgeOut[0]
    RightEdge = EdgeOut[1]

    Inflection = pow(10, y[RightEdge])
    RealRight = SortedBarcodeCount[SortedBarcodeCount >= Inflection]
    if FitBounds is None:
        NewKeep = np.arange(LeftEdge, RightEdge, 1)
    else:
        NewKeep = math.log10(FitBounds[0]) < y < math.log10(FitBounds[1])

    if RightEdge - LeftEdge >= 4:
        fit = interpolate.splrep(x[NewKeep], y[NewKeep], s=0)
        # FitValues = pow(10, fit(x[NewKeep]))
        Dev1 = interpolate.splev(x[NewKeep], fit, der=1)
        Dev2 = interpolate.splev(x[NewKeep], fit, der=2)
        Curvature = Dev2 / (pow(1 + pow(Dev1, 2), 1.5))
        # print("Curvature", Curvature)
        Knee = pow(10, y[np.argmin(Curvature)])
        KneeRank = np.argmin(Curvature)
        PlotKnee = y[np.argmin(Curvature)]
    else:
        Knee = pow(10, y[NewKeep[0]])
        PlotKnee = y[NewKeep[0]]

    plt.plot(RankArray[:len(RealRight)], SortedBarcodeCount[:len(RealRight)], color='orange')
    plt.plot(RankArray[len(RealRight):], SortedBarcodeCount[len(RealRight):],  color='blue')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("CellFilter.png")
    return Inflection


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:r:o:c:w", ["input=", "raw=", "output=", "countcol=", "wholetable="])
    except getopt.GetoptError:
        print(sys.argv[0] + "-i InputConfigFile -o OutputConfigFile -r RawBarcode -c CountColumn -w")
        sys.exit(2)

    InputFile = ''
    OutputPrefix = time.strftime("%b%d", time.localtime()).upper()
    UserPrefix = ''
    skip = 1
    RawBarcodeFile = ''
    GenerateWholeTable = False

    for opt, arg in opts:
        if opt == 'h':
            print(sys.argv[0] + "-i InputConfigFile -o OutputConfigFile -r RawBarcode -c CountColumn")
            sys.exit()
        elif opt in ("-i", "--input"):
            InputFile = arg
        elif opt in ("-o", "--output"):
            UserPrefix = arg
        elif opt in ("-c", "--countcol"):
            skip = int(arg)
        elif opt in ("-r", "--raw"):
            RawBarcodeFile = arg
        elif opt in ("-w", "--wholetable"):
            GenerateWholeTable = True
        else:
            print("Unrecognized parameter " + opt + ", ignoring...\n")

    if InputFile == '':
        raise Exception("No input cluster provided.\n")

    if UserPrefix != '':
        OutputPrefix = UserPrefix
    RawBarcodeFlag = False
    RawCellSet = set()
    if RawBarcodeFile == '':
        sys.stderr.write("[warning] No raw barcode provided, will directly output barcode based on count.")
    else:
        RawBarcodeFileIn = open(RawBarcodeFile, "r", encoding='utf-8')
        RawBarcodeFlag = True
        LineCount = 0
        while True:
            line = RawBarcodeFileIn.readline().strip()
            if not line:
                break
            line = line.split()
            RawID = '|'.join(line)
            RawCellSet.add(RawID)
            LineCount = LineCount + 1
            if LineCount % 10000 == 0:
                sys.stderr.write("[processed] " + str(LineCount) + "\n")
        RawBarcodeFileIn.close()

    BarcodeFragmentList = []
    FragCountList = []
    IDList = []
    FIN = open(InputFile, "r", encoding='utf-8')
    LineCount = 0
    DoubletFilterList = []

    print("Reading fragments")
    FOUT = open(OutputPrefix + "Cell", "w")
    while True:
        line = FIN.readline().strip()
        if not line:
            break
        if line[0] == '#' or line[0] == '@':
            continue
        else:
            InfoArray = line.split("\t")
            if len(InfoArray) > skip + 1:
                BarcodeID = "|".join(InfoArray[0:skip])
                if RawBarcodeFlag:
                    if BarcodeID in RawCellSet:
                        FragCount = int(InfoArray[skip])
                        BarcodeEvent = TranspositionEvent()
                        BarcodeEvent.FragCount = FragCount
                        BarcodeEvent.Coord = line
                        FragCountList.append(FragCount)
                        BarcodeEvent.ID = LineCount
                        BarcodeFragmentList.append(BarcodeEvent)
                else:
                    FragCount = int(InfoArray[skip])
                    BarcodeEvent = TranspositionEvent()
                    BarcodeEvent.FragCount = FragCount
                    BarcodeEvent.Coord = line
                    FragCountList.append(FragCount)
                    BarcodeEvent.ID = LineCount
                    BarcodeFragmentList.append(BarcodeEvent)
                LineCount = LineCount + 1
        if LineCount % 10000 == 0:
            sys.stderr.write("[processed] " + str(LineCount) + "\n")
    print("Removing empty droplet")
    BarcodeCountCutoff = barcodeRanks(M=FragCountList)
    print(BarcodeCountCutoff)
    for Event in BarcodeFragmentList:
        if Event.FragCount >= BarcodeCountCutoff:
            InfoArray = Event.Coord.split("\t")
            for i in range(Event.FragCount):
                Event.Start.add(InfoArray[i * 3 + skip + 1] + "-" + InfoArray[i * 3 + skip + 2])
                Event.End.add(InfoArray[i * 3 + skip + 1] + "-" + InfoArray[i * 3 + skip + 3])
            DoubletFilterList.append(Event)
    print("Calculating similarity")
    print(len(DoubletFilterList))
    ComparePairList = combinations_with_replacement(DoubletFilterList, 2)
    Arr = np.empty((len(DoubletFilterList), len(DoubletFilterList)))
    CompareResult = pd.DataFrame(Arr, columns=range(len(DoubletFilterList)), index=range(len(DoubletFilterList)))
    if GenerateWholeTable:
        for ComparePair in ComparePairList:
            EventCount = compareTranspositionEvent(ComparePair)
            CompareResult.iloc[ComparePair[0].ID, ComparePair[1].ID] = EventCount
            CompareResult.iloc[ComparePair[1].ID, ComparePair[0].ID] = EventCount

        CompareResult.columns = IDList
        CompareResult.index = IDList
        CompareResult.to_csv(OutputPrefix + "PDLet.result.csv")

        for index, row in CompareResult.iterrows():
            if index == row.idxmax():
                FOUT.write(IDList[index] + "\n")


if __name__ == '__main__':
    main()
