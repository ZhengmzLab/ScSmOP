import sys
import getopt


PeakList = {}

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:s", ["input=", "output=", "countcol=", "stdout"])
except getopt.GetoptError:
    print(sys.argv[0] + "-i InputConfigFile -o OutputConfigFile")
    sys.exit(2)

InputFile = ''
UserPrefix = ''
OutputPrefix = ''

for opt, arg in opts:
    if opt == 'h':
        print(sys.argv[0] + "-i InputConfigFile -o OutputConfigFile -s")
        sys.exit()
    elif opt in ("-i", "--input"):
        InputFile = arg
    elif opt in ("-o", "--output"):
        UserPrefix = arg
    elif opt in ("-s", "--stdout"):
        StdOut = True


if UserPrefix != '':
    OutputPrefix = UserPrefix + ".bed"

if OutputPrefix != '':
    sys.stdout = open(OutputPrefix, "w")

if InputFile == '':
    FI = sys.stdin
else:
    FI = open(InputFile, "r")

while True:
    line = FI.readline().strip()
    if not line:
        break

    line = line.split()
    CurPeak = line[5] + "-" + line[6] + "-" + line[7]
    CurCell = line[3]

    if not PeakList.get(CurPeak):
        PeakList.update({CurPeak: dict(PeakCoord=[line[5], line[6], line[7]], CellList={CurCell: 1})})
    else:
        PeakCellCount = PeakList.get(CurPeak)
        if not PeakCellCount["CellList"].get(CurCell):
            PeakCellCount["CellList"].update({CurCell: 1})
        else:
            PeakCellCount["CellList"][CurCell] += 1

for Peak in PeakList.values():
    CurPeakCoord = '\t'.join(Peak["PeakCoord"])
    for Cell, Count in Peak["CellList"].items():
        print(CurPeakCoord + "\t" + Cell + "\t" + str(Count))
