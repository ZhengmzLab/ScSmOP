import sys

RNA_IDB = open(sys.argv[1], "r")
ATAC_IDB = open(sys.argv[2], "r")

ATAC_File = open(sys.argv[3], "r")

Column = sys.argv[4]

TranslatedFile = open(sys.argv[3] + "Syn", "w")

RNA_BarList = {}
while True:
    line = RNA_IDB.readline().strip()
    if not line:
        break
    line = line.split()
    RNABarcode = line[0]
    RNA_Cell = "CELL_" + line[1]

    # print(RNABarcode)
    if not RNA_BarList.get(RNABarcode):
        RNA_BarList.update({RNABarcode: RNA_Cell})
    else:
        PreCell = RNA_BarList.get(RNABarcode)
        if PreCell != RNA_Cell:
            print("Same barcode assigned to different cell\n")


ATAC_BarList = {}
TranslateList = {}
while True:
    line = ATAC_IDB.readline().strip()
    if not line:
        break
    line = line.split()
    ATACBarcode = line[0]
    ATAC_Cell = "CELL_" + line[1]

    if RNA_BarList.get(ATACBarcode):
        TranslateList.update({ATAC_Cell: RNA_BarList.get(ATACBarcode)})

while True:
    line = ATAC_File.readline().strip()
    if not line:
        break
    line = line.split()
    LineColumn = len(line)
    ACell = line.pop(int(Column) - 1)
    if not TranslateList.get(ACell):
        # print(ACell + " not presented in RNA-seq. Output as ATAC-" + ACell)
        BCell = "ATAC-" + ACell
    else:
        BCell = TranslateList.get(ACell)
    line.insert(int(Column) - 1, BCell)
    BLine = "\t".join(line)
    TranslatedFile.write(BLine + "\n")

