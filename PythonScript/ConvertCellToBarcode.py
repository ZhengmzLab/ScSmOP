import sys

Column = 0

RNA_IDB = open(sys.argv[1], "r")
Barcode_File = open(sys.argv[2], "r")

if len(sys.argv) > 3:
    Column = int(sys.argv[3])
FOUT = open("ToBarcode.txt", "w")

RNA_BarList = {}

while True:
    line = RNA_IDB.readline().strip()
    if not line:
        break
    line = line.split()
    RNABarcode = line[0]
    RNA_Cell = "CELL_" + line[1]

    