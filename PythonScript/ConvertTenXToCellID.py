import sys

Column = 0

RNA_IDB = open(sys.argv[1], "r")
Barcode_File = open(sys.argv[2], "r")

if len(sys.argv) > 3:
    Column = int(sys.argv[3])
FOUT = open("ToCellID.txt", "w")

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

while True:
    line = Barcode_File.readline().strip()
    if not line:
        break
    line = line.split()
    RawBarcode = line[Column]
    print(RawBarcode)
    if RNA_BarList.get(RawBarcode):
        CellID = RNA_BarList.get(RawBarcode)
    else:
        print("Barcode do not appear in IDB file")

