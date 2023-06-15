import sys

IDBFile = open(sys.argv[1], "r")
ConfigFile = open(sys.argv[2], "r")

FIN = open(sys.argv[3], "r")
FOUT = open("Out", "w")

BarcodeOrder = 0
BarcodeSet = {}
while True:
    line = ConfigFile.readline().strip()
    if not line:
        break

    line = line.split()[0]
    if not BarcodeSet.get(line):
        BarcodeSet.update({line: BarcodeOrder})
        BarcodeOrder = BarcodeOrder + 1

BarcodeOrderIDB = {}
while True:
    line = IDBFile.readline().strip()
    if not line:
        break

    line = line.split()
    BarcodeOrderInIDB = int(line[0])
    CellInIDB = "CELL_" + line[1]
    if not BarcodeOrderIDB.get(BarcodeOrderInIDB):
        BarcodeOrderIDB.update({BarcodeOrderInIDB: CellInIDB})

while True:
    line = FIN.readline().strip()
    if not line:
        break

    line = line.split()[0]
    if not BarcodeSet.get(line):
        print("Barcode " + line + " not presented in Config file.")
        LineOut = "-\n"
    else:
        Order = BarcodeSet.get(line)
        if not BarcodeOrderIDB.get(Order):
            print("Barcode " + line + " not presented in IDB file.")
            LineOut = "-\n"
        else:
            Cell = BarcodeOrderIDB.get(Order)
            LineOut = Cell + "\n"
    FOUT.write(LineOut)

