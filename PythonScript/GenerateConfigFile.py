import sys
import getopt
import json
import time

# Used to generate customized config file.

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:", ["input=", "output="])
except getopt.GetoptError:
    print(sys.argv[0] + "-i InputConfigFile -o OutputConfigFile")
    sys.exit(2)

InputFile = ''
OutputPrefix = time.strftime("%b%d", time.localtime()).upper()
UserPrefix = ''
for opt, arg in opts:
    if opt == 'h':
        print(sys.argv[0] + "-i InputConfigFile -o OutputConfigFile")
        sys.exit()
    elif opt in ("-i", "--input"):
        InputFile = arg
    elif opt in ("-o", "--output"):
        UserPrefix = arg

if UserPrefix != '':
    OutputPrefix = UserPrefix + "_config.json"
else:
    OutputPrefix = OutputPrefix + "_config.json"
print(InputFile, UserPrefix, OutputPrefix)

with open(InputFile, "r", encoding="utf-8") as INFile:
    content = INFile.read()

ConfigJson = json.loads(content)

BarcodeType = ConfigJson["barcode type"]


for Barcode in BarcodeType:
    if BarcodeType[Barcode].get("WHITE LIST") is not None:
        WhiteList = []
        if BarcodeType[Barcode].get("WHITE LIST") == '':
            WhiteList.append("NNNNN")
            print("Barcode " + Barcode + "not specified, will be treated with NNNNN")
        else:
            with open(BarcodeType[Barcode].get("WHITE LIST")) as BarcodeFile:
                while True:
                    Line = BarcodeFile.readline()
                    if not Line:
                        break
                    BarcodeContent = Line.strip().split()
                    if len(BarcodeContent) == 1:
                        WhiteList.append(BarcodeContent[0])
                    elif len(BarcodeContent) == 2:
                        WhiteList.append({BarcodeContent[0]: BarcodeContent[1]})
        # print("WHITE : ")
        # print(WhiteList)
        ConfigJson["barcode type"][Barcode].update({"WHITE LIST": list(WhiteList)})

with open(OutputPrefix, "w") as FOUT:
    json.dump(ConfigJson, FOUT, indent=4)