import sys
import getopt
import json
import time
import gzip

# Used to generate customized config file.

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:w:", ["input=", "output=", "whitelist="])
except getopt.GetoptError:
    print(sys.argv[0] + "-i InputConfigFile -w Whitelist -o OutputConfigFile")
    sys.exit(2)

InputFile = ''
OutputPrefix = time.strftime("%b%d", time.localtime()).upper()
UserPrefix = ''
WhiteBarcodeList=''
for opt, arg in opts:
    if opt == 'h':
        sys.exit()
    elif opt in ("-i", "--input"):
        InputFile = arg
    elif opt in ("-o", "--output"):
        UserPrefix = arg
    elif opt in ("-w", "--whitelist"):
        WhiteBarcodeList = arg

if UserPrefix != '':
    OutputPrefix = UserPrefix + "_config.json"
else:
    OutputPrefix = OutputPrefix + "_config.json"

with open(InputFile, "r", encoding="utf-8") as INFile:
    content = INFile.read()

ConfigJson = json.loads(content)
BarcodeType = ConfigJson["barcode type"]

print("Generating configuration file for " + OutputPrefix + " from " + WhiteBarcodeList)
Barcode = "BC"
if WhiteBarcodeList.endswith(".gz"):
    BarcodeFile = gzip.open(WhiteBarcodeList, "r")
    NeedDecode = True
else:
    BarcodeFile = open(WhiteBarcodeList, "r", encoding='utf-8')
    NeedDecode = False

WhiteList = []
while True:
    if NeedDecode:
        Line = BarcodeFile.readline().decode('utf-8')
    else:
        Line = BarcodeFile.readline()
    if not Line:
        break
    BarcodeContent = Line.strip().split()
    if len(BarcodeContent) == 1:
        WhiteList.append(BarcodeContent[0])
    elif len(BarcodeContent) == 2:
        WhiteList.append({BarcodeContent[0]: BarcodeContent[1]})
ConfigJson["barcode type"][Barcode].update({"WHITE LIST": list(WhiteList)})



with open(OutputPrefix, "w") as FOUT:
    json.dump(ConfigJson, FOUT, indent=4)


