import sys, os
import getopt
import time


try:
    opts, args = getopt.getopt(sys.argv[1:], "i:o:c:a:s", ["input=", "output=", "countcol=", "additional=", "stdout"])
except getopt.GetoptError:
    print(sys.argv[0] + "-i InputConfigFile -o OutputConfigFile -c CountColumn")
    sys.exit(2)

InputFile = ''
OutputPrefix = time.strftime("%b%d", time.localtime()).upper()
UserPrefix = ''
FragmentCountCol = ''
AddField = ''
StdOut = False

for opt, arg in opts:
    if opt == 'h':
        print(sys.argv[0] + "-i InputConfigFile -o OutputConfigFile")
        sys.exit()
    elif opt in ("-i", "--input"):
        InputFile = arg
    elif opt in ("-o", "--output"):
        UserPrefix = arg
    elif opt in ("-c", "--countcol"):
        FragmentCountCol = arg
    elif opt in ("-a", "--additional"):
        AddField = arg
    elif opt in ("-s", "--stdout"):
        StdOut = True

AddList = []
if AddField != '':
    AddList = AddField.split(',')

if UserPrefix != '':
    OutputPrefix = UserPrefix + ".bed"
else:
    OutputPrefix = OutputPrefix + ".bed"

if not StdOut:
    sys.stdout = open(OutputPrefix, "w")

with open(InputFile, "r", encoding="utf-8") as INFile:
    LineCount = 0
    while True:
        line = INFile.readline().strip()
        LineCount = LineCount + 1
        if not line:
            break
        if FragmentCountCol == '':
            if LineCount > 5:
                print("Count column do not specified or exists in the inline bed file")
                exit(2)
                break
            if line.startswith("@CN"):
                Content = line.split()
                for i in range(len(Content)):
                    if Content[i].lower() == "count":
                        FragmentCountCol = i - 1
        if line.startswith("#") or line.startswith("@"):
            continue

        CPXInfo = line.split("\t")
        FragmentCountCol = int(FragmentCountCol)
        Comment = CPXInfo[:FragmentCountCol]
        FragmentCount = int(CPXInfo[FragmentCountCol])
        Frags = CPXInfo[(FragmentCountCol + 1):]
        for i in range(FragmentCount):
            NewLine = Frags[i * 3] + "\t" + Frags[i * 3 + 1] + "\t" + Frags[i * 3 + 2] + "\t"
            for j in range(len(Comment)):
                NewLine = NewLine + Comment[j]
            for j in AddList:
                NewLine = NewLine + "\t" + Frags[3 * FragmentCount + (int(j) * i)]
            NewLine = NewLine
            print(NewLine)

        LineCount = LineCount + 1