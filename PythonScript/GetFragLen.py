import sys

FI = open(sys.argv[1], "r")
FOUT = open(sys.argv[1] + "FragLen.stat", "w")

Count = 0
while True:
    line = FI.readline()
    if not line:
        break
    if line.startswith("@") or line.startswith("#"):
        continue
    word = line.split("\t")
    Complex = word[0]
    FragCount = int(word[1])
    for i in range(FragCount):
        NewLine = word[i * 3 + 2] + "\t" + word[i * 3 + 3] + "\t" + word[i * 3 + 4] + "\t" + Complex + "-" + str(
            FragCount) + "-" + str(i) + "\t" + str(int(word[i * 3 + 4]) - int(word[i * 3 + 3])) + "\n"
        FOUT.writelines(NewLine)
    Count = Count + 1
    if Count % 10000 == 0:
        sys.stderr.write("[Processed] " + str(Count) + "\n")
