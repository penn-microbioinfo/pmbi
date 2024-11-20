import sys

compl = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        }

with open(sys.argv[1], 'r') as ss:
    num = 0
    for line in ss:
        if num < 2:
            print(line.strip())
        else:
            spl = line.strip().split(",")
            rc = ""
            for char in spl[2]:
                rc+=compl[char]
            spl[2] = rc[::-1]
            print(",".join(spl))
        num+=1


