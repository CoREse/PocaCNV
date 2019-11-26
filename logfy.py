import sys
import math

ifile=open(sys.argv[1],"r")

data=[]
for line in ifile:
    data.append(float(line.split()[-1]))

data2=data.copy()
data2.sort()
median=data2[int(len(data2)/2)]

first=True
for d in data:
    if not first:
        print("\n",end="")
    first=False
    if d==0:
        d=1e-3
    print("%s"%(math.log2(d/median)),end="")

ifile.close()