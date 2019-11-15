from matplotlib import pyplot
import sys

ifile = open(sys.argv[1],"r")
data=[]
for line in ifile:
    if len(line)==0 or line[0]=="#":
        continue
    data.append(float(line.split()[-1]))
pyplot.plot(data)
pyplot.show()