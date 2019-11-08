from matplotlib import pyplot

ifile = open("data/rdtest.txt","r")
data=[]
for line in ifile:
    if len(line)==0 or line[0]=="#":
        continue
    data.append(float(line.split()[2]))
pyplot.plot(data)
pyplot.show()