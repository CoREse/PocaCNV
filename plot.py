from matplotlib import pyplot
import sys

ifile = open(sys.argv[1],"r")
data=[]
Average=0.0
for line in ifile:
    if len(line)==0 or line[0]=="#":
        continue
    data.append(float(line.split()[-1]))
    Average+=data[-1]
Average/=len(data)
D=0.0
Pos=[]
for i in range(len(data)):
    D+=(data[i]-Average)**2
    if data[i]>10:
        data[i]=10
    Pos.append(i+1)
D/=len(data)
print("A:%s, D:%s"%(Average,D))
pyplot.plot(data)
#pyplot.plot(Pos,data,".")
#pyplot.bar(len(data),data)
pyplot.show()