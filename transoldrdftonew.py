import sys
ThisContig=None
ThisData=[]
first=True
with open(sys.argv[1]) as infile:
    ThisContig=None
    ThisData=[]
    for line in infile:
        if line[0]=="#":
            continue
        sl=line.split()
        ContigName=sl[0]
        if ThisContig==None:
            ThisContig=ContigName
        data=sl[2]
        if ContigName!=ThisContig:
            if not first:
                print("\n",end="")
            first=False
            print("#%s %s"%(ThisContig,len(ThisData)),end="")
            for i in range(len(ThisData)):
                print("\n%s"%(ThisData[i]),end="")
            ThisData=[]
            ThisData.append(data)
            ThisContig=ContigName
        else:
            ThisData.append(data)

if not first:
    print("\n",end="")
first=False
print("#%s %s"%(ThisContig,len(ThisData)),end="")
for i in range(len(ThisData)):
    print("\n%s"%(ThisData[i]),end="")