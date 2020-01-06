from matplotlib import pyplot as plt
import sys
import json

ifile = open(sys.argv[1],"r")
WindowSize=100
BorderQ=0.1
MinBorder=50
MaxBorder=200

data=[]
Average=0.0
MRDRs={}
Ci=0
Chr=None
for line in ifile:
    if len(line)==0 or line[:2]=="##":
        continue
    if line[0]=="#":
        sl=line[1:].split()
        Chr=sl[0]
        #print("Reading %s..."%Chr,file=sys.stderr)
        Length=int(sl[1])
        MRDRs[Chr]=[0]*Length
        Ci=0
        continue
    MRDRs[Chr][Ci]=float(line)
    Ci+=1
ifile.close()
OtherData=None
SampleName=None
if len(sys.argv)>2:
    i=1
    while i< len(sys.argv):
        a=sys.argv[i]
        if a=="-SN":
            i+=1
            SampleName=sys.argv[i]
        else:
            ODFN=sys.argv[i]
        i+=1
    ODF=open(ODFN,"r")
    OtherData=json.load(ODF)
    ODF.close()
try:
    WindowSize=OtherData["RDWindowSize"]
except:
    pass

def show(data,WS,WE,Border=0,Title=None,OD=None,OtherData=None,Chr=None):
    #print(data[WS:WE])
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.vlines(range(WS,WE),[0]*(WE-WS),data[WS:WE],color="red")
    ax1.vlines(range(WS-Border,WS),[0]*Border,data[WS-Border:WS])
    ax1.vlines(range(WE,WE+Border),[0]*Border,data[WE:WE+Border])
    if Title!=None:
        ax1.set_title(Title)
    ax1.grid(True)
    if OD!=None:
        ax2=fig.add_subplot(212)
        data2=OD
        ax2.vlines(range(WS,WE),[0]*(WE-WS),data2[WS:WE],color="red")
        ax2.vlines(range(WS-Border,WS),[0]*Border,data2[WS-Border:WS])
        ax2.vlines(range(WE,WE+Border),[0]*Border,data2[WE:WE+Border])
        ax2.set_title("RD Standards")
        if SampleName!=None and OtherData!=None and Chr!=None:
            if "Segmentation" in OtherData.keys():
                Last=WS-Border
                LSeg=None
                for seg in OtherData["Segmentation"][Chr][SampleName]:
                    if WS-Border<=seg[0]<WE+Border:
                        ax1.plot((Last,seg[0]),(seg[1],seg[1]),color="blue")
                        Last=seg[0]
                        print(seg)
                        #ax2.vlines(seg[0],0,data2[seg[0]],color="blue")
                    elif seg[0]>=WE:
                        LSeg=seg
                        break
                if Last!=WS-Border:
                    if LSeg!=None:
                        print(LSeg)
                        ax1.plot((Last,WE+Border),(LSeg[1],LSeg[1]),color="blue")
    plt.show()

def help():
    print("""
-h                  show this help
show chr start end
  or chr:start-end  show the mrdr data
set  key=value      change parameters
     windowsize=100
q                   quit
"""
    ,file=sys.stderr)

while 1:
    print("Enter command:",end="",file=sys.stderr)
    Command=input()
    if Command=="q":
        exit(0)
    elif Command[:4]=="show":
        try:
            if ':' in Command:
                sc=[]
                sl=Command.split()
                sc.append(sl[0])
                sl=sl[1].split(":")
                sc.append(sl[0])
                sl=sl[1].split("-")
                sc+=sl
            else:
                sc=Command.split()
            Chr=sc[1]
            Start=int(sc[2])
            End=int(sc[3])
            WStart=int(Start/WindowSize)
            WEnd=int(End/WindowSize)
            WLength=WEnd-WStart
            Border=int(WLength*BorderQ)
            Border=max(MinBorder,Border)
            Border=min(MaxBorder,Border)
            print("Showing windows %s:%s-%s(WindowSize=%s)"%(Chr,WStart,WEnd,WindowSize),file=sys.stderr)
            OD=None
            if OtherData!=None:
                OD=OtherData["RDWindowStandards"][Chr]
            Title="MRDR data(%s:%s-%s)"%(Chr,Start,End)
            show(MRDRs[Chr],WStart,WEnd,Border,Title,OD,OtherData,Chr)
        except Exception as e:
            print(e,file=sys.stderr)
            help()
    elif Command[:3]=="set":
        for S in Command.split()[1:]:
            k=S.split("=")[0]
            if k.upper()=="WINDOWSIZE":
                v=int(S.split("=")[1])
                WindowSize=v
    else:
        help()