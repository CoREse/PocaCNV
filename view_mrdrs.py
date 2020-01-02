from matplotlib import pyplot as plt
import sys

ifile = open(sys.argv[1],"r")
WindowSize=100
BorderQ=0.1
MinBorder=10
MaxBorder=100

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

def show(data,WS,WE):
    #print(data[WS:WE])
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.vlines(range(WS,WE),[0]*(WE-WS),data[WS:WE])
    ax1.grid(True)
    ax1.axhline(0, color='black',lw=2)
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
            Start=int(int(sc[2])/WindowSize)
            End=int(int(sc[3])/WindowSize)
            Length=End-Start
            Border=int(Length*BorderQ)
            Border=max(MinBorder,Border)
            Border=min(MaxBorder,Border)
            print("Showing windows %s:%s-%s(WindowSize=%s)"%(Chr,Start,End,WindowSize),file=sys.stderr)
            show(MRDRs[Chr],Start-Border,End+Border)
        except:
            help()
    elif Command[:3]=="set":
        for S in Command.split()[1:]:
            k=S.split("=")[0]
            if k.upper()=="WINDOWSIZE":
                v=int(S.split("=")[1])
                WindowSize=v
    else:
        help()