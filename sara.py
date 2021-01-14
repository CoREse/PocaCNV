from array import array
import math

def D(Data,h,t):#t<len(Data)-1
    Left=max(0,t-h+1)
    LC=t-Left+1
    Right=min(len(Data),t+h+1)
    RC=Right-t-1
    if LC!=h or RC!=h:
        return 0
    LSum=0
    RSum=0
    for k in range(Left, t+1):
        LSum+=Data[k]
    for k in range(t+1,Right):
        RSum+=Data[k]
    return (LSum-RSum)/h

def getScans(Data,h):#faster way
    Scans=array("d",[0]*(len(Data)))
    LSum=0
    RSum=0
    for k in range(0,h):
        LSum+=Data[k]
    for k in range(h,h+h):
        RSum+=Data[k]
    Scans[h-1]=(LSum-RSum)/h
    for i in range(h,len(Data)-h):
        Scans[i]=(h*Scans[i-1]+2*Data[i]-Data[i-h]-Data[i+h])/h
    for i in range(len(Scans)):
        Scans[i]=abs(Scans[i])
    return Scans


def getLocalMaximizers(Data,h=None,Scan=D):#Maximizer cut to [t-h+1,t],[t+1,t+h], transform into (End,Ave) should use t+1, default h=log(n)
    if h==None:
        h=int(math.log2(len(Data)))
    #Scans=array("d",[0]*(len(Data)))
    #for i in range(len(Scans)-1):
    #    Scans[i]=abs(Scan(Data,h,i))
    Scans=getScans(Data,h)
    LocalMaximizers=[]
    LastMaxI=-1
    for i in range(h-1,len(Data)-h-1):
        Max=-1
        if LastMaxI>i-h:
            if Scans[i+h]>Scans[LastMaxI]:
                MaxI=i+h
                Max=Scans[MaxI]
                LastMaxI=MaxI
                LocalMaximizers.append((MaxI,Max))
            else:
                MaxI=LastMaxI
                Max=Scans[MaxI]
            continue
        for t in range(i-h+1,i+h+1):
            if Scans[t]>Max:
                Max=Scans[t]
                MaxI=t
        if LastMaxI!=MaxI:
            LastMaxI=MaxI
            LocalMaximizers.append((MaxI,Max))
    return LocalMaximizers

def getSigmaSquareEstimation(Data,CutOffs):
    SquareSum=0
    Last=0
    for End,Ave in CutOffs:
        for i in range(Last,End):
            SquareSum+=Data[i]**2
    return SquareSum/len(Data)

def SaRa(Data,h=None,lmda=0.3):
    LocalMaximizers=getLocalMaximizers(Data,h)
    #LocalMaximizers.sort(key=lambda p:p[1],reverse=True)
    CutOffs=[]
    for i in range(len(LocalMaximizers)):
        if LocalMaximizers[i][1]>lmda:
            CutOffs.append((LocalMaximizers[i][0]+1,0))
    return CutOffs

#Data=[1,2,3,4,5,6,7,2,3,4,1,1,3,4,6,2,4,7,9,3,1,2,4,6,7,4,3,2,1,2,3,4,5,6,7,2,3,4,1,1,3,4,6,2,4,7,9,3,1,2,4,6,7,4,3,2,1,2,3,4,5,6,7,2,3,4,1,1,3,4,6,2,4,7,9,3,1,2,4,6,7,4,3,2,1,2,3,4,5,6,7,2,3,4,1,1,3,4,6,2,4,7,9,3,1,2,4,6,7,4,3,2,1,2,3,4,5,6,7,2,3,4,1,1,3,4,6,2,4,7,9,3,1,2,4,6,7,4,3,2]
#h=10
#print(getLocalMaximizers(Data,h))