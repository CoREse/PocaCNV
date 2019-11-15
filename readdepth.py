import globals
from globals import *
from utils import *
import sys

class RDInterval:
    def __init__(self,WBegin,WEnd,ARD):
        self.WBegin=WBegin
        self.WEnd=WEnd
        self.AverageRD=ARD
        self.Begin=WBegin*RDWindowSize
        self.End=WEnd*RDWindowSize
        self.SupportedSVType=None#0:del, 1:insertion, 2:dup
        if self.AverageRD<0.1:
            self.CN=0
        elif self.AverageRD<0.75:
            self.CN=1
        else:
            self.CN=int(self.AverageRD+0.5)*2-1
        if self.CN==0 or self.CN==1:
            self.SupportedSVType=0
        elif self.CN>2:
            self.SupportedSVType=2

def sigDiff(RDRs,i,CurrentRunRatio):
    DupLine=1.4
    DelLine=0.75
    DupCut=0.4
    DelBackLine=0.95
    if CurrentRunRatio==1:
        if RDRs[i]>DupLine and RDRs[i-1]>DupLine:
            return True
        if RDRs[i]<DelLine and RDRs[i-1]<DelLine:
            return True
    elif CurrentRunRatio<1:
        if RDRs[i]>DelBackLine and RDRs[i-1]>DelBackLine:
            return True
    elif CurrentRunRatio>1:
        if abs(RDRs[i]-CurrentRunRatio)>DupCut and abs(RDRs[i-1]-CurrentRunRatio)>DupCut:
            return True
    return False


def makeRDIntervals(MixedRDRs):
    EvidenceIntervals=[]
    for i in range(len(MixedRDRs)):#for each sample
        CurrentRunStart=0
        CurrentRunRatio=1
        for j in range(1,len(MixedRDRs[i])):
            if sigDiff(MixedRDRs[i],j,CurrentRunRatio):
                if CurrentRunRatio!=1:
                    LastInterval=RDInterval(CurrentRunStart,j-1,CurrentRunRatio)
                    if LastInterval.CN!=2:
                        EvidenceIntervals.append(LastInterval)
                CurrentRunRatio=(MixedRDRs[i][j]+MixedRDRs[i][j-1])/2
                CurrentRunStart=j-1
            elif CurrentRunRatio!=1:
                CurrentRunRatio=(CurrentRunRatio*(j-CurrentRunStart)+MixedRDRs[i][j])/(j-CurrentRunStart+1)
        if CurrentRunRatio!=1:
            LastInterval=RDInterval(CurrentRunStart,len(MixedRDRs[i]),CurrentRunRatio)
            if LastInterval.CN!=2:
                EvidenceIntervals.append(LastInterval)
    return EvidenceIntervals

def makeRDICandidates(RDIs):
    Candidates=[]
    for i in RDIs:
        e=Evidence()
        e.setData(1,i)
        Candidates.append(Candidate([e]))
    Candidates=combineCandidateSets(Candidates,[])
    return Candidates

def partDiff(V1,V2):
    DiffRatio=1.3
    DiffValue=0.8
    if V1==0 or V2==0:
        if abs(V1-V2)>DiffValue:
            return True
    elif V1/V2>DiffRatio or V2/V1>DiffRatio or abs(V1-V2)>DiffValue:
        return True
    return False

def diverify(RDRs,SList,Previous):
    if len(SList)<1:
        return False
    for i in range(len(SList)):
        if i==0 or len(SList)==1:
            Mean=0
            if len(SList)==1:
                Mean=RDRs[0][SList[0]]
            else:
                for j in range(SList[0],SList[1]):
                    Mean+=RDRs[0][j]
                Mean/=SList[1]-SList[0]
            if not partDiff(Previous,Mean):
                return False
        elif i==len(SList)-1:
            Mean=0
            for j in range(SList[i-1],SList[i]):
                Mean+=RDRs[0][j]
            Mean/=SList[i]-SList[i-1]
            if not partDiff(Mean,RDRs[0][SList[i]]):
                return False
        else:
            PMean=0
            NMean=0
            for j in range(SList[i-1],SList[i]):
                PMean+=RDRs[0][j]
            PMean/=SList[i]-SList[i-1]
            if partDiff(PMean,RDRs[0][SList[i]]):
                return True
            for j in range(SList[i],SList[i+1]):
                NMean+=RDRs[0][j]
            NMean/=SList[i+1]-SList[i-1]
            if not partDiff(PMean,NMean):
                return False
    return True

def partition(RDRs):
    Separator=[]
    Previous=0
    LastI=0
    for i in range(1,len(RDRs[0])):
        #if partDiff(RDRs[0][i-1],RDRs[0][i]):
        Previous+=RDRs[0][i-1]
        if diverify(RDRs,[i],Previous/(i-LastI)):
            Separator.append(i)
            Previous=0
            LastI=i
    print(len(Separator), diverify(RDRs,Separator,0),file=sys.stderr)
    i=1
    Previous=0
    while i< len(Separator):
        Poped=Separator.pop(i)
        if diverify(RDRs,Separator,Previous):
            continue
        Separator.insert(i,Poped)
        i+=1
    print(len(Separator), diverify(RDRs,Separator,0),file=sys.stderr) 
    exit(0)

def analyzeRD(RDWindows,WindowsN,OccurredWindowsN,NormalizationOnly=False):
    RDWindowAverages=[0]*WindowsN
    #RDWindowSums=[0]*WindowsN
    SampleN=len(RDWindows)
    SampleSums=[0]*SampleN
    SampleAverages=RDWindows[-1]
    if SampleN<10:#WR,SR组合不太科学，假如几个样本同时有某个变异，那么很可能无法检测出这个变异，样本很多时不如直接用WR（基于变异占少数的假设），但假如变异本身不罕见就又有问题了
        for i in range(WindowsN):
            for j in range(SampleN):
                #RDWindowSums[i]+=RDWindows[j][i]
                SampleSums[j]+=RDWindows[j][i]
            #RDWindowAverages[i]=RDWindowSums[i]/SampleN
        for j in range(SampleN):
            SampleAverages[j]=SampleSums[j]/OccurredWindowsN
        MixedRDRs=[[0]*WindowsN]*SampleN#Mixed Read depth rate
        ZCount=0
        for i in range(OccurredWindowsN):
            if RDWindows[0][int(RefStartPos[RefInd["chr22"]]/100)+i]==0:
                ZCount+=1
        print(OccurredWindowsN, ZCount, SampleSums[0],SampleAverages[0])
        for i in range(SampleN):
            for j in range(WindowsN):
                SR=(RDWindows[i][j]/SampleAverages[i]) if SampleAverages[i]!=0 else 0
                MixedRDRs[i][j]=SR
                #WR=(RDWindows[i][j]/RDWindowAverages[j]) if RDWindowAverages[j]!=0 else 0
                #MixedRDRs[i][j]=(SR-WR)/SampleN+WR
    else:
        for i in range(WindowsN):
            for j in range(SampleN):
                #RDWindowSums[i]+=RDWindows[j][i]
                SampleSums[j]+=RDWindows[j][i]
            #RDWindowAverages[i]=RDWindowSums[i]/SampleN
        for j in range(SampleN):
            SampleAverages[j]=SampleSums[j]/OccurredWindowsN
        MixedRDRs=[[0]*WindowsN]*SampleN#Mixed Read depth rate
        for i in range(SampleN):
            for j in range(WindowsN):
                SR=(RDWindows[i][j]/SampleAverages[i]) if SampleAverages[i]!=0 else 0
                WR=(RDWindows[i][j]/RDWindowAverages[j]) if RDWindowAverages[j]!=0 else 0
                MixedRDRs[i][j]=(SR-WR)/SampleN+WR
    rdtestfile=open("data/rdtest.txt","w")
    for i in range(OccurredWindowsN):
        print("%d %f"%(i,MixedRDRs[0][i]),file=rdtestfile)
    rdtestfile.close()
    exit(0)
    ''' SIGMA(MixedRDRs)=NM, so, ERD=NM/NM=1
    ERD=0
    for i in range(WindowsN):
        for j in range(SampleN):
            ERD+=MixedRDRs[j][i]
    ERD/=OccurredWindowsN*SampleN
    #standardization
    for i in range(SampleN):
        for j in range(WindowsN):
            MixedRDRs[i][j]/=ERD
    ERD=1.0
    '''
    globals.ERD=1.0#ERD
    globals.MixedRDRs=MixedRDRs
    if NormalizationOnly:
        return MixedRDRs
        
    #partition(MixedRDRs)

    RDICandidates=makeRDICandidates(makeRDIntervals(MixedRDRs))
    return RDICandidates

import pysam
def writeRDData(RDWindows,ReferenceFile,SampleNames):
    ReferenceFile:pysam.FastaFile
    for i in range(len(SampleNames)):
        rdfile=open("data/rd%s.txt"%(SampleNames[i]),"w")
        print("#Length:%d"%(len(RDWindows[0])),end="",file=rdfile)
        Chri=0
        Chr=0
        for j in range(len(RDWindows[i])):
            if getTidByCord(j*100)!=Chr:
                Chr=getTidByCord(j*100)
                Chri=0
            print("\n%s %s %s"%(ReferenceFile.references[Chr], Chri, RDWindows[i][j]),end="",file=rdfile)
            Chri+=1
        rdfile.close()
    return

def readRDData(RDWindows, SampleNames, FileName):
    SampleName=FileName.split("\\")[-1].split("/")[-1][2:-4]
    SampleNames.append(SampleName)
    DataFile=open(FileName,"r")
    Length=0
    for line in DataFile:
        if line[0]=='#':
            if line.split(":")[0]=="#Length":
                Length=int(line.split(":")[1])
        break
    DataFile.seek(0)
    if Length==0:
        RDWindows.append([])
        for line in DataFile:
            if line[0]=='#':
                continue
            RDWindows[-1].append(int(line.split()[-1]))
    else:
        RDWindows.append([0]*Length)
        i=0
        for line in DataFile:
            if line[0]=='#':
                continue
            RDWindows[-1][i]=int(line.split()[-1])
            i+=1
            if i>=Length:
                break
    DataFile.close()
    return