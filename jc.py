import pysam
import sys
import time
from globals import *

class Evidence:
    CombinePercentage=0.75
    def __init__(self):
        self.Data=[]
        self.Spread=(0,0)
    
    def addData(self,Type, Data):#Type: 0: DR cluster, 
        self.Data.append((Type,Data))
        self.calculateSpread()
    
    def calculateSpread(self):
        self.Spread=None
        for Datum in self.Data:
            if Datum[0]==0:
                if self.Spread==None:
                    self.Spread=(Datum[1].Start,Datum[1].End)
                else:
                    self.Spread=(min(self.Spread[0],Datum[1].Start),max(self.Spread[1],Datum[1].End))
    def combineEvidence(self, other):#return: 0: not combine, 1: combine, -1: not overlap
        if self.Spread[0]<other.Spread[1] and self.Spread[0]>other.Spread[0]:
            Overlap=min(other.Spread[1]-self.Spread[0],self.Spread[1]-self.Spread[0])
        elif self.Spread[1]>other.Spread[0] and self.Spread[1]<other.Spread[0]:
            Overlap=self.Spread[1]-other.Spread[0]
        elif other.Spread[0]>self.Spread[0] and other.Spread[0]<self.Spread[1]:#self covers other
            Overlap=other.Spread[1]-other.Spread[0]
        else:
            Overlap=0
            return -1
        MinLength=min(self.Spread[1]-self.Spread[0],other.Spread[1]-other.Spread[0])
        Overlap=Overlap/MinLength if MinLength!=0 else 0
        if Overlap>=CombinePercentage:
            self.Data+=other.Data
            self.Spread=(min(self.Spread[0],other.Spread[0]),max(self.Spread[1],other.Spread[1]))
            return 1
        return 0


class PairInfo:
    LMapped=0x1
    LReversed=0x2
    RMapped=0x4
    RReversed=0x8
    TransChrs=0x10
    Crossed=0x20#2 reads crosses another
    Contained=0x40#left read contains right one
    def hasFlags(self,Flags):
        return self.Flags & Flags
    def __init__(self, Read1:pysam.AlignedSegment, Read2:pysam.AlignedSegment):
        self.Flags=0#low to high:LMapped,LReversed,RMapped,RReversed,TransChrs,Crossed
        self.LStart=0
        self.LEnd=0
        self.RStart=0
        self.REnd=0
        self.NMapped=0
        if Read1.is_unmapped:
            if Read2.is_unmapped:
                pass
            else:
                R2S=Read2.reference_start+RefStartPos[Read2.tid]
                R2E=Read2.reference_end+RefStartPos[Read2.tid]
                self.LStart=R2S
                self.LEnd=R2E
                self.LLength=Read2.query_length
                self.RStart=R2S
                self.REnd=R2E
                self.RLength=Read2.query_length
                self.NMapped=1
                self.Flags|=PairInfo.LMapped|PairInfo.RMapped
                if Read2.is_reverse:
                    self.Flags|=PairInfo.LReversed|PairInfo.RReversed
        else:
            R1S=Read1.reference_start+RefStartPos[Read1.tid]
            R1E=Read1.reference_end+RefStartPos[Read1.tid]
            if Read2.is_unmapped:
                self.LStart=R1S
                self.LEnd=R1E
                self.LLength=Read1.query_length
                self.RStart=R1S
                self.REnd=R1E
                self.RLength=Read1.query_length
                self.NMapped=1
                self.Flags|=PairInfo.LMapped|PairInfo.RMapped
                if Read1.is_reverse:
                    self.Flags|=PairInfo.LReversed|PairInfo.RReversed
            else:
                R2S=Read2.reference_start+RefStartPos[Read2.tid]
                R2E=Read2.reference_end+RefStartPos[Read2.tid]
                if Read1.reference_name!=Read2.reference_name:
                    self.Flags|=PairInfo.TransChrs
                if R1S<R2S:
                    self.LStart=R1S
                    self.LEnd=R1E
                    self.LLength=Read1.query_length
                    self.RStart=R2S
                    self.REnd=R2S
                    self.RLength=Read2.query_length
                    if Read1.is_reverse:
                        self.Flags|=PairInfo.LReversed
                    if Read2.is_reverse:
                        self.Flags|=PairInfo.RReversed
                else:
                    self.LStart=R2S
                    self.LEnd=R2E
                    self.LLength=Read2.query_length
                    self.RStart=R1S
                    self.REnd=R1S
                    self.RLength=Read1.query_length
                    if Read1.is_reverse:
                        self.Flags|=PairInfo.RReversed
                    if Read2.is_reverse:
                        self.Flags|=PairInfo.LReversed
                if self.REnd<=self.LEnd:
                    self.Flags|=PairInfo.Contained
                if self.RStart<self.LEnd:
                    self.Flags|=PairInfo.Crossed
                self.NMapped=2
        self.Start=self.LStart
        self.End=max(self.LEnd,self.REnd)
        self.InsertionSize=self.End-self.Start

def isDiscordant(pi):
    if pi.NMapped!=2 or pi.InsertionSize<MedianInsertionSize-3*ISSD or pi.InsertionSize>MedianInsertionSize+3*ISSD:
        return True
    return False

def linked(Pair1, Pair2):
    if Pair1.NMapped==Pair2.NMapped:
        if Pair1.hasFlags(PairInfo.Crossed|PairInfo.Contained) or Pair2.hasFlags(PairInfo.Crossed|PairInfo.Contained):
            return False
        if abs(Pair1.Start-Pair2.Start)<3*ISSD and abs(Pair1.End-Pair2.End)<3*ISSD\
            and (Pair1.NMapped!=2 or (Pair1.NMapped==2 and Pair1.LEnd<=Pair2.RStart and Pair2.LEnd<=Pair1.RStart)):
            return True
        return False
    else:
        return False
def cluster(Reads):
    Clusters=[]
    Edges=[]
    print("building edges...",file=sys.stderr)
    for i in range(len(Reads)):
        Edges.append([])
        for j in range(i+1,len(Reads)):
            if linked(Reads[i],Reads[j]):
                Edges[i].append(j)
            if abs(Reads[i].Start-Reads[j].Start)>=3*ISSD:
                break
    print("edges built, clustering...",file=sys.stderr)
    Current=[]
    CurrentI=[]
    Clustered=[False]*len(Reads)
    CS=0
    CE=0
    for i in range(len(Reads)):
        if Clustered[i]:
            continue
        CurrentI.append(i)
        for k in CurrentI:
            for j in Edges[k]:
                if Clustered[j]:
                    continue
                NotLinked=False
                for l in range(len(CurrentI)):
                    if not linked(Reads[CurrentI[l]],Reads[j]):
                        NotLinked=True
                        break
                if NotLinked:
                    continue
                CurrentI.append(j)
                Clustered[j]=True
        for j in CurrentI:
            Current.append(Reads[j])
        CS=Current[0].Start
        CE=Current[0].End
        for r in Current:
            CS=min(CS,r.Start)
            CE=max(CE,r.End)
        Clusters.append((CS,CE,Current))
        Current=[]
        CurrentI=[]
    print("clustered.",file=sys.stderr)
    return Clusters

def getTidByCord(Cordinate):
    i=1
    while i<len(RefStartPos):
        if Cordinate<RefStartPos[i]:
            return i-1
        i+=1
    return i-1

def getScore(E):
    Score=0
    U=0
    N=0
    for d in E:
        if d[0]==0:
            for c in d[1]:
                for r in Cluster[2]:
                    U+=r.InsertionSize
                N+=len(c[2])
    return 0

def callSV(ReferenceFile,Cluster):
    SV=""
    U=0
    if Cluster[2][0].NMapped!=2:
        return ""
    if len(Cluster[2])<0:
        return ""
    for r in Cluster[2]:
        U+=r.InsertionSize
    U=(float(U)/len(Cluster[2])-MedianInsertionSize)/ISSD*float(len(Cluster[2]))**0.5#U~N(0,1)
    BKL=Cluster[0]
    BKR=Cluster[1]
    for r in Cluster[2]:
        BKL=max(BKL,r.LEnd)
        BKR=min(BKR,r.RStart)
    if U>10:
        SV+=ReferenceFile.references[getTidByCord(Cluster[0])]+":"+str(1+Cluster[0]-RefStartPos[getTidByCord(Cluster[0])])+"-"+ReferenceFile.references[getTidByCord(Cluster[1])]+":"+str(1+Cluster[1]-RefStartPos[getTidByCord(Cluster[1])])
        if Cluster[2][0].InsertionSize<MedianInsertionSize:
            SV+=", INS"
        else:
            SV+=", DEL"
        SV+=", Breakpoint:[%s,%s]"%(ReferenceFile.references[getTidByCord(BKL)]+":"+str(1+BKL-RefStartPos[getTidByCord(BKL)]),ReferenceFile.references[getTidByCord(BKR)]+":"+str(1+BKR-RefStartPos[getTidByCord(BKR)]))
    return SV

def makeDRCEvidences(DRClusters):
    Evidences=[]
    EI=-1
    for i in range(len(DRClusters)):
        if DRClusters[i]!=None:
            Temp=Evidence()
            Temp.addData(0,DRClusters[i])
            Evidences.append(Temp)
            EI+=1
        else:
            continue
        for j in range(i+1,len(DRClusters)):
            if DRClusters[j]!=None:
                Temp=Evidence()
                Temp.addData(0,DRClusters[j])
            else:
                continue
            CResult=Evidences[EI].combineEvidence(Temp)
            if CResult==-1:
                break
            elif CResult==1:
                DRClusters[j]=None
    return Evidences

ReferenceFile=pysam.FastaFile(sys.argv[1])
PosCount=0
for tid in range(ReferenceFile.nreferences):
    RefStartPos.append(PosCount)
    PosCount+=ReferenceFile.lengths[tid]
RefLength=PosCount
WindowSize=100
WindowsN=int(RefLength/WindowSize)+1
RDWindows=[]

SampleIndex=0
for i in range(2,len(sys.argv)):
    SamFile=pysam.AlignmentFile(sys.argv[i],"rb",reference_filename=sys.argv[1])
    RDWindows.append([0]*WindowsN)
    for read in SamFile:
        read: pysam.AlignedSegment
        if read.is_paired:
            if read.query_name in UnpairedReads:
                Pair=PairInfo(UnpairedReads.pop(read.query_name),read)
                if isDiscordant(Pair):
                    DiscordantReads.append(Pair)
                RDWindows[SampleIndex][int((Pair.Start+Pair.End)/2/100)]+=1
            else:
                UnpairedReads[read.query_name]=read
    SamFile.close()
    SampleIndex+=1
RDWindowsAverage=[0]*WindowsN
RDWindowsSum=[0]*WindowsN
SampleN=len(RDWindows)
SampleSum=[0]*SampleN
SampleAverage=[0]*SampleN
for i in range(WindowsN):
    for j in range(SampleN):
        RDWindowsSum[i]+=RDWindows[j][i]
        SampleSum[j]+=RDWindows[j][i]
    RDWindowsAverage[i]=RDWindowsSum[i]/SampleN
for j in range(SampleN):
    SampleAverage[j]=SampleSum[j]/WindowsN
DiscordantReads.sort(key=lambda r:r.Start)
DRClusters=cluster(DiscordantReads)

#Evidences=makeDRCEvidences()

#analyzeRD
#analyzeDR
#combineEvidence

for c in DRClusters:
    SV=callSV(ReferenceFile,c)
    if SV!="":
        print(SV)

ReferenceFile.close()