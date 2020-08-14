from utils import *
import globals as g
import pysam
import sys

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
    def __init__(self, Read1:pysam.AlignedSegment, Read2:pysam.AlignedSegment,SampleIndex):
        self.Flags=0#low to high:LMapped,LReversed,RMapped,RReversed,TransChrs,Crossed
        self.LStart=0
        self.LEnd=0
        self.RStart=0
        self.REnd=0
        self.NMapped=0
        self.SampleIndex=SampleIndex
        if Read1.is_unmapped:
            if Read2.is_unmapped:
                pass
            else:
                R2S=Read2.reference_start+g.RefStartPos[Read2.tid]
                R2E=Read2.reference_end+g.RefStartPos[Read2.tid]
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
            R1S=Read1.reference_start+g.RefStartPos[Read1.tid]
            R1E=Read1.reference_end+g.RefStartPos[Read1.tid]
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
                R2S=Read2.reference_start+g.RefStartPos[Read2.tid]
                R2E=Read2.reference_end+g.RefStartPos[Read2.tid]
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
                self.Flags|=PairInfo.LMapped|PairInfo.RMapped
        self.Start=self.LStart
        self.End=max(self.LEnd,self.REnd)
        self.InsertionSize=self.End-self.Start
    
    def isDiscordant(pi,conditionF=None):
        if conditionF==None:
            if pi.NMapped!=2 or pi.InsertionSize<MedianInsertionSize-3*g.ISSD or pi.InsertionSize>MedianInsertionSize+3*g.ISSD:
                return True
        else:
            return conditionF(pi)
        return False

    def linked(Pair1, Pair2):
        if Pair1.classify()!=Pair2.classify():
            return False
        if Pair1.NMapped==Pair2.NMapped:
            if Pair1.hasFlags(PairInfo.Crossed|PairInfo.Contained) or Pair2.hasFlags(PairInfo.Crossed|PairInfo.Contained):
                return False
            if abs(Pair1.Start-Pair2.Start)<3*g.ISSD and abs(Pair1.End-Pair2.End)<3*g.ISSD\
                and (Pair1.NMapped!=2 or (Pair1.NMapped==2 and Pair1.LEnd<=Pair2.RStart and Pair2.LEnd<=Pair1.RStart)):
                return True
            return False
        else:
            return False
    
    def classify(self):#return: 0: normal, 1: single mapped, 2: no mapped, 3: short insertion size, 4 long insertion size
        if self.NMapped==1:
            return 1
        elif self.NMapped==0:
            return 2
        elif self.InsertionSize<MedianInsertionSize-3*g.ISSD:
            return 3
        elif self.InsertionSize>MedianInsertionSize+3*g.ISSD:
            return 4
        return 0

def cluster(Reads):
    Clusters=[]
    Edges=[]
    print("building edges...",file=sys.stderr)
    for i in range(len(Reads)):
        Edges.append([])
        for j in range(i+1,len(Reads)):
            if Reads[i].linked(Reads[j]):
                Edges[i].append(j)
            if abs(Reads[i].Start-Reads[j].Start)>=3*g.ISSD:
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
                    if not Reads[CurrentI[l]].linked(Reads[j]):
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

def makeDRCCandidates(DRClusters):
    Candidates=[]
    CI=-1
    for i in range(len(DRClusters)):
        if DRClusters[i]!=None:
            NewE=Evidence(0,DRClusters[i][2],DRClusters[i][0],DRClusters[i][1])
            Temp=Candidate([NewE])
            Candidates.append(Temp)
            CI+=1
        else:
            continue
        for j in range(i+1,len(DRClusters)):
            if DRClusters[j]!=None:
                NewE=Evidence(0,DRClusters[j][2],DRClusters[j][0],DRClusters[j][1])
                Temp=Candidate([NewE])
            else:
                continue
            CResult=Candidates[CI].mergeCandidate(Temp)
            if CResult==-1:
                break
            elif CResult==1:
                DRClusters[j]=None
    return Candidates