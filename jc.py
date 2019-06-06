import pysam
import sys
import time
from globals import *
class PairInfo:
    def __init__(self, Read1:pysam.AlignedSegment, Read2:pysam.AlignedSegment):
        self.TransChr=False
        if Read1.is_unmapped:
            if Read2.is_unmapped:
                self.Start=0
                self.End=0
                self.Interval=0
                self.NMapped=0
            else:
                R2S=Read2.reference_start+RefStartPos[Read2.tid]
                R2E=Read2.reference_end+RefStartPos[Read2.tid]
                self.Start=R2S
                self.End=R2E
                self.Interval=0
                self.NMapped=1
        else:
            R1S=Read1.reference_start+RefStartPos[Read1.tid]
            R1E=Read1.reference_end+RefStartPos[Read1.tid]
            if Read2.is_unmapped:
                self.Start=R1S
                self.End=R1E
                self.Interval=0
                self.NMapped=1
            else:
                R2S=Read2.reference_start+RefStartPos[Read2.tid]
                R2E=Read2.reference_end+RefStartPos[Read2.tid]
                if Read1.reference_name!=Read2.reference_name:
                    self.TransChr=True
                self.Start=min(R1S,R2S)
                self.End=max(R1E,R2E)
                self.Interval=self.End-self.Start#what if R1 expand over R2?
                self.NMapped=2

def analysePair(Read1: pysam.AlignedSegment, Read2: pysam.AlignedSegment):
    pi=PairInfo(Read1,Read2)
    if isDiscordant(pi):
        return pi
    else:
        return None

def isDiscordant(pi):
    if pi.NMapped!=2 or pi.Interval<AverageInsertSize-3*ISSD or pi.Interval>AverageInsertSize+3*ISSD:
        return True
    return False

def cluster(Reads):
    Clusters=[]
    i=0
    Current=[]
    CS=0
    CE=0
    while i<len(Reads):
        if Reads[i]==None:
            i+=1
            continue
        cr=Reads[i]
        if Current==[]:
            Current.append(cr)
            CS=cr.Start
            CE=cr.End
            Reads[i]=None
            i+=1
        else:
            if cr.Start>CS+3*ISSD or cr.End>CE+3*ISSD:
                Clusters.append((CS,CE,Current))
                Current=[]
                i=0
                continue
            if cr.End> CE:
                    CE=cr.End
            Current.append(cr)
            Reads[i]=None
            i+=1
    return Clusters

def getTidByCord(Cordinate):
    i=1
    while i<len(RefStartPos):
        if Cordinate<RefStartPos[i]:
            return i-1
        i+=1
    return i-1

def callSV(Cluster):
    SV=""
    if len(Cluster[2])>=5:
        SV+=SamFile.get_reference_name(getTidByCord(Cluster[0]))+":"+str(1+Cluster[0]-RefStartPos[getTidByCord(Cluster[0])])+"-"+SamFile.get_reference_name(getTidByCord(Cluster[1]))+":"+str(1+Cluster[1]-RefStartPos[getTidByCord(Cluster[1])])
        if Cluster[2][0].Interval<AverageInsertSize:
            SV+=", INS"
        else:
            SV+=", DEL"
    return SV

SamFile=pysam.AlignmentFile(sys.argv[1],"rb",reference_filename=sys.argv[2])
PosCount=0
for tid in range(SamFile.nreferences):
    RefStartPos.append(PosCount)
    PosCount+=SamFile.lengths[tid]
for read in SamFile:
    read: pysam.AlignedSegment
    if read.is_paired:
        if read.query_name in UnpairedReads:
            Pair=analysePair(UnpairedReads.pop(read.query_name),read)
            if Pair!=None:
                DiscordantReads.append(Pair)
        else:
            UnpairedReads[read.query_name]=read
print(len(DiscordantReads))
DiscordantReads.sort(key=lambda r:r.Start)
DRClusters=cluster(DiscordantReads)
for c in DRClusters:
    SV=callSV(c)
    if SV!="":
        print(SV)
SamFile.close()
