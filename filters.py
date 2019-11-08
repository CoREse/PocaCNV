import pysam
from globals import *
from utils import *
from readpair import PairInfo

ExcludedAreas=[]
ExcludedAreasByContig={}
def readExcludedAreas(DataFile,ContigNameOccurred=None):
    for i in range(len(RefInd)):
        ExcludedAreasByContig[i]=[]
    for line in DataFile:
        sl=line.split()
        if len(sl)!=4:
            continue
        Contig=sl[0]
        try:
            Start=int(sl[1])+RefStartPos[RefInd[Contig]]
            End=int(sl[2])+RefStartPos[RefInd[Contig]]
            ExcludedAreasByContig[RefInd[Contig]].append((Start,End))
            if ContigNameOccurred!=None:
                if not ContigNameOccurred[Contig]:
                    continue
        except:
            continue
        ExcludedAreas.append((Start,End))
        
def filtExcludedAreas(Intervals):
    FilteredIntervals=[]
    for Int in Intervals:
        Flag=True
        SE=getSE(Int)
        for EA in ExcludedAreasByContig[getTidByCord((SE[0]+SE[1])/2)]:
            if inclusion(EA,getSE(Int))>0:
                Flag=False
                break
        if Flag:
            FilteredIntervals.append(Int)
    return FilteredIntervals

def getSE(Interval):
    if type(Interval)==PairInfo:
        return (Interval.Start,Interval.End)
    return (Interval.Begin,Interval.End)