import pysam
from utils import *
from readpair import PairInfo
import globals as g

def readExcludedAreas(DataFile,RefFile):
    ExcludedAreasByContig={}
    for tid in range(RefFile.nreferences):
        ExcludedAreasByContig[RefFile.references[tid]]=[]
    for line in DataFile:
        sl=line.split()
        if len(sl)!=4:
            continue
        Contig=sl[0]
        Start=int(sl[1])
        End=int(sl[2])
        ExcludedAreasByContig[Contig].append((Start,End))
    return ExcludedAreasByContig
        
def filtExcludedAreas(Intervals,ExcludedAreasByContig,Contig):
    FilteredIntervals=[]
    ExcludedAreas=ExcludedAreasByContig[Contig.Name]
    for Int in Intervals:
        Flag=True
        SE=getSE(Int)
        for EA in ExcludedAreas:
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