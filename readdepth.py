import globals as g
from utils import *
import sys
from contig import *
import statistics
from array import array
import rpy2.robjects as robjects
from scipy.stats import poisson
import gc
import math
import multiprocessing as mp
import subprocess
from sara import SaRa
from drpprocessing import *
from cllprocessing import *

#import pyximportcpp; pyximportcpp.install()
from rdprocessing import *

Deploid=set(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"\
    ,"CHR1","CHR2","CHR3","CHR4","CHR5","CHR6","CHR7","CHR8","CHR9","CHR10","CHR11","CHR12","CHR13","CHR14","CHR15","CHR16","CHR17","CHR18","CHR19","CHR20","CHR21","CHR22"])
Haploid=set(["Y","CHRY"])
def ploidyProcessing(TheContig):
    if TheContig.Name.upper() in Deploid:
        return
    if TheContig.Name.upper() in Haploid:
        for i in range(len(TheContig.RDWindows)):
            TheContig.Ploidies[i]=1
        return
    for i in range(len(TheContig.RDWindows)):
        if TheContig.ContigSampleReadCounts[i]<g.SampleReadCount[i]/TheContig.Genome.GenomeLength*TheContig.NLength*0.75:
            TheContig.Ploidies[i]=1

def cn2filter(Interval,TheContig,Confidence=None):
    if Interval.mu==None:
        mu,mus=Interval.calcMuMus(TheContig)
    else:
        mu,mus=(Interval.mu,Interval.mus)
    if Confidence==None:
        Confidence=g.Parameters.CN2FilterConfidence
    qint=poisson.interval(Confidence,mu)
    if qint[0]<mus<qint[1]:
        return False
    return True

def cn2likely(Interval,TheContig,SP=None):
    if Interval.mu==None:
        mu,mus=Interval.calcMuMus(TheContig,SP)
    else:
        mu,mus=(Interval.mu,Interval.mus)
    cd=poisson.cdf(mus,mu)
    return min(cd,1-cd)
def mulikely(mu,mus):#在mu的情况下，出现小于mus或者大于mus的概率
    cd=poisson.cdf(mus,mu)
    return min(cd,1-cd)
"""
def getSampleSum(TheContig, SampleI, WBegin, WEnd):
    SampleRD=0
    for j in range(WBegin,WEnd):
        SampleRD+=TheContig.RDWindows[SampleI][j]
    return SampleRD

#should be irrelevant to ploidy if we use local, since all windows in the same sample have the same ploidies.
def getSP(TheContig, WBegin, WEnd, NSD=3, MinimumTake=0.8, local=g.StatLocal):#get rd sum and sample read count sum
    SRS=0
    SRC=0
    if WEnd<=WBegin:
        return (0,0)
    SampleN=len(TheContig.SampleNames)
    SampleRDs=[0]*SampleN
    Stat=g
    if local:
        Stat=TheContig
    #start with middle p
    #this seems better
    '''SamplePs=[0]*SampleN
    for i in range(SampleN):
        SampleRDs[i]=getSampleSum(TheContig,i,WBegin,WEnd)
        SamplePs[i]=SampleRDs[i]/Stat.SampleReadCount[i]
    SamplePs.sort()
    EstimatedP=SamplePs[int(SampleN/2)]
    SRS=EstimatedP
    SRC=1
    if SampleN%2==0:
        EstimatedP+=SamplePs[int(SampleN/2)-1]
        EstimatedP/=2'''
    #start with whole
    for i in range(SampleN):
        SampleRDs[i]=getSampleSum(TheContig,i,WBegin,WEnd)
        SRS+=SampleRDs[i]
        SRC+=Stat.SampleReadCount[i]
    EstimatedP=SRS/SRC
    RemovedSet=set()
    LastP=0
    SampleSTDs=[0]*SampleN
    while LastP!=EstimatedP:
        RemovedSet=set()
        for i in range(SampleN):
            SampleSTDs[i]=EstimatedP*Stat.SampleReadCount[i]
            if abs(SampleRDs[i]-SampleSTDs[i])>NSD*SampleSTDs[i]**0.5:
                RemovedSet.add(i)
        if len(RemovedSet)>SampleN*MinimumTake:
            break
        SRS=0
        SRC=0
        for i in range(SampleN):
            if i not in RemovedSet:
                SRS+=SampleRDs[i]
                SRC+=Stat.SampleReadCount[i]
        LastP=EstimatedP
        EstimatedP=SRS/SRC
    #the efficient way
    '''for i in range(WBegin,WEnd):
        SRS+=g.StatisticalRDWindowSums[i]
        SRC+=g.StatisticalReadCounts[i]
    SRC/=WEnd-WBegin'''
    #the efficient way 2
    '''for i in range(WBegin,WEnd):
        SRS+=g.RDWindowSums[i]
    SRC=g.AllReadCount'''
    return (SRS,SRC)"""

def scanConZero(SampleMRDRs):
    Begin=None
    AsZero=0.25
    Ints=[]
    Ave=0
    for i in range(len(SampleMRDRs)):
        R=SampleMRDRs[i]
        if R==0:
            if Begin==None:
                Begin=i
                Ave=0
            else:
                Ave+=R
        elif R<AsZero:
            Ave+=R
        else:
            if Begin==None:
                continue
            else:
                Ave/=i-Begin
                Ints.append(((Begin,i),Ave))
                Begin=None
    return Ints

def makeSampleRDEvidences(SampleData,SampleI,SampleName,Ploidy,RDWindowSize=None,SampleSegData=None,SampleDataAcc=None):
    print(gettime()+"segmenting %s..."%SampleName,file=sys.stderr)
    #SampleMRDRs=g.MyGenome.get(ContigName).MixedRDRs[SampleI]
    sys.stderr.flush()
    #for i in range(len(SampleMRDRs)):
    #    print(SampleMRDRs[i])
    #exit(0)
    if SampleSegData==None:
        SampleSegData=SampleData
    SampleIntervals=[]
    CutOffs=segmentation(SampleSegData)
    Last=0
    if SampleDataAcc==None:
        for End,Ave in CutOffs:
            Ave=0
            for i in range(Last,End):
                Ave+=SampleData[i]
            Ave/=End-Last
            SampleIntervals.append(RDInterval(SampleI,Last,End,Ave,Ploidy,None,RDWindowSize))
            Last=End
    else:
        for End,Ave in CutOffs:
            Ave=SampleDataAcc[End]-SampleDataAcc[Last]
            Ave/=End-Last
            SampleIntervals.append(RDInterval(SampleI,Last,End,Ave,Ploidy,None,RDWindowSize))
            Last=End
    SampleIntervals.sort(key=lambda I:I.WBegin)
    SampleEvidences=[]
    for I in SampleIntervals:
        if I.CN!=I.Ploidy:
            e=Evidence()
            e.setData(1,I)
            SampleEvidences.append(e)
    del SampleData
    del SampleSegData
    del SampleIntervals
    return SampleEvidences

def makeRDEvidences(Data,TheContig,SegData=None,DataAcc=None):#because robjects.r is singleton, use multiprocessing instead of multithreading #seems it's rpy2 that consumes much memory
    if SegData==None:
        SegData=Data
    if g.RParallel:
        CutOffs=dnacopy_cbs_multi(SegData)
        #Intervals=[]
        Evidences=[]
        for i in range(len(Data)):
            Ploidy=TheContig.Ploidies[i]
            Intervals=[]
            Evidences.append([])
            Last=0
            for End,Ave in CutOffs[i]:
                Ave=0
                for j in range(Last,End):
                    Ave+=Data[i][j]
                Ave/=End-Last
                Intervals.append(RDInterval(i,Last,End,Ave,Ploidy,None,g.RDWindowSize))
                Last=End
            Intervals.sort(key=lambda I:I.WBegin)
            for I in Intervals:
                if I.CN!=I.Ploidy:
                    e=Evidence()
                    e.setData(1,I)
                    Evidences[i].append(e)
        #TheContig.Intervals=Intervals
        return Evidences
    if g.ThreadN==1 or len(Data[0])<10000:#process cost is big
        Evidences=[None]*len(Data)
        for i in range(len(Data)):
            Evidences[i]=makeSampleRDEvidences(Data[i],i,g.SampleNames[i],TheContig.Ploidies[i],g.RDWindowSize,SegData[i],DataAcc[i] if DataAcc!=None else None)
    else:
        ctx=mp.get_context("fork")
        pool=ctx.Pool(g.ThreadN)
        args=[]
        for i in range(len(Data)):
            args.append((Data[i],i,g.SampleNames[i],TheContig.Ploidies[i],g.RDWindowSize,SegData[i],DataAcc[i] if DataAcc!=None else None))
        addPool(pool)
        Evidences=pool.starmap(makeSampleRDEvidences,args)
        print(gettime()+"Intervals for %s made. "%(TheContig.Name)+getMemUsage(),file=sys.stderr)
        delPool()
        pool.close()
    #TheContig.Intervals=Intervals
    return Evidences

def segmentation(data):
    #return SaRa(data)
    return dnacopy_cbs(data)

script=None
def dnacopy_cbs(data):
    datastring=""
    first=True
    for d in data:
        if not first:
            datastring+=","
        first=False
        #datastring+="%.7s"%(d)#math.log2(d/2.0 if d!=0 else sys.float_info.min))
        datastring+="%.7s"%(math.log2(d/2.0 if d!=0 else sys.float_info.min))
    robjects.r("rm(list = ls(all.names=TRUE))")
    robjects.r("rddata<-data.frame(mrd=c(%s))"%datastring)
    global script
    if script==None:
        sf=open("dnacopy_cbs.r","r")
        script=str(sf.read())
        sf.close()
    x=robjects.r(script)
    datavec=x[1]
    ends=datavec[3]
    means=datavec[5]
    CutOffs=[]
    for i in range(len(ends)):
        CutOffs.append((ends[i],means[i]))
    return CutOffs

def dnacopy_cbs_multi(data):
    robjects.r("rm(list = ls(all.names=TRUE))")
    robjects.r('rddatalist<-vector("list",%s)'%len(data))
    robjects.r("ThreadN<-%s"%g.ThreadN)
    print(gettime()+"Transfering data to R...",file=sys.stderr)
    for i in range(len(data)):
        s=data[i]
        datastring=""
        first=True
        for d in s:
            if not first:
                datastring+=","
            first=False
            #datastring+="%.7s"%(d)#math.log2(d/2.0 if d!=0 else sys.float_info.min))
            datastring+="%.7s"%(math.log2(d/2.0 if d!=0 else sys.float_info.min))
        robjects.r("rddatalist[[%s]]<-data.frame(mrd=c(%s))"%(i+1,datastring))
    global script
    if script==None:
        sf=open("dnacopy_cbs_multi.r","r")
        script=str(sf.read())
        sf.close()
    print(gettime()+"Segmenting...",file=sys.stderr)
    xs=robjects.r(script)
    print(gettime()+"Segmented! Making intervals...",file=sys.stderr)
    CutOffs=[]
    for x in xs:
        CutOffs.append([])
        datavec=x[1]
        ends=datavec[3]
        means=datavec[5]
        for i in range(len(ends)):
            CutOffs[-1].append((ends[i],means[i]))
    return CutOffs

def extractEvidences(Intervals):
    warn("Making evidences...")
    Evidences=[]
    for S in Intervals:
        for I in S:
            if I.CN!=I.Ploidy:
                e=Evidence()
                e.setData(1,I)
                Evidences.append(e)
    return Evidences

def makeRDICandidates(Evidences):
    warn("Making candidates...")
    Candidates=[]
    for s in Evidences:
        for e in s:
            Candidates.append(Candidate([e]))
    Candidates=combineCandidateSets(Candidates,[])
    for c in Candidates:
        c.unifyEvidences()
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

def getNormalRD_overall(TheContig,SampleI,WindowI):
    return 0 if g.StatisticalReadCounts[WindowI]==0 else g.StatisticalRDWindowSums[WindowI]*(g.SampleReadCount[SampleI]/g.StatisticalReadCounts[WindowI])

def getNormalRD_contig(TheContig,SampleI,WindowI):
    return 0 if g.StatisticalReadCounts[WindowI]==0 else g.StatisticalRDWindowSums[WindowI]*(TheContig.SampleReadCount[SampleI]/g.StatisticalReadCounts[WindowI])

if g.StatLocal:
    def getNormalRD(TheContig,SampleI,WindowI):
        return getNormalRD_contig(TheContig,SampleI,WindowI)
else:
    def getNormalRD(TheContig,SampleI,WindowI):
        return getNormalRD_overall(TheContig,SampleI,WindowI)

def getNormalRDDirect(StatisticalReadCounts, StatisticalRDWindowSums,SampleReadCount, SampleI, WindowI):
    return 0 if StatisticalReadCounts[WindowI]==0 else StatisticalRDWindowSums[WindowI]*(SampleReadCount[SampleI]/StatisticalReadCounts[WindowI])

def getIntervalNormalRD(TheContig,SampleI,WindowB,WindowE):
    return getIntervalNormalRD_overall(TheContig,SampleI,WindowB,WindowE)

def getIntervalNormalRD_overall(TheContig,SampleI,WindowB,WindowE):
    SP=getSP(TheContig,WindowB,WindowE)
    return 0 if SP[1]==0 else g.SampleReadCount[SampleI]/SP[1]*SP[0]

def getIntervalNormalRD_contig(TheContig,SampleI,WindowB,WindowE):
    SP=getSP(TheContig,WindowB,WindowE)
    return 0 if SP[1]==0 else TheContig.ContigSampleReadCounts[SampleI]/SP[1]*SP[0]

def analyzeRD(RDWindows,WindowsN,TheContig,NormalizationOnly=False):
    print(gettime()+"processing %s RD data..."%TheContig.Name,file=sys.stderr)
    #PPRDWindowAverages=[0]*WindowsN
    PP=0.9
    #RDWindowMedians=[0]*WindowsN
    RDWindowSums=array("d",[0]*WindowsN)
    g.RDWindowSums=RDWindowSums
    SampleN=len(RDWindows)
    #RDWindowStandards=PPRDWindowAverages
    #WindowsP=[0]*WindowsN#use maximal-likelyhood estimate to estimate poisson(np)'s p, with sequencing reads number as n. As a result, the p=(sum of rd of all samples)/(sum of n)
    #furthermore, the standard rd of sample i should be p*ni=(sum of rd)*((ni)/(sum of n))
    AllReadCount=0#sum of n
    StatisticalReadCounts=array("d",[0]*(WindowsN+1))#[0]*WindowsN#remove the least and last samples
    StatisticalRDWindowSums=array("d",[0]*(WindowsN+1))#[0]*WindowsN
    g.StatisticalRDWindowSums=StatisticalRDWindowSums
    g.StatisticalReadCounts=StatisticalReadCounts
    for i in range(SampleN):
        AllReadCount+=g.SampleReadCount[i]
    g.AllReadCount=AllReadCount

    MixedRDRs=[]
    for i in range(SampleN):
        MixedRDRs.append(array("f",[0]*WindowsN))#Mixed Read depth rate
    
    if g.StatLocal:
        SampleReadCount=TheContig.SampleReadCount
    else:
        SampleReadCount=g.SampleReadCount
    
    RDWindowsAcc=[]
    for i in range(SampleN):
        RDWindowsAcc.append(array("f",[0]*(WindowsN+1)))
    RDWindowStandardsAcc=array("f",[0]*(WindowsN+1))
    MixedRDRsAcc=[]
    for i in range(SampleN):
        MixedRDRsAcc.append(array("f",[0]*(WindowsN+1)))
    TheContig.MixedRDRs=MixedRDRs
    #TheContig.RDWindowStandards=RDWindowStandards
    TheContig.RDWindowsAcc=RDWindowsAcc
    TheContig.MixedRDRsAcc=MixedRDRsAcc

    processingRD(RDWindows, SampleN, WindowsN, MixedRDRs, RDWindowSums, RDWindowsAcc, MixedRDRsAcc, StatisticalReadCounts, StatisticalRDWindowSums, SampleReadCount, TheContig.Ploidies)
    
    g.ERD=1.0#ERD
    #g.MixedRDRs=MixedRDRs
    #TheContig.RDWindowStandardsAcc=RDWindowStandardsAcc
    #TheContig.MRMedians=MRMedians
    if NormalizationOnly:
        return MixedRDRs
    MRDEvidences=makeRDEvidences(MixedRDRs,TheContig,DataAcc=TheContig.MixedRDRsAcc)
    CLLEvidences=makeRDEvidences(MixedRDRs,TheContig,TheContig.AverageClipLengths,DataAcc=TheContig.MixedRDRsAcc)
    Evidences=combineEvidences(MRDEvidences,CLLEvidences,TheContig)
    EvidenceN=0
    for i in range(SampleN):
        EvidenceN+=len(Evidences[i])
    warn("Number of evidences:%s"%EvidenceN)
    Evidences=processEvidencesWithDRPs(Evidences,TheContig)
    Evidences=processEvidencesWithCLL(Evidences,TheContig)
    RDICandidates=makeRDICandidates(Evidences)
    SegSD=False
    if SegSD:
        SDCandidates=getSDCandidates(TheContig)
        RDICandidates=combineCandidateSets(RDICandidates,SDCandidates)
    return RDICandidates