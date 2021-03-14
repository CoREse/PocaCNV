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
        #if Ploidy==1:
            #for j in range(len(TheContig.RDWindows[i])):
            #    TheContig.RDWindows[i][j]*=2
            #TheContig.ContigSampleReadCounts[i]*=2

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
def mulikely(mu,mus):
    cd=poisson.cdf(mus,mu)
    return min(cd,1-cd)

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
    return (SRS,SRC)

class RDInterval:
    def __init__(self,Sample,WBegin,WEnd,ARD,Ploidy,TheContig=None,RDWindowSize=None):#multiprocessing will make g.RDWindowSize to default value
        self.WBegin=WBegin
        self.WEnd=WEnd
        self.AverageRD=ARD
        self.Sample=Sample
        self.Ploidy=Ploidy
        self.TheContig=TheContig
        if RDWindowSize==None:
            self.RDWindowSize=g.RDWindowSize
        else:
            self.RDWindowSize=RDWindowSize
        self.refresh()
    def refresh(self):
        self.Begin=self.WBegin*self.RDWindowSize
        self.End=self.WEnd*self.RDWindowSize
        self.SupportedSVType=None#0:del, 1:insertion, 2:dup
        self.mu=None
        self.mus=None
        if 1.8<self.AverageRD<2.2:
            self.CN=2
        elif 1.5<=self.AverageRD<=1.8:
            self.CN=1
        elif 2.2<=self.AverageRD<=1.5:
            self.CN=3
        else:
            self.CN=int(self.AverageRD+0.5)
        #if self.CN==0 or self.CN==1:
        if self.CN<self.Ploidy:
            self.SupportedSVType=0
        elif self.CN>self.Ploidy:
            self.SupportedSVType=2
    def setContig(TheContig):
        self.TheContig=TheContig
    def calcMuMus(self,TheContig=None,SP=None, local=g.StatLocal):#SP[0]=RDSum, SP[1]=Sample Read Count Sum(g.AllReadCount)
        if TheContig==None:
            TheContig=self.TheContig
        if TheContig==None:
            raise Exception("No contig given.")
        if SP==None:
            SP=getSP(TheContig,self.WBegin,self.WEnd,local)
        Stat=g
        if local:
            Stat=TheContig
        mu=SP[0]*(Stat.SampleReadCount[self.Sample]/SP[1])
        mus=TheContig.RDWindowsAcc[self.Sample][self.WEnd]-TheContig.RDWindowsAcc[self.Sample][self.WBegin]
        mu=int(mu+0.5)
        mus=int(mus+0.5)
        self.mu=mu
        self.mus=mus
        return mu,mus

    def calcMuMusOld(self,TheContig=None,RDS=None,RDSAcc=None):
        if TheContig==None:
            TheContig=self.TheContig
        if TheContig==None:
            raise Exception("No contig given.")
        RDSData=TheContig.RDWindowStandards
        if RDS!=None:
            RDSData=RDS
        length=self.WEnd-self.WBegin
        mu=0
        mus=0
        if RDSAcc==None:
            if RDS==None:
                RDSAcc=TheContig.RDWindowStandardsAcc
        mus=TheContig.RDWindowsAcc[self.Sample][self.WEnd]-TheContig.RDWindowsAcc[self.Sample][self.WBegin]
        if RDSAcc==None:
            for i in range(self.WBegin,self.WEnd):
                mu+=RDSData[i]
                #mus+=TheContig.RDWindows[self.Sample][i]
            #mus+=TheContig.MixedRDRs[self.Sample][i]/2.0*RDSData[i]
        else:
            mu=RDSAcc[self.WEnd]-RDSAcc[self.WBegin]
        mu=int(mu+0.5)
        mus=int(mus+0.5)
        self.mu=mu
        self.mus=mus
        return mu,mus

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

def makeSampleRDIntervals(SampleMRDRs,SampleI,SampleName,Ploidy,RDWindowSize=None):
    print(gettime()+"segmenting %s..."%SampleName,file=sys.stderr)
    sys.stderr.flush()
    SampleIntervals=[]
    CutOffs=segmentation(SampleMRDRs)
    Last=0
    for End,Ave in CutOffs:
        Ave=0
        for i in range(Last,End):
            Ave+=SampleMRDRs[i]
        Ave/=End-Last
        SampleIntervals.append(RDInterval(SampleI,Last,End,Ave,Ploidy,None,RDWindowSize))
        Last=End
    #CZInts=scanConZero(SampleMRDRs)
    #SampleIntervals.sort(key=lambda I:I.WBegin)
    #for I in CZInts:
    #    TempInt=RDInterval(SampleI,I[0][0],I[0][1],I[1],None,RDWindowSize)
    #    findNBindOrInsert(TempInt,SampleIntervals)
    SampleIntervals.sort(key=lambda I:I.WBegin)
    del SampleMRDRs
    return SampleIntervals

def findNBindOrInsert(Int,Ints):#Ints should be sorted
    BindPercentage=0.8
    Binded=False
    IntLength=Int.WEnd-Int.WBegin
    IntDis=IntLength*(1-BindPercentage)
    for I in Ints:
        Over=calcOverlap(Int.WBegin,Int.WEnd,I.WBegin,I.WEnd)
        if Over<0 and Int.WBegin-I.WBegin>IntDis:
            break
        if Over/(Int.WEnd-Int.WBegin)>=BindPercentage and Over/(I.WEnd-I.WBegin)>=BindPercentage:
            I.WBegin=int((Int.WBegin+I.WBegin)/2)
            Binded=True
    if not Binded:
        Ints.append(Int)

def makeRDIntervals(MixedRDRs,TheContig):#because robjects.r is singleton, use multiprocessing instead of multithreading #seems it's rpy2 that consumes much memory
    if g.ThreadN==1 or len(MixedRDRs[0])<10000:#process cost is big
        Intervals=[None]*len(MixedRDRs)
        for i in range(len(MixedRDRs)):
            Intervals[i]=makeSampleRDIntervals(MixedRDRs[i],i,g.SampleNames[i],TheContig.Ploidies[i],g.RDWindowSize)
    else:
        ctx=mp.get_context("spawn")
        pool=ctx.Pool(g.ThreadN)
        args=[]
        for i in range(len(MixedRDRs)):
            args.append((MixedRDRs[i],i,g.SampleNames[i],TheContig.Ploidies[i],g.RDWindowSize))
        addPool(pool)
        Intervals=pool.starmap(makeSampleRDIntervals,args)
        print(gettime()+"Intervals for %s made. "%(TheContig.Name)+getMemUsage(),file=sys.stderr)
        delPool()
        pool.terminate()
    TheContig.Intervals=Intervals
    return Intervals

def getSDCandidates(TheContig):
    SDData=TheContig.RDWindowStandards
    SDAve=statistics.mean(SDData)
    if SDAve==0:
        return []
    SDRs=[]
    SDAves=[SDAve]*len(SDData)
    for d in SDData:
        SDRs.append(d/SDAve)
    CutOffs=segmentation(SDRs)
    TheContig.SDSegments=[]
    Last=0
    SDIntervals=[]
    for End,Ave in CutOffs:
        Ave=0
        for i in range(Last,End):
            Ave+=SDRs[i]
        Ave/=End-Last
        TheContig.SDSegments.append((End,Ave))
        if Ave==0 or 0.6<Ave<1.4:
            continue
        for s in range(len(g.SampleNames)):
            SAve=0
            for i in range(Last,End):
                SAve+=TheContig.MixedRDRs[s][i]
            SAve/=End-Last
            SInt=RDInterval(s,Last,End,SAve,TheContig.Ploidies[s])
            if cn2likely(SInt,TheContig,SDAves)<1-g.Parameters.CN2FilterConfidence:
                SDIntervals.append(SInt)
        Last=End
    Evidences=[]
    for I in SDIntervals:
        e=Evidence()
        e.setData(1,I)
        Evidences.append(e)
    SampleCandidates=makeRDICandidates(Evidences)
    return SampleCandidates

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

def extractEvidences(Intervals):
    Evidences=[]
    for S in Intervals:
        for I in S:
            if I.CN!=I.Ploidy:
                e=Evidence()
                e.setData(1,I)
                Evidences.append(e)
    return Evidences

def getExtendables(I,OI):#return None or list with at most 2 (begin, end)
    Over=calcOverlap(I.WBegin,I.WEnd,OI.WBegin,OI.WEnd)
    ILength=I.WEnd-I.WBegin
    OLength=OI.WEnd-OI.WBegin
    Ex=[]
    if Over==-1:
        Endure=(ILength)*0.05
        if Endure<1:
            Endure=1
        if Endure>10:
            Endure=10
        Endure=int(Endure)
        if I.WEnd<OI.WBegin and OI.WBegin-I.WEnd<=Endure:
            return [(I.WEnd,OI.WEnd)]
        if OI.WEnd<I.WBegin and I.WBegin-OI.WEnd<=Endure:
            return [(OI.WBegin,I.WBegin)]
        return None
    if Over>=OLength:
        return None
    if Over>=ILength:
        Edges=[I.WBegin,I.WEnd,OI.WBegin,OI.WEnd]
        Edges.sort()
        if Edges[1]>Edges[0]:
            Ex.append((Edges[0],Edges[1]))
        if Edges[3]>Edges[2]:
            Ex.append((Edges[2],Edges[3]))
        if len(Ex)==0:
            return None
        return Ex
    if I.WEnd<OI.WEnd:
        return [(I.WEnd,OI.WEnd)]
    return [(OI.WBegin,I.WBegin)]
    return None

def verifyExtendables(I,Extendables,TheContig):
    if Extendables==None or len(Extendables)==0:
        return False
    if I.CN==2:
        return False
    if I.mus==None:
        mu,mus=I.calcMuMus(TheContig)
    else:
        mu,mus=(I.mu,I.mus)
    muo=0
    for e in Extendables:
        for i in range(e[0],e[1]):
            muo+=TheContig.MixedRDRs[I.Sample][i]/2.0*TheContig.RDWindowStandards[i]
    if mulikely(mu,muo)<mulikely(mu,mus)*1.3 and ((mu<mus and mu<muo) or(mu>mus and mu>muo)):
        return True
    return False

def extend(I,Extendables):
    if Extendables==None or len(Extendables)==0:
        return I.WBegin,I.WEnd
    if len(Extendables)==1:
        if Extendables[0][0]==I.WEnd:
            return I.WBegin,Extendables[0][1]
        return Extendables[0][0],I.WEnd
    return Extendables[0][0],Extendables[1][1]

def extendIntervals(Intervals,TheContig):#extend intervals by examine other samples intervals
    for s in range(len(Intervals)):
        for i in range(len(Intervals[s])):
            I=Intervals[s][i]
            ExtendBegin=None
            ExtendEnd=None
            for SI in Intervals:
                for OI in SI:
                    Extendables=getExtendables(I,OI)
                    if verifyExtendables(I,Extendables,TheContig):
                        ExBegin,ExEnd=extend(I,Extendables)
                        if ExtendBegin==None:
                            ExtendBegin=ExBegin
                        else:
                            ExtendBegin=min(ExtendBegin,ExBegin)
                        if ExtendEnd==None:
                            ExtendEnd=ExEnd
                        else:
                            ExtendEnd=max(ExtendEnd,ExEnd)
            if ExtendBegin!=None and ExtendEnd!=None:
                Intervals[s][i].WBegin=ExtendBegin
                Intervals[s][i].WEnd=ExtendEnd
                Ave=0
                for j in range(ExtendBegin,ExtendEnd):
                    Ave+=TheContig.MixedRDRs[s][j]
                Ave/=ExtendEnd-ExtendBegin
                Intervals[s][i].AverageRD=Ave
                Intervals[s][i].refresh()
    return Intervals

def extendEvidences(Evidences,TheContig):#extend intervals by examine other samples intervals
    return Evidences#obsolete
    print(gettime()+"Extending evidences(%s)..."%(len(Evidences)),file=sys.stderr)
    Evidences.sort(key=lambda e:e.Data.WBegin)
    S=0#where now.WBegin > all WEnd before s
    for s in range(len(Evidences)):
        I=Evidences[s].Data
        for i in range(S,s):
            if I.WBegin>Evidences[i].Data.WBegin+10:
                S=i+1
            else:
                break
        ExtendBegin=None
        ExtendEnd=None
        MayBreak=False
        if (s+1)%10000==0:
            print(gettime()+"%s evidences extended."%(s+1),file=sys.stderr)
        for s in range(S,len(Evidences)):
            E=Evidences[s]
            OI=E.Data
            if OI.WBegin>I.WEnd:
                MayBreak=True
            Extendables=getExtendables(I,OI)
            if Extendables==None and MayBreak:
                break
            if verifyExtendables(I,Extendables,TheContig):
                ExBegin,ExEnd=extend(I,Extendables)
                if ExtendBegin==None:
                    ExtendBegin=ExBegin
                else:
                    ExtendBegin=min(ExtendBegin,ExBegin)
                if ExtendEnd==None:
                    ExtendEnd=ExEnd
                else:
                    ExtendEnd=max(ExtendEnd,ExEnd)
        if ExtendBegin!=None and ExtendEnd!=None:
            Evidences[s].Data.WBegin=ExtendBegin
            Evidences[s].Data.WEnd=ExtendEnd
            Ave=0
            for j in range(ExtendBegin,ExtendEnd):
                Ave+=TheContig.MixedRDRs[I.Sample][j]
            Ave/=ExtendEnd-ExtendBegin
            Evidences[s].Data.AverageRD=Ave
            Evidences[s].Data.refresh()
            Evidences[s].analyzeData()
    return Evidences

def makeRDICandidates(Evidences):
    Candidates=[]
    for e in Evidences:
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

def unifyRD(RDWindows):
    for s in range(len(RDWindows)):
        for i in range(len(RDWindows[s])):
            RDWindows[s][i]/=g.SequenceDepthRatio[s]

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

def getIntervalNormalRDOld(TheContig, SampleI, WindowB, WindowE,PP=0.9):#for a window interval
    SampleN=len(TheContig.RDWindows)
    WindowN=len(TheContig.RDWindows[0])
    SampleSums=[]
    for i in range(SampleN):
        Sum=0
        for j in range(WindowN):
            Sum+=TheContig.RDWindows[i][j]
        SampleSums.append((Sum,i))
    SampleSums.sort(key=lambda s:s[0])
    Sum=0
    De=int(SampleN*0.5*(1-PP))
    StatReadCount=0
    for s in SampleSums[De:SampleN-De]:
        Sum+=s[0]
        StatReadCount+=g.SampleReadCount[s[1]]
    return Sum*(g.SampleReadCount[SampleI]/StatReadCount)

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
    RDWindowAverages=[0]*WindowsN
    PPRDWindowAverages=[0]*WindowsN
    PP=0.9
    Mean2Adjust=False
    Smooth=1
    #RDWindowMedians=[0]*WindowsN
    RDWindowSums=[0]*WindowsN
    g.RDWindowSums=RDWindowSums
    SampleN=len(RDWindows)
    SampleSums=[0]*SampleN
    SampleAverages=[0]*SampleN
    RDWindowStandards=PPRDWindowAverages
    #WindowsP=[0]*WindowsN#use maximal-likelyhood estimate to estimate poisson(np)'s p, with sequencing reads number as n. As a result, the p=(sum of rd of all samples)/(sum of n)
    #furthermore, the standard rd of sample i should be p*ni=(sum of rd)*((ni)/(sum of n))
    AllReadCount=0#sum of n
    StatisticalReadCounts=[0]*WindowsN#remove the least and last samples
    StatisticalRDWindowSums=[0]*WindowsN
    for i in range(SampleN):
        AllReadCount+=g.SampleReadCount[i]
    g.AllReadCount=AllReadCount
    #for i in range(WindowsN):
    #    WindowsP[i]=RDWindowSums[i]/AllReadCount
    #use estimate p to ultilize the data from different sequencing coverage
    #unifyRD(RDWindows)#rd/ratio doesnt' obey same dist that coverage/ratio obeys(different variance), so it's better all samples from the same sequence coverage
    if SampleN<2:#WR,SR组合不太科学，假如几个样本同时有某个变异，那么很可能无法检测出这个变异，样本很多时不如直接用WR（基于变异占少数的假设），但假如变异本身不罕见就又有问题了
        for i in range(WindowsN):
            for j in range(SampleN):
                #RDWindowSums[i]+=RDWindows[j][i]
                SampleSums[j]+=RDWindows[j][i]
            #RDWindowAverages[i]=RDWindowSums[i]/SampleN
        MixedRDRs=[]
        for j in range(SampleN):
            SampleAverages[j]=SampleSums[j]/WindowsN
            MixedRDRs.append(array("f",[0]*WindowsN))
        #MixedRDRs=[[0]*WindowsN]*SampleN#Mixed Read depth rate, THIS is AWFUL! this will make SampleN copy of array objects, sharing the same memory space!
        """
        ZCount=0
        for i in range(OccurredWindowsN):
            if RDWindows[0][int(RefStartPos[RefInd["chr22"]]/100)+i]==0:
                ZCount+=1
        print(OccurredWindowsN, ZCount, SampleSums[0],SampleAverages[0])
        """
        for i in range(SampleN):
            for j in range(WindowsN):
                SR=(RDWindows[i][j]/SampleAverages[i]) if SampleAverages[i]!=0 else 0
                MixedRDRs[i][j]=SR
                #WR=(RDWindows[i][j]/RDWindowAverages[j]) if RDWindowAverages[j]!=0 else 0
                #MixedRDRs[i][j]=(SR-WR)/SampleN+WR
    else:
        WindowSamples=[]
        for j in range(SampleN):
            WindowSamples.append((0,0))
        for i in range(WindowsN):
            for j in range(SampleN):
                RDWindowSums[i]+=RDWindows[j][i]
                SampleSums[j]+=RDWindows[j][i]
                #WindowSamples[j]=(RDWindows[j][i],j)
            RDWindowAverages[i]=RDWindowSums[i]/SampleN
            #RDWindowMedians[i]=statistics.median(WindowSamples)
            #WindowSamples.sort(key=lambda s:s[0])
            #De=int(SampleN*0.5*(1-PP))
            SP=getSP(TheContig,i,i+1)
            StatisticalRDWindowSums[i]=SP[0]
            StatisticalReadCounts[i]=SP[1]
            #for s in WindowSamples[De:SampleN-De]:
            #    StatisticalRDWindowSums[i]+=s[0]
            #    StatisticalReadCounts[i]+=g.SampleReadCount[s[1]]
            #PPRDWindowAverages[i]=statistics.mean(WindowSamples[De:SampleN-De])
        g.StatisticalRDWindowSums=StatisticalRDWindowSums
        g.StatisticalReadCounts=StatisticalReadCounts
        AllZeroLeft=-1
        AllZeroRight=WindowsN
        Left1=False
        Right1=False
        for i in range(WindowsN):
            if not Right1:
                j=WindowsN-i
                if j>=0:
                    if RDWindowSums[i]!=0:
                        Right1=True
                    else:
                        AllZeroRight=j
            if not Left1:
                if RDWindowSums[i]!=0:
                    Left1=True
                else:
                    AllZeroLeft=i
            if Right1 and Left1:
                break
        SampleSumAverage=0
        for j in range(SampleN):
            SampleAverages[j]=SampleSums[j]/WindowsN
            SampleSumAverage+=SampleSums[j]
        SampleSumAverage/=SampleN
        #SampleSequenceDepthRatio=[0]*SampleN
        #SampleMedians=[0]*SampleN
        #print(SampleSumAverage, SampleSums, SampleMedians, TheContig.Name)
        #for j in range(SampleN):
        #    SampleSequenceDepthRatio[j]=SampleSums[j]/SampleSumAverage if SampleSumAverage!=0 else 1//already unified
            #SampleMedians[j]=statistics.median(RDWindows[j])
        
        MixedRDRs=[]
        for i in range(SampleN):
            MixedRDRs.append(array("f",[0]*WindowsN))#Mixed Read depth rate
        for i in range(SampleN):
            for j in range(WindowsN):
                S0Value=0
                if j<=AllZeroLeft or j>=AllZeroRight:
                    S0Value=1
                #WR=(RDWindows[i][j]/RDWindowStandards[j]) if RDWindowStandards[j]!=0 else S0Value
                Standard=getNormalRD(TheContig,i,j)
                WR=RDWindows[i][j]/Standard if Standard!=0 else S0Value
                if Standard==0 and RDWindows[i][j]!=0:
                    WR=RDWindows[i][j]/RDWindowAverages[j]
                MixedRDRs[i][j]=WR
                #if j==10000:
                #    print("SRC:%s,SRS:%s,SR:%s,Standard:%s,WR:%s"%(g.StatisticalReadCounts[j],g.StatisticalRDWindowSums[j],g.SampleReadCount[i],Standard,WR),file=sys.stderr)
                #MixedRDRs[i][j]=(WR/SampleSequenceDepthRatio[i]) if SampleSequenceDepthRatio[i]!=0 else 0
                '''
                SR=(RDWindows[i][j]/SampleAverages[i]) if SampleAverages[i]!=0 else 0
                MixedRDRs[i][j]=(SR-WR)/SampleN+WR
                '''
        #MRMedians=[0]*SampleN
        #for j in range(SampleN):
        #    MRMedians[j]=statistics.median(MixedRDRs[j])
        for i in range(SampleN):
            for j in range(WindowsN):
                #MixedRDRs[i][j]=((MixedRDRs[i][j]/MRMedians[i])*2.0) if MRMedians[i]!=0 else 0#standardization to make median 2.0#Nonesense
                MixedRDRs[i][j]*=2.0#diploid
                #if RDWindowStandards[j]==0:
                    #MixedRDRs[i][j]=2#if windows average is 0, we consider here is not valuable, so mark as normal(CN=2)
        if Smooth>1:
            if Smooth%2==0:
                Smooth+=1
            for i in range(SampleN):
                Smoother=array("f",[0]*WindowsN)
                for j in range(WindowsN):
                    SH=int(Smooth/2)
                    SB=0 if j<SH else j-SH
                    SE=WindowsN if j+SH>WindowsN else j+SH
                    for k in range(SB,SE):
                        Smoother[j]+=MixedRDRs[i][k]
                    Smoother[j]/=Smooth
                MixedRDRs[i]=Smoother
        if Mean2Adjust:
            #sf=open("data/SampleCountMeanData.txt","w")
            #first=True
            for i in range(SampleN):
                #if not first:
                #    print("\n",end="",file=sf)
                #first=False
                Mean=0
                for j in range(WindowsN):
                    Mean+=MixedRDRs[i][j]
                Mean/=WindowsN
                ARatio=2.0/Mean
                for j in range(WindowsN):
                    MixedRDRs[i][j]*=ARatio
                #print("%s: Count:%s, Mean:%s"%(TheContig.SampleNames[i],g.SampleReadCount[i],Mean),end="",file=sf)
            #sf.close()
            #exit(0)
    """
    rdtestfile=open("data/rdtest.txt","a")
    for i in range(WindowsN):
        print("%d %f"%(i,MixedRDRs[0][i]),file=rdtestfile)
    rdtestfile.close()
    exit(0)
    """
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
    '''sf=open("data/SampleCountMeanDataContigCount.txt","w")
    first=True
    for i in range(SampleN):
        if not first:
            print("\n",end="",file=sf)
        first=False
        Mean=0
        for j in range(WindowsN):
            Mean+=MixedRDRs[i][j]
        Mean/=WindowsN
        print("%s: Count:%s, Mean:%s, ContigCount:%s"%(TheContig.SampleNames[i],g.SampleReadCount[i],Mean,TheContig.SampleReadCount[i]),end="",file=sf)
    sf.close()
    exit(0)'''
    #sf=open("data/HG00404_w100_standard.mrd","w")
    #first=True
    #for j in range(WindowsN):
    #    if not first:
    #        print("\n",end="",file=sf)
    #    first=False
    #    print("%s,%s"%(MixedRDRs[0][j],getNormalRD(TheContig,0,j)),end="",file=sf)
    #sf.close()
    #exit(0)
    RDWindowStandardsAcc=array("f",[0]*(WindowsN+1))
    Sum=0
    for i in range(WindowsN):
        RDWindowStandardsAcc[i]=Sum
        Sum+=RDWindowStandards[i]
    RDWindowStandardsAcc[WindowsN]=Sum
    RDWindowsAcc=[]
    for i in range(SampleN):
        RDWindowsAcc.append(array("f",[0]*(WindowsN+1)))
        Sum=0
        for j in range(WindowsN):
            RDWindowsAcc[i][j]=Sum
            Sum+=RDWindows[i][j]
        RDWindowsAcc[i][WindowsN]=Sum
    g.ERD=1.0#ERD
    #g.MixedRDRs=MixedRDRs
    TheContig.MixedRDRs=MixedRDRs
    TheContig.RDWindowStandards=RDWindowStandards
    TheContig.RDWindowsAcc=RDWindowsAcc
    TheContig.RDWindowStandardsAcc=RDWindowStandardsAcc
    #TheContig.MRMedians=MRMedians
    if NormalizationOnly:
        return MixedRDRs
    
    #partition(MixedRDRs)

    RDICandidates=makeRDICandidates(extractEvidences(makeRDIntervals(MixedRDRs,TheContig)))
    SegSD=False
    if SegSD:
        SDCandidates=getSDCandidates(TheContig)
        RDICandidates=combineCandidateSets(RDICandidates,SDCandidates)
    return RDICandidates

import pysam
def writeRDData(mygenome,SampleNames,SampleReadCount):
    for i in range(len(SampleNames)):
        writeSampleRDData(mygenome,SampleNames[i],i,SampleReadCount[i],g.RDWindowSize,g.SamplePaths[i])
    return

def getRDFPath(SamPath,Path="data"):
    SamFileName=SamPath.split("/")[-1].split("\\")[-1]
    return "%s/%s.rdf"%(Path,SamFileName)

def writeSampleRDData(mygenome, SampleName, SampleI, SampleReadCount, WindowSize, SamplePath, Path="data"):
    #SamFileName=g.SamplePaths[SampleI].split("/")[-1].split("\\")[-1]
    rdfile=open(getRDFPath(SamplePath,Path),"w")
    print("##Sample:%s"%SampleName,end="",file=rdfile)
    print("\n##SampleReadCount:%s\n##WindowSize:%s"%(SampleReadCount,WindowSize),end="",file=rdfile)
    for c in mygenome.Contigs:
        print("\n#%s %s"%(c.Name,c.Length),end="",file=rdfile)
        for j in range(len(c.RDWindows[SampleI])):
            print("\n%d"%(c.RDWindows[SampleI][j]),end="",file=rdfile)
    rdfile.close()
    return

def writeMixedRDData(mygenome,SampleNames):
    #ReferenceFile:pysam.FastaFile
    for i in range(len(SampleNames)):
        SamFileName=g.SamplePaths[i].split("/")[-1].split("\\")[-1]
        rdfile=open("data/%s.mrd"%(SamFileName),"w")
        print("##Sample:%s"%SampleNames[i],end="",file=rdfile)
        for c in mygenome.Contigs:
            print("\n#%s %s"%(c.Name,c.Length),end="",file=rdfile)
            for j in range(len(c.MixedRDRs[i])):
                print("\n%.8s"%(c.MixedRDRs[i][j]),end="",file=rdfile)
        rdfile.close()
    return

def readRDData(mygenome, SampleNames, FileName):
    SampleName=FileName.split("\\")[-1].split("/")[-1][:-4]
    if "rd"==SampleName[:2]:
        SampleNameS=SampleName.split("rd")
        SampleName=""
        for s in SampleNameS[1:]:
            SampleName+=s
    DataFile=open(FileName,"r")
    ContigName=None
    mygenome.addSample(SampleName)
    Skip=False
    Windows=None
    ConI=0
    ReadCount=0
    SampleReadCount=None
    for line in DataFile:
        if line[0]=='#':
            if line[1]=="#":
                sl=line[2:].split(":",maxsplit=1)
                if sl[0]=="Sample":
                    SampleName=sl[1].strip()
                elif sl[0]=="SampleReadCount":
                    SampleReadCount=int(sl[1].strip())
            else:
                sl=line[1:].split()
                ContigName=sl[0]
                if not mygenome.hasContig(ContigName):
                    Skip=True
                    continue
                Length=int(sl[1])
                TheContig=mygenome.get(ContigName)
                Windows=TheContig.RDWindows[-1]
                ContigSampleReadCounts=TheContig.ContigSampleReadCounts
                Skip=False
                ConI=0
            continue
        if Skip:
            continue
        #try:
        #    Windows[ConI]=int(line)
        #except IndexError as e:
        #    if len(Windows)<=ConI:
        #        print("data exceed contig %s's capacity(data no.%s, len of %s:%s)!"%(ContigName,ConI,ContigName,len(Windows)),file=sys.stderr)
        #    else:
        #        raise e
        Windows[ConI]=float(line)
        ReadCount+=Windows[ConI]
        ContigSampleReadCounts[-1]+=Windows[ConI]
        ConI+=1
    SampleNames.append(SampleName)
    mygenome.changeSampleName(-1,SampleName)
    if SampleReadCount==None:
        g.SampleReadCount.append(ReadCount)
    else:
        g.SampleReadCount.append(SampleReadCount)
    DataFile.close()
    return