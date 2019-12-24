import globals as g
from utils import *
import sys
from contig import *
import statistics
from array import array
import rpy2.robjects as robjects
from scipy.stats import poisson

def cn2filter(Interval,TheContig,Confidence=None):
    if Interval.mu==None:
        mu,mus=Interval.calcMuMus(TheContig)
    else:
        mu,mus=(Interval.mu,Interval.mus)
    if Confidence==None:
        Confidence=g.CN2FilterConfidence
    qint=poisson.interval(Confidence,mu)
    if qint[0]<mus<qint[1]:
        return False
    return True

class RDInterval:
    def __init__(self,Sample,WBegin,WEnd,ARD,TheContig=None):
        self.WBegin=WBegin
        self.WEnd=WEnd
        self.AverageRD=ARD
        self.Begin=WBegin*g.RDWindowSize
        self.End=WEnd*g.RDWindowSize
        self.Sample=Sample
        self.SupportedSVType=None#0:del, 1:insertion, 2:dup
        self.mu=None
        self.mus=None
        self.TheContig=TheContig
        if 1.8<self.AverageRD<2.2:
            self.CN=2
        elif 1.5<=self.AverageRD<=1.8:
            self.CN=1
        elif 2.2<=self.AverageRD<=1.5:
            self.CN=3
        else:
            self.CN=int(self.AverageRD+0.5)
        if self.CN==0 or self.CN==1:
            self.SupportedSVType=0
        elif self.CN>2:
            self.SupportedSVType=2
    def setContig(TheContig):
        self.TheContig=TheContig
    def calcMuMus(self,TheContig=None):
        if TheContig==None:
            TheContig=self.TheContig
        if TheContig==None:
            raise Exception("No contig given.")
        length=self.WEnd-self.WBegin
        mu=0
        mus=0
        for i in range(self.WBegin,self.WEnd):
            mu+=TheContig.RDWindowStandards[i]
            mus+=TheContig.MixedRDRs[self.Sample][i]/2.0*TheContig.RDWindowStandards[i]
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

def makeSampleRDIntervals(SampleMRDRs,SampleI,SampleIntervals,SampleName):
    print(gettime()+"segmenting %s..."%SampleName,file=sys.stderr)
    CutOffs=segmentation(SampleMRDRs)
    Last=0
    for End,Ave in CutOffs:
        SampleIntervals.append(RDInterval(SampleI,Last,End,Ave))
        Last=End

def makeRDIntervals(MixedRDRs):#because robjects.r is singleton, use multiprocessing instead of multithreading
    Intervals=[None]*len(MixedRDRs)
    if g.ThreadN==1 or len(MixedRDRs[0])<10000:#process cost is big
        for i in range(len(MixedRDRs)):
            Intervals[i]=[]
            makeSampleRDIntervals(MixedRDRs[i],i,Intervals[i])
    else:
        manager=g.Manager
        pool=g.Pool
        args=[]
        #MMixedRDRs=manager.list(MixedRDRs)
        for i in range(len(MixedRDRs)):
            Intervals[i]=manager.list()
            args.append((MixedRDRs[i],i,Intervals[i],g.SampleNames[i]))
        pool.starmap(makeSampleRDIntervals,args)
    return Intervals

def segmentation(data):
    return dnacopy_cbs(data)

script=None
def dnacopy_cbs(data):
    datastring=""
    first=True
    for d in data:
        if not first:
            datastring+=","
        first=False
        datastring+="%.7s"%d
    robjects.r("rddata=data.frame(mrd=c(%s))"%datastring)
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
            if I.CN!=2:
                e=Evidence()
                e.setData(1,I)
                Evidences.append(e)
    return Evidences
            

def makeRDICandidates(Evidences):
    Candidates=[]
    for e in Evidences:
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

def unifyRD(RDWindows):
    for s in range(len(RDWindows)):
        for i in range(len(RDWindows[s])):
            RDWindows[s][i]/=g.SequenceDepthRatio[s]

def analyzeRD(RDWindows,WindowsN,TheContig,NormalizationOnly=False):
    print(gettime()+"processing %s RD data..."%TheContig.Name,file=sys.stderr)
    RDWindowAverages=[0]*WindowsN
    RDWindowMedians=[0]*WindowsN
    RDWindowSums=[0]*WindowsN
    SampleN=len(RDWindows)
    SampleSums=[0]*SampleN
    SampleAverages=[0]*SampleN
    RDWindowStandards=RDWindowMedians
    unifyRD(RDWindows)#rd/ratio doesnt' obey same dist that coverage/ratio obeys(different variance), so it's better all samples from the same sequence coverage
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
        WindowSamples=[0]*SampleN
        for i in range(WindowsN):
            for j in range(SampleN):
                RDWindowSums[i]+=RDWindows[j][i]
                SampleSums[j]+=RDWindows[j][i]
                WindowSamples[j]=RDWindows[j][i]
            RDWindowAverages[i]=RDWindowSums[i]/SampleN
            RDWindowMedians[i]=statistics.median(WindowSamples)
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
                WR=(RDWindows[i][j]/RDWindowStandards[j]) if RDWindowStandards[j]!=0 else 0
                MixedRDRs[i][j]=WR
                #MixedRDRs[i][j]=(WR/SampleSequenceDepthRatio[i]) if SampleSequenceDepthRatio[i]!=0 else 0
                '''
                SR=(RDWindows[i][j]/SampleAverages[i]) if SampleAverages[i]!=0 else 0
                MixedRDRs[i][j]=(SR-WR)/SampleN+WR
                '''
        MRMedians=[0]*SampleN
        for j in range(SampleN):
            MRMedians[j]=statistics.median(MixedRDRs[j])
        for i in range(SampleN):
            for j in range(WindowsN):
                #MixedRDRs[i][j]=((MixedRDRs[i][j]/MRMedians[i])*2.0) if MRMedians[i]!=0 else 0#standardization to make median 2.0#Nonesense
                MixedRDRs[i][j]*=2.0#diploid
                if RDWindowStandards[j]==0:
                    MixedRDRs[i][j]=2#if windows average is 0, we consider here is not valuable, so mark as normal(CN=2)
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
    g.ERD=1.0#ERD
    g.MixedRDRs=MixedRDRs
    TheContig.MixedRDRs=MixedRDRs
    TheContig.RDWindowStandards=RDWindowStandards
    TheContig.MRMedians=MRMedians
    if NormalizationOnly:
        return MixedRDRs
        
    #partition(MixedRDRs)

    RDICandidates=makeRDICandidates(extractEvidences(makeRDIntervals(MixedRDRs)))
    return RDICandidates

import pysam
def writeRDData(mygenome,SampleNames):
    for i in range(len(SampleNames)):
        writeSampleRDData(mygenome,SampleNames[i],i)
    return

def writeSampleRDData(mygenome, SampleName, SampleI):
    rdfile=open("data/rd%s.rdf"%(SampleName),"w")
    first=True
    for c in mygenome.Contigs:
        if not first:
            print("\n",end="",file=rdfile)
        first=False
        print("#%s %s"%(c.Name,c.Length),end="",file=rdfile)
        for j in range(len(c.RDWindows[SampleI])):
            print("\n%d"%(c.RDWindows[SampleI][j]),end="",file=rdfile)
    rdfile.close()
    return

def writeMixedRDData(mygenome,ReferenceFile,SampleNames):
    ReferenceFile:pysam.FastaFile
    for i in range(len(SampleNames)):
        rdfile=open("data/mrd%s.mrd"%(SampleNames[i]),"w")
        first=True
        for c in mygenome.Contigs:
            if not first:
                print("\n",end="",file=rdfile)
            first=False
            print("#%s %s"%(c.Name,c.Length),end="",file=rdfile)
            for j in range(len(c.MixedRDRs[i])):
                print("\n%.8s"%(c.MixedRDRs[i][j]),end="",file=rdfile)
        rdfile.close()
    return

def readRDData(mygenome, SampleNames, FileName):
    SampleName=FileName.split("\\")[-1].split("/")[-1][:-4]
    SampleNameS=SampleName.split("rd")
    SampleName=""
    for s in SampleNameS[1:]:
        SampleName+=s
    SampleNames.append(SampleName)
    DataFile=open(FileName,"r")
    ContigName=None
    mygenome.addSample(SampleName)
    Skip=False
    Windows=None
    ConI=0
    ReadCount=0
    for line in DataFile:
        if line[0]=='#':
            sl=line[1:].split()
            ContigName=sl[0]
            if not mygenome.hasContig(ContigName):
                Skip=True
                continue
            Length=int(sl[1])
            Windows=mygenome.get(ContigName).RDWindows[-1]
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
        ConI+=1
    g.SampleReadCount.append(ReadCount)
    DataFile.close()
    return