import globals as g
import time
import os
import psutil

def gettime():
    return time.strftime("[%Y.%m.%d,%H:%M:%S]",time.localtime())

def getMemUsage():
    vms=0
    rss=0
    for p in g.Processes:
        rss+=p.memory_info().rss
        vms+=p.memory_info().vms
    return "Memroy usage:%.5sgb(rss),%.5sgb(vms)."%(rss/1024/1024/1024,vms/1024/1024/1024)

def getPid(i):
    return os.getpid()

def addPool(ThePool):
    Pids=ThePool.map(getPid,range(g.ThreadN))
    for pid in Pids:
        g.Processes.append(psutil.Process(pid))

def delPool():
    g.Processes=g.Processes[:1]

def calclulateSequenceDepthRatio(SequenceDepths=None):
    if SequenceDepths==None:
        data=g.SampleReadCount
    else:
        data=SequenceDepths
    Ave=0
    for c in data:
        Ave+=c
    Ave/=len(data)
    if Ave!=0:
        for c in data:
            g.SequenceDepthRatio.append(c/Ave)
    else:
        for c in data:
            g.SequenceDepthRatio.append(1)

def calcOverlap(B1,E1,B2,E2):
    if B1<=E2 and B1>=B2:
        Overlap=min(E2-B1,E1-B1)
    elif E1>=B2 and E1<=E2:
        Overlap=E1-B2
    elif B2>=B1 and B2<=E1:#1 covers 2
        Overlap=E2-B2
    else:
        Overlap=-1
    return Overlap

def inclusion(In1,In2):#return 0: not overlapped, 1: overlapped, 2: In1 include In2, 3: In2 inlcude In1, 4: identical
    if In1==In2:
        return 4
    elif In1[0]>=In2[0] and In1[1]<=In2[1]:
        return 3
    elif In1[0]<=In2[0] and In1[1]>=In2[1]:
        return 2
    elif In1[0]>=In2[0] and In1[0]<=In2[1] or In1[1]>=In2[0] and In1[1]<=In2[1]:
        return 1
    return 0

def getTidByCord(Cordinate):
    i=1
    while i<len(g.RefStartPos):
        if Cordinate<g.RefStartPos[i]:
            return i-1
        i+=1
    return i-1

#Candidate can have evidences from different sapmles, evidence will only be one sample's evidence
class Candidate:
    CombinePercentage=0.97
    CombineRange=300
    RangePercentage=0.03
    MinRange=100
    MaxRange=10000
    CombineMode=2#0: overlap > CombinePercentage of both Candidate, 1: each evidence has begin and end range within CombineRange with each other, 2: like one, but with flexible CombineRange due to length, 3: like one, but will average the begin and end(by the evidence)
    AveRange=True
    AveBreak=False
    def __init__(self,Es=[]):
        self.SVType=None#0: deletion, 1: insertion, 2:dup, 3:CNV
        self.Evidences=Es
        self.Begin=0xffffffff
        self.End=0
        self.BreakLeft=0
        self.BreakRight=0xffffffff
        self.BeginRange=(0,0)
        self.EndRange=(0,0)
        self.deductSVType()
        self.calculateSpread()
        self.CombineRange=None
    
    def calculateSpread(self):
        SumLeft=0
        SumRight=0
        BN=0
        if Candidate.CombineMode==3:
            SBegin=0
            SEnd=0
            for e in self.Evidences:
                SBegin+=e.Begin
                SEnd+=e.End
            self.Begin=int(SBegin/len(self.Evidences))
            self.End=int(SEnd/len(self.Evidences))
            self.BreakLeft=self.Begin
            self.BreakRight=self.End
        else:
            for e in self.Evidences:
                self.Begin, self.End=(min(self.Begin,e.Begin),max(self.End,e.End))
                if e.Type==0:
                    for pi in e.Data:
                        if Candidate.AveBreak:
                            SumLeft+=pi.LEnd
                            SumRight+=pi.RStart
                            BN+=1
                        else:
                            self.BreakLeft=max(self.BreakLeft,pi.LEnd)
                            self.BreakRight=min(self.BreakRight,pi.RStart)
                if e.Type==1:
                    if Candidate.AveBreak:
                        SumLeft+=e.Data.Begin
                        SumRight+=e.Data.End
                        BN+=1
                    else:
                        self.BreakLeft=max(self.BreakLeft,e.Data.Begin)
                        self.BreakRight=min(self.BreakRight,e.Data.End)
                if Candidate.AveBreak:
                    self.BreakLeft=SumLeft/BN
                    self.BreakRight=SumRight/BN

            if Candidate.CombineMode==1 or Candidate.CombineMode==2:
                self.BeginRange=(self.Begin,self.Begin)
                self.EndRange=(self.End,self.End)
                for e in self.Evidences:
                    if e.Type==0:
                        for pi in e.Data:
                            self.BeginRange=(min(self.BeginRange[0],pi.LStart),max(self.BeginRange[1],pi.LEnd))
                            self.EndRange=(min(self.EndRange[0],pi.RStart),max(self.EndRange[1],pi.REnd))
                    if e.Type==1:
                        self.BeginRange=(min(self.BeginRange[0],e.Data.Begin),max(self.BeginRange[1],e.Data.Begin))
                        self.EndRange=(min(self.EndRange[0],e.Data.End),max(self.EndRange[1],e.Data.End))
        if self.AveRange:
            SBegin=0
            SEnd=0
            for e in self.Evidences:
                SBegin+=e.Begin
                SEnd+=e.End
            self.Begin=int(SBegin/len(self.Evidences))
            self.End=int(SEnd/len(self.Evidences))

    def deductSVType(self):
        if len(self.Evidences)!=0:
            self.SVType=3#self.Evidences[0].SupportedSVType

    def addEvidence(self,e):
        if len(Evidences)==0:
            Evidences.append(e)
            self.deductSVType()
            self.calculateSpread()
    
    def mergeCandidate(self,other):#return: 0: not combine, 1: combine, -1: not overlap
        Overlap=calcOverlap(self.Begin,self.End,other.Begin,other.End)
        if Overlap==-1:
            return -1
        if self.SVType!=other.SVType:
            return 0
        if calcOverlap(self.BreakLeft,self.BreakRight,other.BreakLeft,other.BreakRight)<=0:#can't be zero, othersize will generate 0 length cnv
            return 0
        #MinLength=min(self.End-self.Begin,other.End-other.Begin)
        ToCombine=False
        if Candidate.CombineMode==0 or Candidate.CombineMode==3:
            Length=max(self.End-self.Begin,other.End-other.Begin)
            Overlap=Overlap/Length if Length!=0 else 0
            if Overlap>=Candidate.CombinePercentage:
                ToCombine=True
        else:
            if Candidate.CombineMode==1:
                CombineRange=Candidate.CombineRange
            else:
                if self.CombineRange==None:
                    self.CombineRange=min(int((self.End-self.Begin)*Candidate.RangePercentage),Candidate.MaxRange)
                    self.CombineRange=max(Candidate.MinRange,self.CombineRange)
                CombineRange=self.CombineRange
            if self.BeginRange[0]>=other.BeginRange[1]-CombineRange and self.BeginRange[1]<=other.BeginRange[0]+CombineRange and self.EndRange[0]>=other.EndRange[1]-CombineRange and self.EndRange[1]<=other.EndRange[0]+CombineRange:
                ToCombine=True
        if ToCombine:
            self.Evidences+=other.Evidences
            self.calculateSpread()
            return 1
        return 0
    
    def unifyEvidences(self):
        for e in self.Evidences:
            e.Begin=self.Begin
            e.End=self.End
            if Candidate.CombineMode==2:
                e.Begin=self.BreakLeft
                e.End=self.BreakRight
            if e.Type==1:
                e.Data.Begin=self.Begin
                e.Data.End=self.End
                e.Data.WBegin=int(self.Begin/g.RDWindowSize)
                e.Data.WEnd=int(self.End/g.RDWindowSize)
                if e.Data.WEnd<=e.Data.WBegin:
                    e.Data.WEnd=e.Data.WBegin+1

class Evidence:
    CombinePercentage=0.9
    def __init__(self,Type=0,Data=None,Begin=0xffffffff,End=0):
        self.Type=Type#Type: 0: DR cluster, Data is array of PairInfo;1: RD Interval, Data is a RDInterval
        self.Data=Data
        self.Begin=Begin
        self.End=End
        self.analyzeData()
    
    def setData(self, Type, Data):
        self.Type=Type
        self.Data=Data
        self.analyzeData()
    
    def calculateSpread(self):
        if self.Type==0:
            for pi in self.Data:
                self.Begin, self.End=(min(self.Begin,pi.Start),max(self.End,pi.End))
        elif self.Type==1:
            self.Begin=self.Data.Begin
            self.End=self.Data.End
    
    def analyzeData(self):
        if self.Data==None:
            return
        self.calculateSpread()
        self.SupportedSVType=None
        if self.Type==0:
            if self.Data!=None and len(self.Data)!=0:
                ClassOfD0=self.Data[0].classify()
                if ClassOfD0==3:
                    self.SupportedSVType=1
                elif ClassOfD0==4:
                    self.SupportedSVType=0
                self.Sample=self.Data[0].SampleIndex
        elif self.Type==1:
            self.SupportedSVType=self.Data.SupportedSVType
            self.Sample=self.Data.Sample

    #do not use this
    def combineEvidence(self, other):#DO NOT USE THIS! return: 0: not combine, 1: combine, -1: not overlap
        if self.Spread[0]<other.Spread[1] and self.Spread[0]>other.Spread[0]:
            Overlap=min(other.Spread[1]-self.Spread[0],self.Spread[1]-self.Spread[0])
        elif self.Spread[1]>other.Spread[0] and self.Spread[1]<other.Spread[0]:
            Overlap=self.Spread[1]-other.Spread[0]
        elif other.Spread[0]>self.Spread[0] and other.Spread[0]<self.Spread[1]:#self covers other
            Overlap=other.Spread[1]-other.Spread[0]
        else:
            Overlap=0
            return -1
        #MinLength=min(self.Spread[1]-self.Spread[0],other.Spread[1]-other.Spread[0])
        if sefl.Type==1:
            return 0
        MaxLength=max(self.Spread[1]-self.Spread[0],other.Spread[1]-other.Spread[0])
        Overlap=Overlap/MaxLength if MaxLength!=0 else 0
        if Overlap>=Evidence.CombinePercentage:
            self.Data+=other.Data
            self.Spread=(min(self.Spread[0],other.Spread[0]),max(self.Spread[1],other.Spread[1]))
            return 1
        return 0

def combineCandidateSets(Candidates1, Candidates2):
    NewCandidateSet=[]
    Candidates1+=Candidates2
    Candidates1.sort(key=lambda c:c.Begin)
    ni=-1
    for i in range(len(Candidates1)):
        if Candidates1[i]==None:
            continue
        NewCandidateSet.append(Candidates1[i])
        ni+=1
        Candidates1[i]=None
        for j in range(i+1,len(Candidates1)):
            if Candidates1[j]==None:
                continue
            mr=NewCandidateSet[ni].mergeCandidate(Candidates1[j])
            if mr==-1:
                break
            elif mr==1:
                Candidates1[j]=None
    return NewCandidateSet

