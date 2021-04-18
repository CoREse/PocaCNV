import globals
g=globals
from utils import *
import math
from scipy.stats import poisson
from consts import *
import sys
from variants import CNV
from readdepth import cn2likely
import torch

def conditionP(O,Ei):
    return (Ei**O)*(math.e**(-Ei))/math.factorial(int(O))

import CGetRDScores
def getRDScores(Candidates,TheContig,ThreadN=1):
    print(gettime()+"Getting Scores...",file=sys.stderr)
    Scores=[]
    if g.EDataName!=None:
        EDF=open("%s/%s%s_%s_%s_EData.txt"%(g.EDataPath, g.EDataName,len(TheContig.SampleNames),g.RDWindowSize,TheContig.Name),"w")
    else:
        EDF=None
    #model=torch.load("ScoringTrain/Model032/Model032.pickle")
    if g.StatLocal:
        SampleReadCount=TheContig.SampleReadCount
    else:
        SampleReadCount=g.SampleReadCount
    Scores=CGetRDScores.CGetRDScores(Candidates,TheContig,ThreadN,SampleReadCount)
    print(gettime()+"Filting NN...",file=sys.stderr)
    model=torch.load("data/ScoringTrainModelData")
    for i in range(len(Candidates)):
        #Scores.append(getScore(Candidates[i],TheContig))
        passNN(Candidates[i],TheContig,Scores[i],(EDF,i),model)
    if EDF!=None:
        EDF.close()
    return Scores

def printEData(SegFileNNumber, TheContig, SiblingCount, E,CScore):
    SegNum=SegFileNNumber[1]
    ContigLength=TheContig.NLength
    #SiblingCount
    Start=E.Begin
    End=E.End
    Mu=E.Data.mu
    MuS=E.Data.mus
    PassConfidence=E.PassConfidence
    CN=E.Data.CN
    Confidence=E.Confidence
    #CScore=float(sl[10])
    Length=End-Start
    StartPortion=Start/ContigLength
    EndPortion=End/ContigLength
    HasSibling=1 if SiblingCount>0 else 0
    MultiSibling=1 if SiblingCount>1 else 0
    SampleCount=len(TheContig.SampleNames)
    SiblingRatio=SiblingCount/SampleCount
    SDRPN=len(E.SupportedDRPs)
    ContigSampleReadCoverage=TheContig.ContigSampleReadCounts[E.Sample]/TheContig.NLength
    SDRPRatio=SDRPN/ContigSampleReadCoverage if ContigSampleReadCoverage!=0 else 0
    HasSDRP=1 if SDRPN>0 else 0
    HasMultiSDRP=1 if SDRPN>1 else 0
    if (SegFileNNumber[0]!=None):
        print("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"%(SegFileNNumber[1],TheContig.NLength,SiblingCount,SiblingRatio,g.SampleNames[E.Sample],E.Begin,E.End,E.Data.mu,E.Data.mus,E.PassConfidence,E.Data.CN,E.Confidence,CScore,HasSDRP,HasMultiSDRP,SDRPRatio),file=SegFileNNumber[0])
    #DataItem=[SiblingRatio,HasSibling,MultiSibling,Length,StartPortion,EndPortion,Mu,MuS,PassConfidence,CN,Confidence,CScore,CNPriors[int(CN)] if int(CN)<len(CNPriors) else 0]
    DataItem=[SiblingRatio,HasSibling,MultiSibling,Length,StartPortion,EndPortion,Mu,MuS,PassConfidence,CN,Confidence,CScore,CNPriors[int(CN)] if int(CN)<len(CNPriors) else 0,HasSDRP,HasMultiSDRP,SDRPRatio]
    #print(DataItem,file=sys.stderr)
    #return torch.Tensor([SegFileNNumber[1],TheContig.NLength,SiblingCount,E.Begin,E.End,E.Data.mu,E.Data.mus,E.PassConfidence,E.Data.CN,E.Confidence,CScore,ChrNo,CNPriors[int(CN)] if int(CN)<len(CNPriors) else 0])
    return torch.Tensor(DataItem)

def getRDScore(C, TheContig):
    Score=0
    P=1
    CN2L=0
    for e in C.Evidences:
        e.PassConfidence=0
        mu=0
        mus=0
        if e.Data.mu==None:
            mu,mus=e.Data.calcMuMus(TheContig)
        else:
            mu,mus=(e.Data.mu,e.Data.mus)
        e.PassConfidence=1-cn2likely(e.Data,TheContig)
        CN2L=max(CN2L,1-e.PassConfidence)
        #v=e.Data.AverageRD/2.0*TheContig.MRMedians[e.Sample]*mu-mu#mu=lambda0*length, let averagerd*lambda0 be lambda
        #v=e.Data.AverageRD/2.0*mu-mu#mu=lambda0*length, let averagerd*lambda0 be lambda
        #mus=e.Data.AverageRD/2.0*mu
        eCN=e.Data.CN
        MP=0
        MCN=0
        CNPN=len(CNPriors)-1
        Pmus=0
        for i in range(CNPN):
            Pmus+=CNPriors[i]*poisson.pmf(mus,int(mu*i/2))
        if Pmus==0:
            MCN=eCN
            MP=1
        else:
            for CN in range(max(0,eCN-1),eCN+2):
                if CN>CNPN:
                    CN=CNPN
                Pmuscn=poisson.pmf(mus,int(mu*CN/2))*CNPriors[CN]
                Pd=Pmuscn/Pmus if Pmus!=0 else 0
                if Pd>MP:
                    MP=Pd
                    MCN=CN
            #print(Pd,Pmuscn,Pmus, CN, poisson.pmf(mus,int(mu*CN/2)),file=sys.stderr)
        e.Data.CN=MCN
        e.Confidence=MP
        #print(e.Confidence)
        #if e.Confidence>g.SampleConfidenceThreshold:
        #    P*=1-MP
        '''
        qint=poisson.interval(0.99,mu)#(nlambda-k(nlambda)^0.5,nlambda+k(nlambda)^0.5)
        if qint[0]<v<qint[1]:
            Score+=1
        qint=poisson.interval(0.999,mu)
        if qint[0]<v<qint[1]:
            Score+=0.5
        qint=poisson.interval(0.9999,mu)
        if qint[0]<v<qint[1]:
            Score+=0.5
        qint=poisson.interval(0.99999,mu)
        if qint[0]<v<qint[1]:
            Score+=1
        '''
    return 1-CN2L

def passNN(C,TheContig,Score,SegmentFileNNumber,Model):
    for e in C.Evidences:
        Traits=printEData(SegmentFileNNumber,TheContig,len(C.Evidences),e,Score)
        YP=Model(Traits)
        #print(Traits,YP,file=sys.stderr)
        #exit(0)
        if YP[1]>YP[0]:
            e.NN=1
        else:
            e.NN=0

def getScore(C,TheContig):
    return getRDScore(C,TheContig)
    Score=0
    RDScore=0
    U=0
    N=0
    for e in C.Evidences:
        if e.Type==0:
            if e.Data[0].NMapped!=2:
                continue
            for r in e.Data:
                U+=r.InsertionSize
            N+=len(e.Data)
        if e.Type==1:
            RDScore+=5
    Score+=RDScore
    #read-depth break interval score
    MeanRD=0
    IntervalWN=int(C.BreakRight/g.RDWindowSize)+1-int(C.BreakLeft/g.RDWindowSize)
    SampleN=len(TheContig.MixedRDRs)
    for i in range(int(C.BreakLeft/g.RDWindowSize),int(C.BreakRight/g.RDWindowSize)):
        MeanRD+=TheContig.MixedRDRs[SampleN-1][i]
    MeanRD/=IntervalWN
    if MeanRD>100:#reduce calculation
        MeanRD=100
    ECN=[]
    for i in range(11):
        ECN.append(globals.ERD*i/2)
    Support=0
    Other=0
    if C.SVType==0:
        Support=conditionP(MeanRD,ECN[0])*CNPriors[0]+conditionP(MeanRD,ECN[0])*CNPriors[1]
        for i in range(2,11):
            Other+=conditionP(MeanRD,ECN[i])*CNPriors[i]
    elif C.SVType==2:
        for i in range(3,11):
            Support+=conditionP(MeanRD,ECN[i])*CNPriors[i]
        for i in range(1,3):
            Other+=conditionP(MeanRD,ECN[i])*CNPriors[i]
    if Other!=0:
        Support/=Other
    if Support>1:
        Score+=200

    if N==0:
        return Score
    U=((float(U)/N)-g.MedianInsertionSize)/g.ISSD*float(N)**0.5
    Score+=U
    return Score

def getBreak(C):
    return (C.BreakLeft,C.BreakRight)

def getSVType(C):
    if C.SVType==0:
        return "DEL"
    elif C.SVType==1:
        return "INS"
    elif C.SVType==2:
        return "DUP"
    return "NONE"

def prefilters(C):
    if C.BreakRight<C.BreakLeft:
        return False
    return True

def getInvolvedSamples(C):
    Ss=set()
    for e in C.Evidences:
        Ss.add(e.Sample)
    return Ss

def callSV(ReferenceFile,C,TheContig,Score=None):
    SV=""
    if prefilters(C)!=True:
        return SV
    SVType=getSVType(C)
    if Score==None:
        Score=getScore(C,TheContig)
    BKL,BKR=getBreak(C)
    if Score>g.Parameters.ScoreThreshold:
        Alleles=set()
        Samples=[]
        Occured=set()
        for E in C.Evidences:
            #if E.Confidence<=g.Parameters.SampleConfidenceThreshold:
            #    continue
            if E.NN!=1:
                continue
            Ploidy=E.Data.Ploidy
            if E.Data.CN==Ploidy:
                continue
            if E.Data.CN<Ploidy:
                Alleles.add(0)
            else:
                Alleles.add(E.Data.CN-(1 if Ploidy==2 else 0))
        Alleles=list(Alleles)
        if len(Alleles)==0:
            return SV
        Alleles.sort()
        for E in C.Evidences:
            if E.Sample in Occured:# or E.Confidence<=g.Parameters.SampleConfidenceThreshold:
                continue
            Occured.add(E.Sample)
            SA=(0,0)
            if E.Data.CN==0:
                SA=(1,1)
            elif E.Data.CN==1:
                SA=(1,0)
            else:
                for i in range(len(Alleles)):
                    if E.Data.CN-1==Alleles[i]:
                        SA=(i+1,0)
                        break
            Samples.append((E.Sample,SA,E.Confidence))
        Samples.sort(key=lambda s:s[0])
        for i in range(len(Alleles)):
            Alleles[i]="<CN%s>"%Alleles[i]
        SV=CNV(BKL+1,BKR+1,False,Alleles,Samples,Score)
    return SV

def callSV_old(ReferenceFile,C,TheContig):
    SV=""
    if prefilters(C)!=True:
        return SV
    SVType=getSVType(C)
    Score=getScore(C,TheContig)
    BKL,BKR=getBreak(C)
    InvolvedSamples=getInvolvedSamples(C)
    InvolvedSamples=list(InvolvedSamples)
    InvolvedSamples.sort()
    if Score>0.999:
        SV+=TheContig.Name+":"+str(1+C.Begin)+"-"+TheContig.Name+":"+str(1+C.End)
        SV+=", "+SVType
        SV+=", Breakpoint:[%s,%s]"%(TheContig.Name+":"+str(1+BKL),TheContig.Name+":"+str(1+BKR))
        SV+=", Sample(s):"
        for SI in InvolvedSamples:
            SV+=" "+g.SampleNames[SI]
        SV+=", Score:%s"%Score
    return SV
