from globals import *
import globals
from utils import *
import math
from scipy.stats import poisson

def conditionP(O,Ei):
    return (Ei**O)*(math.e**(-Ei))/math.factorial(int(O))

def getRDScore(C, TheContig):
    Score=0
    for e in C.Evidences:
        mu=0
        length=int(e.End/globals.RDWindowSize)-int(e.Begin/globals.RDWindowSize)
        for i in range(int(e.Begin/globals.RDWindowSize),int(e.End/globals.RDWindowSize)):
            mu+=TheContig.RDWindowStandards[i]
        mu=int(mu)
        v=e.Data.AverageRD/2.0*TheContig.MRMedians[e.Sample]*mu-mu#mu=lambda0*length, let averagerd*lambda0 be lambda
        qint=poisson.interval(0.9,mu)#(nlambda-k(nlambda)^0.5,nlambda+k(nlambda)^0.5)
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
    return Score

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
    IntervalWN=int(C.BreakRight/RDWindowSize)+1-int(C.BreakLeft/RDWindowSize)
    SampleN=len(TheContig.MixedRDRs)
    for i in range(int(C.BreakLeft/RDWindowSize),int(C.BreakRight/RDWindowSize)):
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
    U=((float(U)/N)-MedianInsertionSize)/ISSD*float(N)**0.5
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

def callSV(ReferenceFile,C,TheContig):
    SV=""
    if prefilters(C)!=True:
        return SV
    SVType=getSVType(C)
    Score=getScore(C,TheContig)
    BKL,BKR=getBreak(C)
    InvolvedSamples=getInvolvedSamples(C)
    InvolvedSamples=list(InvolvedSamples)
    InvolvedSamples.sort()
    if Score>3:
        SV+=TheContig.Name+":"+str(1+C.Begin)+"-"+TheContig.Name+":"+str(1+C.End)
        SV+=", "+SVType
        SV+=", Breakpoint:[%s,%s]"%(TheContig.Name+":"+str(1+BKL),TheContig.Name+":"+str(1+BKR))
        SV+=", Sample(s):"
        for SI in InvolvedSamples:
            SV+=" "+g.SampleNames[SI]
    return SV
