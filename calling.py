from globals import *
from utils import *
def getScore(C):
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

def callSV(ReferenceFile,C):
    SV=""
    SVType=getSVType(C)
    Score=getScore(C)
    BKL,BKR=getBreak(C)
    if Score>15:
        SV+=ReferenceFile.references[getTidByCord(C.Begin)]+":"+str(1+C.Begin-RefStartPos[getTidByCord(C.Begin)])+"-"+ReferenceFile.references[getTidByCord(C.End)]+":"+str(1+C.End-RefStartPos[getTidByCord(C.End)])
        SV+=", "+SVType
        SV+=", Breakpoint:[%s,%s]"%(ReferenceFile.references[getTidByCord(BKL)]+":"+str(1+BKL-RefStartPos[getTidByCord(BKL)]),ReferenceFile.references[getTidByCord(BKR)]+":"+str(1+BKR-RefStartPos[getTidByCord(BKR)]))
    return SV
