from globals import *
from utils import *

class RDInterval:
    def __init__(self,WBegin,WEnd,ARD):
        self.WBegin=WBegin
        self.WEnd=WEnd
        self.AverageRD=ARD
        self.Begin=WBegin*RDWindowSize
        self.End=WEnd*RDWindowSize
        self.SupportedSVType=None#0:del, 1:insertion, 2:dup
        if self.AverageRD<0.1:
            self.CN=0
        elif self.AverageRD<0.75:
            self.CN=1
        else:
            self.CN=int(self.AverageRD+0.5)*2-1
        if self.CN==0 or self.CN==1:
            self.SupportedSVType=0
        elif self.CN>2:
            self.SupportedSVType=2

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


def makeRDIntervals(MixedRDRs):
    EvidenceIntervals=[]
    for i in range(len(MixedRDRs)):#for each sample
        CurrentRunStart=0
        CurrentRunRatio=1
        for j in range(1,len(MixedRDRs[i])):
            if sigDiff(MixedRDRs[i],j,CurrentRunRatio):
                if CurrentRunRatio!=1:
                    LastInterval=RDInterval(CurrentRunStart,j-1,CurrentRunRatio)
                    if LastInterval.CN!=2:
                        EvidenceIntervals.append(LastInterval)
                CurrentRunRatio=(MixedRDRs[i][j]+MixedRDRs[i][j-1])/2
                CurrentRunStart=j-1
            elif CurrentRunRatio!=1:
                CurrentRunRatio=(CurrentRunRatio*(j-CurrentRunStart)+MixedRDRs[i][j])/(j-CurrentRunStart+1)
        if CurrentRunRatio!=1:
            LastInterval=RDInterval(CurrentRunStart,len(MixedRDRs[i]),CurrentRunRatio)
            if LastInterval.CN!=2:
                EvidenceIntervals.append(LastInterval)
    return EvidenceIntervals

def makeRDICandidates(RDIs):
    Candidates=[]
    for i in RDIs:
        e=Evidence()
        e.setData(1,i)
        Candidates.append(Candidate([e]))
    Candidates=combineCandidateSets(Candidates,[])
    return Candidates

def analyzeRD(RDWindows,WindowsN,OccurredWindowsN):
    RDWindowAverages=[0]*WindowsN
    RDWindowSums=[0]*WindowsN
    SampleN=len(RDWindows)
    SampleSums=[0]*SampleN
    SampleAverages=[0]*SampleN
    for i in range(WindowsN):
        for j in range(SampleN):
            RDWindowSums[i]+=RDWindows[j][i]
            SampleSums[j]+=RDWindows[j][i]
        RDWindowAverages[i]=RDWindowSums[i]/SampleN
    for j in range(SampleN):
        SampleAverages[j]=SampleSums[j]/OccurredWindowsN
    MixedRDRs=[[0]*WindowsN]*SampleN#Mixed Read depth rate
    for i in range(SampleN):
        for j in range(WindowsN):
            SR=(RDWindows[i][j]/SampleAverages[i]) if SampleAverages[i]!=0 else 1
            WR=(RDWindows[i][j]/RDWindowAverages[j]) if RDWindowAverages[j]!=0 else 1
            MixedRDRs[i][j]=(SR-WR)/SampleN+WR
    RDICandidates=makeRDICandidates(makeRDIntervals(MixedRDRs))
    return RDICandidates