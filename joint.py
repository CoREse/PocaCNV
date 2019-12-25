from readdepth import *
from utils import *
import statistics
import globals as g

class CandidateCluster:
    ClusterSideDistance=500
    ClusterConcDistance=300
    ContainedRatio=0.95
    def __init__(self,Contig):
        self.Candidates=[]
        self.Begin=0
        self.End=0
        self.SVType=None
        self.Contig=Contig
    def calc(self):
        if len(self.Candidates)==0:
            return
        self.Begin=self.Candidates[0].Begin
        self.End=self.Candidates[0].End
        self.SVType=self.Candidates[0].SVType
        for C in self.Candidates[1:]:
            self.Begin=min(self.Begin,C.Begin)
            self.End=max(self.End,C.End)
        self.Length=self.End-self.Begin
    def take(self,C):
        if self.Candidates==[]:
            self.Candidates.append(C)
            self.calc()
            return True
        if self.SVType!=C.SVType:
            return 0
        OL=calcOverlap(self.Begin,self.End,C.Begin,C.End)
        if OL==-1:
            return -1
            if self.Begin>=C.End:
                if self.Begin-C.End>CandidateCluster.ClusterConcDistance:
                    return -1
            else:
                #self.End<=C.Begin:
                if C.Begin-self.End>CandidateCluster.ClusterConcDistance:
                    return -2
            self.Candidates.append(C)
            self.calc()
            return True
        else:
            SContained=False
            CContained=False
            if abs(self.Length-OL)<=CandidateCluster.ClusterSideDistance or 0.95<OL/self.Length<1/0.95:
                SContained=True
            elif abs(C.End-C.Begin-OL)<=CandidateCluster.ClusterSideDistance or 0.95<OL/C.End-C.Begin<1/0.95:
                CContained=True
            if SContained and CContained:
                self.Candidates.append(C)
                self.calc()
            return True
    
    def makeNew(self):
        BeginEdge=[]
        EndEdge=[]
        '''
        for C in self.Candidates:
            if C.Begin-self.Begin<=CandidateCluster.ClusterSideDistance:
                BeginEdge.append(C.Begin)
            if self.End-C.End<=CandidateCluster.ClusterSideDistance:
                EndEdge.append(C.End)'''
        for C in self.Candidates:
            BeginEdge.append(C.Begin)
            EndEdge.append(C.End)
        WBegin=int(statistics.median(BeginEdge)/g.RDWindowSize)
        WEnd=int(statistics.median(EndEdge)/g.RDWindowSize)
        MixedRDRs=self.Contig.MixedRDRs
        Evidences=[]
        m=0
        while m< len(self.Candidates):
            C=self.Candidates[m]
            k=0
            while k< len(C.Evidences):
                E=C.Evidences[k]
                Sample=E.Sample
                Ave=0
                for i in range(WBegin,WEnd):
                    Ave+=MixedRDRs[Sample][i]
                Ave/=WEnd-WBegin
                Interval=RDInterval(Sample,WBegin,WEnd,Ave,self.Contig)
                if cn2filter(Interval,self.Contig):
                    e=Evidence()
                    e.setData(1,Interval)
                    Evidences.append(e)
                    C.Evidences.pop(k)
                else:
                    k+=1
            if len(C.Evidences)==0:
                self.Candidates.pop(m)
            else:
                m+=1
        if len(Evidences)==0:
            return (None,self.Candidates)
        else:
            NewCan=Candidate(Evidences)
            return(NewCan,self.Candidates)

def uniformlyCombine(Candidates, Contig):
    Candidates.sort(key=lambda c:c.Begin)
    Combined=[]
    while len(Candidates)>0:
        Cluster=CandidateCluster(Contig)
        i=0
        while i< len(Candidates):
            C=Candidates[i]
            Res=Cluster.take(C)
            if Res==True:
                Candidates.pop(i)
            elif Res==-2:
                break
            else:
                i+=1
        i=0
        while i<len(Candidates):#after second time there is very unlikely other candidates that can be taken
            C=Candidates[i]
            Res=Cluster.take(C)
            if Res==True:
                Candidates.pop(i)
            elif Res==-2:
                break
            else:
                i+=1
        NewCandidate,Others=Cluster.makeNew()
        Others.sort(key=lambda c:c.Begin)
        i=0
        j=0
        while i<len(Others):
            while j <len(Candidates):
                if Others[i].Begin<Candidates[j].Begin:
                    Candidates.insert(j,Others[i])
                    break
                j+=1
            i+=1
        if NewCandidate!=None:
            Combined.append(NewCandidate)
        else:
            Combined.append(Candidates[0])
            Candidates.pop(0)
    return Combined

            
    
