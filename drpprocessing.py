import globals as g
from utils import *

#def indexDPRs(TheContig):
#    TheContig.SortedDRPsBySampleI={}
#    for i in range(len(TheContig.RDWindows)):
#        TheContig.SortedDRPsBySampleI[i]=[]
#    for d in TheContig.DRPs:
#        TheContig.SortedDRPsBySampleI[g.SampleNameIndexes[d.SampleName]].append(d)
#    for i in range(len(TheContig.RDWindows)):
#        TheContig.SortedDRPsBySampleI[i].sort(key=lambda d: d.Start)

def processEvidenceWithDRPs(e, DRPs, MaxEnds):
    Start=findIndex(MaxEnds[e.Sample],e.Begin)
    End=findStartIndex(DRPs[e.Sample],e.End)
    for d in DRPs[e.Sample][Start:End]:
        if (d.SupportedVariantType ==0 and e.SupportedSVType> 0) or (d.SupportedVariantType==1 and e.SupportedSVType==0):
            continue 
        if inclusion((e.Begin,e.End),(d.Start,d.End))>2:#d include e or identical
            e.SupportedDRPs.append(d)

def processEvidencesWithDRPs(Es, TheContig):
    warn("processing DRPs...")
    MaxEnds=[]
    for Ds in TheContig.DRPs:
        MaxEnd=0
        MaxEnds.append([])
        for d in Ds:
            if d.End>MaxEnd:
                MaxEnd=d.End
            MaxEnds[-1].append(MaxEnd)
        MaxEnds[-1].sort()
    for s in Es:
        for e in s:
            processEvidenceWithDRPs(e,TheContig.DRPs,MaxEnds)
    return Es

def findStartIndex(DRPs,Value,s=0,e=-1):#return index with value, not necessarily the first
    if e==-1:
        e=len(DRPs)
    if e<=s+1:
        return s
    m=int((e+s)/2)
    if DRPs[m].Start==Value:
        return m
    elif DRPs[m].Start<Value:
        return findStartIndex(DRPs,Value,m+1,e)
    else:
        return findStartIndex(DRPs,Value,s,m)

def findIndex(Values,Value,s=0,e=-1):#return index with value, not necessarily the first
    if e==-1:
        e=len(Values)
    if e<=s+1:
        return s
    m=int((e+s)/2)
    if Values[m]==Value:
        return m
    elif Values[m]<Value:
        return findIndex(Values,Value,m+1,e)
    else:
        return findIndex(Values,Value,s,m)