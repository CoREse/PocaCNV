import globals as g
from utils import *
import multiprocessing as mp

#def indexDPRs(TheContig):
#    TheContig.SortedDRPsBySampleI={}
#    for i in range(len(TheContig.RDWindows)):
#        TheContig.SortedDRPsBySampleI[i]=[]
#    for d in TheContig.DRPs:
#        TheContig.SortedDRPsBySampleI[g.SampleNameIndexes[d.SampleName]].append(d)
#    for i in range(len(TheContig.RDWindows)):
#        TheContig.SortedDRPsBySampleI[i].sort(key=lambda d: d.Start)

def processSampleEvidencesWithDRPs(SEs, DRPs, MaxEnds):
    for e in SEs:
        Start=findIndex(MaxEnds,e.Begin)
        End=findStartIndex(DRPs,e.End)
        for d in DRPs[Start:End]:
            if (d.SupportedVariantType ==0 and e.SupportedSVType> 0) or (d.SupportedVariantType==1 and e.SupportedSVType==0):
                continue 
            if inclusion((e.Begin,e.End),(d.Start,d.End))>2:#d include e or identical
                e.SupportedDRPs.append(d)
    return SEs

def getMaxEnds(DRPs):
    MaxEnd=0
    MaxEnds=[]
    for d in DRPs:
        if d.End>MaxEnd:
            MaxEnd=d.End
        MaxEnds.append(MaxEnd)
    MaxEnds.sort()
    return MaxEnds

def processEvidencesWithDRPs(Es, TheContig):
    warn("processing DRPs...")
    MaxEnds=[]
    SampleN=len(Es)
    ctx=mp.get_context("fork")
    pool=ctx.Pool(g.ThreadN)
    args=[]
    for i in range(SampleN):
        args.append((TheContig.DRPs[i]))
    addPool(pool)
    MaxEnds=pool.map(getMaxEnds,args)
    print(gettime()+"Max Ends made. "+getMemUsage(),file=sys.stderr)
    args=[]
    for i in range(SampleN):
        args.append((Es[i],TheContig.DRPs[i],MaxEnds[i]))
    addPool(pool)
    Es=pool.starmap(processSampleEvidencesWithDRPs,args)
    print(gettime()+"Evidences processed. "+getMemUsage(),file=sys.stderr)
    delPool()
    pool.terminate()
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