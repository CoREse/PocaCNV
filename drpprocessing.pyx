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

def processSampleEvidencesWithDRPs(SEs, DRPs, MaxEnds,SampleName):
    warn("Processing DRPs for %s"%SampleName)
    for e in SEs:
        Start=findIndex(MaxEnds,e.Begin)
        End=findStartIndex(DRPs,e.End)
        e.SupportedDRPCount=0
        for i in range(len(DRPs[Start:End])):
            d=DRPs[i]
            if (d.SupportedVariantType ==0 and e.SupportedSVType> 0) or (d.SupportedVariantType==1 and e.SupportedSVType==0):
                continue 
            if inclusion((e.Begin,e.End),(d.Start,d.End))>2:#d include e or identical
                #e.SupportedDRPs.append(i)
                e.SupportedDRPCount+=1
    return SEs

def getSampleEvidencesDRPCounts(SEs, DRPs, MaxEnds,SampleName):
    warn("Processing DRPs for %s"%SampleName)
    DRPCounts=[0]*len(SEs)
    for k in range(len(SEs)):
        e=SEs[k]
        Start=findIndex(MaxEnds,e.Begin)
        End=findStartIndex(DRPs,e.End)
        for i in range(len(DRPs[Start:End])):
            d=DRPs[i]
            if (d.SupportedVariantType ==0 and e.SupportedSVType> 0) or (d.SupportedVariantType==1 and e.SupportedSVType==0):
                continue 
            if inclusion((e.Begin,e.End),(d.Start,d.End))>2:#d include e or identical
                #e.SupportedDRPs.append(i)
                DRPCounts[k]+=1
    return DRPCounts

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
    delPool()
    pool.terminate()
    #ctx=mp.get_context("fork")
    #pool=ctx.Pool(g.ThreadN)
    #addPool(pool)
    #args=[]
    #for i in range(SampleN):
    #    args.append((Es[i],TheContig.DRPs[i],MaxEnds[i],TheContig.SampleNames[i]))
    #addPool(pool)
    #for i in range(SampleN):
    #    NEs.append(processSampleEvidencesWithDRPs(Es[i],TheContig.DRPs[i],MaxEnds[i],TheContig.SampleNames[i]))
    for i in range(SampleN):
        Counts=getSampleEvidencesDRPCounts(Es[i],TheContig.DRPs[i],MaxEnds[i],TheContig.SampleNames[i])
        for j in range(len(Counts)):
            Es[i][j].SupportedDRPCount=Counts[j]
    #Counts=pool.starmap(getSampleEvidencesDRPCounts,args)
    #for i in range(len(Es)):
    #    for j in range(len(Es[i])):
    #        Es[i][j].SupportedDRPCount=Counts[i][j]
    #        for k in range(len(NEs[i][j].SupportedDRPs)):
    #            NEs[i][j].SupportedDRPs[k]=TheContig.DRPs[i][NEs[i][j].SupportedDRPs[k]]
    print(gettime()+"Evidences processed. "+getMemUsage(),file=sys.stderr)
    return Es

def findStartIndexRecur(DRPs,Value,s=0,e=-1):#return index with value, not necessarily the first
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

def findIndexRecur(Values,Value,s=0,e=-1):#return index with value, not necessarily the first
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

def findStartIndex(DRPs,Value,s=0,e=-1):#return index with value, not necessarily the first
    if e==-1:
        e=len(DRPs)
    while e>s+1:
        m=int((e+s)/2)
        if DRPs[m].Start==Value:
            return m
        elif DRPs[m].Start<Value:
            s=m+1
            continue
        else:
            e=m
            continue
    return s

def findIndex(Values,Value,s=0,e=-1):#return index with value, not necessarily the first
    if e==-1:
        e=len(Values)
    while e>s+1:
        m=int((e+s)/2)
        if Values[m]==Value:
            return m
        elif Values[m]<Value:
            s=m+1
            continue
        else:
            e=m
            continue
    return s