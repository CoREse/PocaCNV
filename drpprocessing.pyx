import globals as g
from utils import *
import multiprocessing as mp
cimport cython
from cython.parallel import parallel,prange
from libc.stdlib cimport malloc,free,qsort

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
        for i in range(Start,End):
            d=DRPs[i]
            if (d.SupportedVariantType ==0 and e.SupportedSVType> 0) or (d.SupportedVariantType==1 and e.SupportedSVType==0):
                continue 
            if inclusion((e.Begin,e.End),(d.Start,d.End))>2:#d include e or identical
                #e.SupportedDRPs.append(i)
                DRPCounts[k]+=1
    return DRPCounts

cdef int CfindIndex(int *Values,int Value,int s,int e) nogil:#return index with value, not necessarily the first
    cdef int m
    while e>s+1:
        m=<int>((e+s)/2)
        if Values[m]==Value:
            return m
        elif Values[m]<Value:
            s=m+1
            continue
        else:
            e=m
            continue
    return s

cdef int Cinclusion(int s1, int e1, int s2, int e2) nogil:#return 0: not overlapped, 1: overlapped, 2: In1 include In2, 3: In2 inlcude In1, 4: identical
    if s1==s2 and e1==e2:
        return 4
    elif s1>=s2 and e1<=e2:
        return 3
    elif s1<=s2 and e1>=e2:
        return 2
    elif s1>=s2 and s1<=e2 or e1>=s2 and e1<=e2:
        return 1
    return 0

cdef int * CgetSampleEvidencesDRPCounts(int * Counts, int NSE, int * SEStarts, int *SEEnds, int *SESVs, int SNDRP, int * SDRPStarts, int * SDRPEnds, int * SDRPSVs, int * SMaxEnds) nogil:
    cdef int i,j
    cdef int Start, End
    for i from 0<=i<NSE:
        Counts[i]=0
        Start=CfindIndex(SMaxEnds,SEStarts[i],0,SNDRP)
        End=CfindIndex(SDRPStarts,SEEnds[i],0,SNDRP)
        for j from Start<=j<End:
            if (SDRPSVs[j] ==0 and SESVs[i]> 0) or (SDRPSVs[j]==1 and SESVs[i]==0):
                continue 
            if Cinclusion(SEStarts[i],SEEnds[i],SDRPStarts[j],SDRPEnds[j])>2:#d include e or identical
                #e.SupportedDRPs.append(i)
                Counts[i]+=1
    return Counts

def getMaxEnds(DRPs):
    MaxEnd=0
    MaxEnds=[]
    for d in DRPs:
        if d.End>MaxEnd:
            MaxEnd=d.End
        MaxEnds.append(MaxEnd)
    #MaxEnds.sort()
    return MaxEnds

def processEvidencesWithDRPs(Es, TheContig):
    warn("processing DRPs...")
    SampleN=len(Es)
    cdef int CSampleN=SampleN
    cdef int ** MaxEnds=<int **>malloc(SampleN*sizeof(int*))
    cdef int i,j
    cdef int NDRP,MaxEnd=0,MaxNE=0

    cdef int** DRPStarts=<int **>malloc(SampleN*sizeof(int*))
    cdef int** DRPEnds=<int **>malloc(SampleN*sizeof(int*))
    cdef int** DRPSVs=<int **>malloc(SampleN*sizeof(int*))
    cdef int** EStarts=<int **>malloc(SampleN*sizeof(int*))
    cdef int** EEnds=<int **>malloc(SampleN*sizeof(int*))
    cdef int** ESVs=<int **>malloc(SampleN*sizeof(int*))
    cdef int * NDRPs=<int *>malloc(SampleN*sizeof(int))
    cdef int * NEs=<int *>malloc(SampleN*sizeof(int))
    for i in range(SampleN):
        NDRPs[i]=len(TheContig.DRPs[i])
        MaxEnds[i]=<int*>malloc(NDRPs[i]*sizeof(int))
        DRPStarts[i]=<int*>malloc(NDRPs[i]*sizeof(int))
        DRPEnds[i]=<int*>malloc(NDRPs[i]*sizeof(int))
        DRPSVs[i]=<int*>malloc(NDRPs[i]*sizeof(int))
        for j in range(NDRPs[i]):
            DRPStarts[i][j]=TheContig.DRPs[i][j].Start
            DRPEnds[i][j]=TheContig.DRPs[i][j].End
            DRPSVs[i][j]=TheContig.DRPs[i][j].SupportedVariantType

        NEs[i]=len(Es[i])
        if MaxNE<NEs[i]:
            MaxNE=NEs[i]
        EStarts[i]=<int*>malloc(NEs[i]*sizeof(int))
        EEnds[i]=<int*>malloc(NEs[i]*sizeof(int))
        ESVs[i]=<int*>malloc(NEs[i]*sizeof(int))
        for j in range(NEs[i]):
            EStarts[i][j]=Es[i][j].Begin
            EEnds[i][j]=Es[i][j].End
            ESVs[i][j]=Es[i][j].SupportedSVType
    for i in range(SampleN):
        MaxEnd=0
        for j in range(NDRPs[i]):
            if MaxEnd<DRPEnds[i][j]:
                MaxEnd=TheContig.DRPs[i][j].End
            MaxEnds[i][j]=MaxEnd
    cdef int * Counts=<int*>malloc(sizeof(int)*MaxNE)
    for i in range(SampleN):
        CgetSampleEvidencesDRPCounts(Counts,NEs[i],EStarts[i],EEnds[i],ESVs[i],NDRPs[i],DRPStarts[i],DRPEnds[i],DRPSVs[i],MaxEnds[i])
        for j in range(NEs[i]):
            Es[i][j].SupportedDRPCount=Counts[j]
            #print(i,j,Es[i][j].SupportedDRPCount)
    free(Counts)
    for i in range(SampleN):
        free(MaxEnds[i])
        free(DRPStarts[i])
        free(DRPEnds[i])
        free(DRPSVs[i])
        free(EStarts[i])
        free(EEnds[i])
        free(ESVs[i])
    free(NDRPs)
    free(NEs)
    free(MaxEnds)
    free(DRPStarts)
    free(DRPEnds)
    free(DRPSVs)
    free(EStarts)
    free(EEnds)
    free(ESVs)
    
    print(gettime()+"Evidences processed. "+getMemUsage(),file=sys.stderr)
    return Es

def processEvidencesWithDRPsOld(Es, TheContig):
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
    pool.close()
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