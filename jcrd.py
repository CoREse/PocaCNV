import pysam
import sys
import time
import globals
from calling import *
from readpair import *
from readdepth import *
from filters import *
from contig import *
from array import array
from report import *
from joint import uniformlyCombine

import os
import psutil
from multiprocessing import Manager,Pool
process = psutil.Process(os.getpid())

g=globals
g.Processes.append(process)
def getPid(i):
    return os.getpid()
if ThreadN!=1:
    globals.Pool=Pool(globals.ThreadN)
    globals.Manager=Manager()
    g.Processes.append(psutil.Process(g.Manager._process.ident))
    Pids=g.Pool.map(getPid,range(g.ThreadN))
    for pid in Pids:
        g.Processes.append(psutil.Process(pid))

def getMemUsage():
    vms=0
    rss=0
    for p in g.Processes:
        rss+=p.memory_info().rss
        vms+=p.memory_info().rss
    return "Memroy usage:%.5sgb(rss),%.5sgb(vms)."%(rss/1024/1024/1024,vms/1024/1024/1024)

WriteRDData=False
Contigs={"chr22"}#if not vacant, contain only those contigs
ExcludeRegionsFileName=None
#ExcludeRegionsFileName="data/mdust-v28-p1.bed"
#ExcludeRegionsFileName="data/exclusion_regions.txt"

print(gettime()+"Joint calling started...", file=sys.stderr)
print(gettime()+"Reading reference...",file=sys.stderr)
ReferenceFile=pysam.FastaFile(sys.argv[1])
PosCount=0
mygenome=Genome(ReferenceFile.filename)
for tid in range(ReferenceFile.nreferences):
    RefInd[ReferenceFile.references[tid]]=tid
    RefStartPos.append(PosCount)
    PosCount+=ReferenceFile.lengths[tid]
    if len(Contigs)!=0:
        if ReferenceFile.references[tid] not in Contigs:
            continue
    c=Contig(ReferenceFile.references[tid],int(ReferenceFile.lengths[tid]/globals.RDWindowSize)+(1 if ReferenceFile.lengths[tid]%globals.RDWindowSize!=0 else 0))
    mygenome.RefID.append(tid)
    mygenome.append(c)
RefLength=PosCount
print(gettime()+"Reference %s read. Length:%s, Contigs:%s."%(sys.argv[1],PosCount,len(mygenome)),file=sys.stderr)
print(gettime()+"Reading samples...",file=sys.stderr)
SampleIndex=0

PairCount=0
ReadCount=0
LCount=0
RCount=0
UnmappedCount=0
SampleNames=[]

for i in range(2,len(sys.argv)):
    if sys.argv[i].split(".")[-1]=="rdf":
        readRDData(mygenome,SampleNames,sys.argv[i])
        #OccurredWindowsN=len(RDWindows[0])
        print(gettime()+"Sample %s read. %s"%(SampleNames[-1],getMemUsage()),file=sys.stderr)
    else:
        SamFile=pysam.AlignmentFile(sys.argv[i],"rb",reference_filename=sys.argv[1])
        SampleNames.append(sys.argv[i].split("/")[-1].split("\\")[-1])
        mygenome.addSample(SampleNames[-1])
        ReadCount=0
        for read in SamFile:
            ReadCount+=1
            read: pysam.AlignedSegment
            if read.is_unmapped:
                UnmappedCount+=1
            else:
                if len(Contigs)!=0:
                    CID=mygenome.getContigID(read.tid)
                    if CID==-1:
                        continue
                    mygenome[CID].RDWindows[SampleIndex][int((int((read.reference_start+read.reference_end)/2))/globals.RDWindowSize)]+=1
                    ReadCount+=1
                else:
                    mygenome[read.tid].RDWindows[SampleIndex][int((int((read.reference_start+read.reference_end)/2))/globals.RDWindowSize)]+=1
                    ReadCount+=1
        SamFile.close()
        g.SampleReadCount.append(ReadCount)
        print(gettime()+"Sample %s read. %s"%(SampleNames[-1],getMemUsage()),file=sys.stderr)
    if WriteRDData:
        print(gettime()+"Storing rd data for %s..."%(SampleNames[-1]),file=sys.stderr)
        writeSampleRDData(mygenome,SampleNames[-1],SampleIndex)
    SampleIndex+=1
#print(ReadCount,PairCount,LCount,RCount,UnmappedCount,file=sys.stderr)
#exit(0)
globals.SampleNames=SampleNames
calclulateSequenceDepthRatio()
'''
if len(RDWindows)>1:
    RDWindows.append([0]*WindowsN)#last sample is the s/"sum sample"/"average sample", because sum sample will significantly influence the WR variable
    SumI=len(RDWindows)-1
    for i in range(WindowsN):
        for j in range(SumI):
            RDWindows[SumI][i]+=RDWindows[j][i]
        RDWindows[SumI][i]/=SumI
'''

if ExcludeRegionsFileName!=None:
    EAFile=open(ExcludeRegionsFileName,"r")
    ExcludedAreasByContig=readExcludedAreas(EAFile,ReferenceFile)
    EAFile.close()
'''
if WriteRDData:
    writeRDData(mygenome,ReferenceFile,SampleNames)
    exit(0)
'''

print("%s"%(getMemUsage()),file=sys.stderr)
print(gettime()+"Samples read, calculating RD data...", file=sys.stderr)

#for c in mygenome:
#    analyzeRD(c.RDWindows,c.Length,c, True)
#    print(gettime()+"%s analyzed. Memory usage: %.6sgb"%(c.Name, process.memory_info().vms/1024/1024/1024),file=sys.stderr)

#print("Memory usage: %.6sgb"%(process.memory_info().vms/1024/1024/1024),file=sys.stderr)
#writeMixedRDData(mygenome,ReferenceFile,SampleNames)
#exit(0)

RDICandidates=[]
for c in mygenome:
    c.RDICandidates=analyzeRD(c.RDWindows,c.Length,c)
    #c.RDICandidates=filtExcludedAreas(c.RDICandidates)
    RDICandidates.append(c.RDICandidates)
    print(gettime()+"%s analyzed. %s"%(c.Name, getMemUsage()),file=sys.stderr)
print("Number of RDI candidates:%d"%(len(RDICandidates)),file=sys.stderr)
print("%s"%(getMemUsage()),file=sys.stderr)

Candidates=RDICandidates
"""
for c in mygenome:
    Candidates.append(combineCandidateSets(DRCandidates,RDICandidates))
    Candidates[-1]=filtExcludedAreas(Candidates[-1])
"""
#analyzeRD
#analyzeDR
#combineEvidence
CCount=0
for cs in Candidates:
    CCount+=len(cs)
print(gettime()+"Number of candidates:%d. Filtering candidates in LPR..."%CCount,file=sys.stderr)
if ExcludeRegionsFileName!=None:
    for i in range(len(Candidates)):
        Candidates[i]=filtExcludedAreas(Candidates[i],ExcludedAreasByContig,mygenome[i])
CCount=0
for cs in Candidates:
    CCount+=len(cs)
print(gettime()+"Number of filtered candidates:%d. Uniformly combining..."%(CCount),file=sys.stderr)
for i in range(len(Candidates)):
    Candidates[i]=uniformlyCombine(Candidates[i],mygenome[i])
CCount=0
for cs in Candidates:
    CCount+=len(cs)
print(gettime()+"Number of uniformly combined candidates:%d. CNV calling..."%(CCount),file=sys.stderr)
reportVCFHeader(sys.stdout)
for i in range(len(mygenome)):
    SVs=[]
    print(gettime()+"Calling CNV for %s"%mygenome[i].Name,file=sys.stderr)
    for C in Candidates[i]:
        SV=callSV(ReferenceFile,C,mygenome[i])
        if SV!="":
            SV.Chrom=mygenome[i].Name
            SVs.append(SV)
    SVs.sort(key=lambda s:s.BreakLeft)
    reportVCF(SVs,ReferenceFile.fetch(mygenome[i].Name),sys.stdout)
ReferenceFile.close()
print(gettime()+"All done.",file=sys.stderr)