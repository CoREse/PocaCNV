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

import os
import psutil
process = psutil.Process(os.getpid())

WriteRDData=False
Contigs={"chr22"}#if not vacant, contain only those contigs

print(gettime()+"Joint calling started...", file=sys.stderr)
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
        print(gettime()+"Sample %s read. Memory usage:%.6sgb."%(SampleNames[-1],process.memory_info().vms/1024/1024/1024),file=sys.stderr)
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
        print(gettime()+"Sample %s read. Memory usage:%.6sgb."%(SampleNames[-1],process.memory_info().vms/1024/1024/1024),file=sys.stderr)
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

#EAFile=open("data/exclusion_regions.txt","r")
#readExcludedAreas(EAFile,ContigNameOccurred)
#EAFile.close()
'''
if WriteRDData:
    writeRDData(mygenome,ReferenceFile,SampleNames)
    exit(0)
'''

print("Memory usage: %.6sgb"%(process.memory_info().vms/1024/1024/1024),file=sys.stderr)
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
    print(gettime()+"%s analyzed. Memory usage: %.6sgb"%(c.Name, process.memory_info().vms/1024/1024/1024),file=sys.stderr)
print("Number of RDI candidates:%d"%(len(RDICandidates)),file=sys.stderr)
print("Memory usage: %.6sgb"%(process.memory_info().vms/1024/1024/1024),file=sys.stderr)

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
print("Number of candidates:%d"%(CCount),file=sys.stderr)
reportVCFHeader(sys.stdout)
for i in range(len(mygenome)):
    SVs=[]
    for C in Candidates[i]:
        SV=callSV(ReferenceFile,C,mygenome[i])
        if SV!="":
            SV.Chrom=mygenome[i].Name
            SVs.append(SV)
    SVs.sort(key=lambda s:s.BreakLeft)
    reportVCF(SVs,ReferenceFile,sys.stdout)
ReferenceFile.close()