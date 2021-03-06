import pysam
import sys
import time
from globals import *
from calling import *
from readpair import *
from readdepth import *
from filters import *

print("Joint calling started...", file=sys.stderr)
ReferenceFile=pysam.FastaFile(sys.argv[1])
PosCount=0
for tid in range(ReferenceFile.nreferences):
    RefInd[ReferenceFile.references[tid]]=tid
    RefStartPos.append(PosCount)
    PosCount+=ReferenceFile.lengths[tid]
RefLength=PosCount
WindowsN=int(RefLength/RDWindowSize)+1
RDWindows=[]

print("Reading samples...",file=sys.stderr)
SampleIndex=0
OccurredWindowsN=0
OccurredContigs={}
ContigNameOccurred={}

PairCount=0
ReadCount=0
LCount=0
RCount=0
UnmappedCount=0
for i in range(2,len(sys.argv)):
    SamFile=pysam.AlignmentFile(sys.argv[i],"rb",reference_filename=sys.argv[1])
    RDWindows.append([0]*WindowsN)
    for read in SamFile:
        ReadCount+=1
        read: pysam.AlignedSegment
        if read.is_unmapped:
            UnmappedCount+=1
        if read.is_paired:
            if read.query_name in UnpairedReads:
                OccurredContigs[read.reference_id]=True
                Pair=PairInfo(UnpairedReads.pop(read.query_name),read,i-2)
                PairCount+=1
                if Pair.isDiscordant():
                    DiscordantReads.append(Pair)
                #RDWindows[SampleIndex][int((Pair.Start+Pair.End)/2/100)]+=1
                if Pair.hasFlags(PairInfo.LMapped):
                    RDWindows[SampleIndex][int((Pair.LStart+Pair.LEnd)/2/RDWindowSize)]+=1
                    LCount+=1
                if Pair.hasFlags(PairInfo.RMapped):
                    RDWindows[SampleIndex][int((Pair.RStart+Pair.REnd)/2/RDWindowSize)]+=1
                    RCount+=1
            else:
                UnpairedReads[read.query_name]=read
    SamFile.close()
    SampleIndex+=1
#print(ReadCount,PairCount,LCount,RCount,UnmappedCount,file=sys.stderr)
#exit(0)

if len(RDWindows)>1:
    RDWindows.append([0]*WindowsN)#last sample is the s/"sum sample"/"average sample", because sum sample will significantly influence the WR variable
    SumI=len(RDWindows)-1
    for i in range(WindowsN):
        for j in range(SumI):
            RDWindows[SumI][i]+=RDWindows[j][i]
        RDWindows[SumI][i]/=SumI

for c in OccurredContigs:
    OccurredWindowsN+=ReferenceFile.lengths[c]/RDWindowSize+1
    ContigNameOccurred[ReferenceFile.references[c]]=True

EAFile=open("data/exclusion_regions.txt","r")
readExcludedAreas(EAFile,ContigNameOccurred)
EAFile.close()

print("Samples read, calculating RD data...", file=sys.stderr)
RDICandidates=analyzeRD(RDWindows,WindowsN,OccurredWindowsN)
RDICandidates=filtExcludedAreas(RDICandidates)
print("Number of RDI candidates:%d"%(len(RDICandidates)),file=sys.stderr)

print("RD data handled, calculating DRP data...", file=sys.stderr)
DiscordantReads=filtExcludedAreas(DiscordantReads)
DiscordantReads.sort(key=lambda r:r.Start)
DRClusters=cluster(DiscordantReads)
print("Number of dr clusters:%d"%(len(DRClusters)),file=sys.stderr)

DRCandidates=makeDRCCandidates(DRClusters)
print("Number of dr candidates:%d"%(len(DRCandidates)),file=sys.stderr)

Candidates=combineCandidateSets(DRCandidates,RDICandidates)
Candidates=filtExcludedAreas(Candidates)

#analyzeRD
#analyzeDR
#combineEvidence
print("Number of candidates:%d"%(len(Candidates)),file=sys.stderr)
for C in Candidates:
    SV=callSV(ReferenceFile,C)
    if SV!="":
        print(SV)

ReferenceFile.close()