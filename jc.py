import pysam
import sys
import time
from globals import *
from calling import *
from readpair import *
from readdepth import *

print("Joint calling started...", file=sys.stderr)
ReferenceFile=pysam.FastaFile(sys.argv[1])
PosCount=0
for tid in range(ReferenceFile.nreferences):
    RefStartPos.append(PosCount)
    PosCount+=ReferenceFile.lengths[tid]
RefLength=PosCount
WindowsN=int(RefLength/RDWindowSize)+1
RDWindows=[]

print("Reading samples...",file=sys.stderr)
SampleIndex=0
OccurredWindowsN=0
OccurredContigs=[]
for i in range(2,len(sys.argv)):
    SamFile=pysam.AlignmentFile(sys.argv[i],"rb",reference_filename=sys.argv[1])
    RDWindows.append([0]*WindowsN)
    for read in SamFile:
        read: pysam.AlignedSegment
        if read.is_paired:
            if read.query_name in UnpairedReads:
                OccurredContigs.append(read.reference_id)
                Pair=PairInfo(UnpairedReads.pop(read.query_name),read)
                if Pair.isDiscordant():
                    DiscordantReads.append(Pair)
                RDWindows[SampleIndex][int((Pair.Start+Pair.End)/2/100)]+=1
            else:
                UnpairedReads[read.query_name]=read
    SamFile.close()
    SampleIndex+=1
for c in OccurredContigs:
    OccurredWindowsN+=ReferenceFile.lengths[c]/RDWindowSize+1

print("Samples read, calculating RD data...", file=sys.stderr)
RDICandidates=analyzeRD(RDWindows,WindowsN,OccurredWindowsN)
print("Number of RDI candidates:%d"%(len(RDICandidates)),file=sys.stderr)

print("RD data handled, calculating DRP data...", file=sys.stderr)
DiscordantReads.sort(key=lambda r:r.Start)
DRClusters=cluster(DiscordantReads)
print("Number of dr clusters:%d"%(len(DRClusters)),file=sys.stderr)

DRCandidates=makeDRCCandidates(DRClusters)
print("Number of dr candidates:%d"%(len(DRCandidates)),file=sys.stderr)

Candidates=combineCandidateSets(DRCandidates,RDICandidates)

#analyzeRD
#analyzeDR
#combineEvidence
print("Number of candidates:%d"%(len(Candidates)),file=sys.stderr)
for C in Candidates:
    SV=callSV(ReferenceFile,C)
    if SV!="":
        print(SV)

ReferenceFile.close()