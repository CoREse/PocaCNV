import pysam
import sys
import time
from globals import *
from calling import *
from readpair import *
from readdepth import *
from filters import *

WriteRDData=False

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
SampleNames=[]
for i in range(2,len(sys.argv)):
    if sys.argv[i].split(".")[-1]=="txt":
        readRDData(RDWindows,SampleNames,sys.argv[i])
        OccurredWindowsN=len(RDWindows[0])
    else:
        SamFile=pysam.AlignmentFile(sys.argv[i],"rb",reference_filename=sys.argv[1])
        SampleNames.append(sys.argv[i].split("/")[-1].split("\\")[-1])
        RDWindows.append([0]*WindowsN)
        for read in SamFile:
            ReadCount+=1
            read: pysam.AlignedSegment
            if read.is_unmapped:
                UnmappedCount+=1
            else:
                RDWindows[SampleIndex][int((int((read.reference_start+read.reference_end)/2)+RefStartPos[read.tid])/RDWindowSize)]+=1
                OccurredContigs[read.reference_id]=True
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
    OccurredWindowsN+=int(ReferenceFile.lengths[c]/RDWindowSize)+1
    ContigNameOccurred[ReferenceFile.references[c]]=True

EAFile=open("data/exclusion_regions.txt","r")
readExcludedAreas(EAFile,ContigNameOccurred)
EAFile.close()

if WriteRDData:
    writeRDData(RDWindows,ReferenceFile,SampleNames)
    exit(0)

print("Samples read, calculating RD data...", file=sys.stderr)
RDICandidates=analyzeRD(RDWindows,WindowsN,OccurredWindowsN)
RDICandidates=filtExcludedAreas(RDICandidates)
print("Number of RDI candidates:%d"%(len(RDICandidates)),file=sys.stderr)

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