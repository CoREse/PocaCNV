import pysam
import sys
import time
from globals import *
from calling import *
from readpair import *
from readdepth import *
from filters import *
from benchmark_lib import *

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
ReadPairs=[]
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
                    ReadPairs.append(Pair)
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

if len(RDWindows)>1:
    RDWindows.append([0]*WindowsN)#last sample is the "sum sample"
    SumI=len(RDWindows)-1
    for i in range(WindowsN):
        for j in range(SumI):
            RDWindows[SumI][i]+=RDWindows[j][i]
for c in OccurredContigs:
    OccurredWindowsN+=ReferenceFile.lengths[c]/RDWindowSize+1
    ContigNameOccurred[ReferenceFile.references[c]]=True

EAFile=open("exclusion_regions.txt","r")
readExcludedAreas(EAFile,ContigNameOccurred)
EAFile.close()

print("Samples read, calculating RD data...", file=sys.stderr)
MixedRDRs=analyzeRD(RDWindows,WindowsN,OccurredWindowsN,True)

#read gold standard
print(OccurredContigs)
contigs=[ReferenceFile.references[key] for key in OccurredContigs.keys()]
print(contigs)
filename="/home/cre/data/HG00514.BIP-unified.vcf.gz"
vcf=pysam.VariantFile(filename,"r")
goldstandard=[]
included={}
if contigs!=None:
    for c in contigs:
        included[c]=True
for record in vcf.fetch():
    try:
        if contigs!=None:
            if not included[record.chrom]:
                continue
        temp=None
        temp=interval(record.chrom,RefStartPos[RefInd[record.chrom]]+record.pos,RefStartPos[RefInd[record.chrom]]+record.stop+1,record.info["SVTYPE"])
        goldstandard.append(temp)
    except:
        pass
vcf.close()
goldstandard.sort(key=lambda g:g.start)
'''
ITableS=[0]*int(RefLength/10000)+1
ITableE=[0]*int(RefLength/10000)+1
SI=0
EI=0
for i in range(len(ITableS)):
    if goldstandard[SI].start>i*10000:
        SI+=1
    if goldstandard[EI].end>i*10000:
        EI+=1
    ITableS[i]=SI-1
    ITableE[i]=EI
'''
print("size of gs:",len(goldstandard))

PairTraits=[]
traitfilename="bamtraits.txt"
traitfile=open(traitfilename,"w")

TotalPair=len(ReadPairs)
DoneCount=0
for p in ReadPairs:
    Trait={}
    Trait["LMapped"]=1 if p.hasFlags(PairInfo.LMapped) else 0
    Trait["LReversed"]=1 if p.hasFlags(PairInfo.LReversed) else 0
    Trait["RMapped"]=1 if p.hasFlags(PairInfo.RMapped) else 0
    Trait["RReversed"]=1 if p.hasFlags(PairInfo.RReversed) else 0
    Trait["TransChrs"]=1 if p.hasFlags(PairInfo.TransChrs) else 0
    Trait["Crossed"]=1 if p.hasFlags(PairInfo.Crossed) else 0
    Trait["Contained"]=1 if p.hasFlags(PairInfo.Contained) else 0
    Trait["LStart"]=p.LStart
    Trait["LEnd"]=p.LEnd
    Trait["RStart"]=p.RStart
    Trait["REnd"]=p.REnd
    Trait["NMapped"]=p.NMapped
    Trait["InsertionSize"]=p.InsertionSize
    Trait["LLength"]=p.LLength
    Trait["RLength"]=p.RLength
    Trait["Start"]=p.Start
    Trait["End"]=p.End
    ARD=0
    for i in range(int(p.Start/RDWindowSize),int(p.End/RDWindowSize)):
        ARD+=MixedRDRs[p.SampleIndex][i]
    ARD/=abs(int(p.End/RDWindowSize)-int(p.Start/RDWindowSize)) if int(p.End/RDWindowSize)-int(p.Start/RDWindowSize)!=0 else 1
    Trait["ARD"]=ARD
    pint=interval("",p.Start,p.End,"")
    Label="NONE"
    GSS=0#ITableE[p.Start/10000]
    GSE=len(goldstandard)#ITableS[p.End/10000]
    for i in range(GSS,GSE):
        gs=goldstandard[i]
        l=min(pint.end-pint.start,gs.end-gs.start)
        if pint.end-pint.start==0:
            Label="NONE"
            break
        if l==0:
            continue
        if gs.start>pint.end:
            break
        ol=calcOverlap(pint.start,pint.end,gs.start,gs.end)
        if ol/l>0.1:
            Label=gs.type
            break
    Trait["Label"]=Label
    pt=str(Trait)
    pt=pt.replace("'",'"')
    print(pt,file=traitfile)
    DoneCount+=1
    if DoneCount%10000==0:
        print("Readpairs:%s/%s"%(DoneCount,TotalPair),file=sys.stderr)
traitfile.close()

ReferenceFile.close()