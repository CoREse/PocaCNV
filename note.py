import pysam
import sys
SamFile=pysam.AlignmentFile(sys.argv[1],"rb",reference_filename=sys.argv[2])
n=0
UnpairedReads={}
DiscordantReads=[]
i=0
for read in SamFile:
    read: pysam.AlignedSegment
    if read.is_paired:
        if UnpairedReads.count(read.query_name)==0:
            UnpairedReads[read.query_name]=read
        elif UnpairedReads.count(read.query_name)==1:
            Pair=analysePair(UnpairedReads.pop(read.query_name),read)
            if Pair!=None:
                DiscordantReads.append(Pair)
SamFile.close()
FastaFile=pysam.FastaFile(sys.argv[2])
#print(FastaFile.references)
FastaFile.close()
print(n)