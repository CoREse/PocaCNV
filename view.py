import sys
import pysam

fa=pysam.FastaFile(sys.argv[1])

sl=sys.argv[2].split(":")
if len(sl)>1:
    print(fa.fetch(sl[0],int(sl[1].split("-")[0]), int(sl[1].split("-")[1])))
else:
    print(fa.fetch(sys.argv[2],int(sys.argv[3]), int(sys.argv[4])))