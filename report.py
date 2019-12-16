import globals as g
from consts import DefaultHeader
def reportVCF(SVs,RefSeq,OutFile,ReportHeader=False):
    if ReportHeader:
        reportVCFHeader(OutFile)
    for SV in SVs:
        ##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
        INFO=""
        Ref=RefSeq[SV.BreakLeft-1]
        Alt=""
        fal=True
        for a in SV.Alleles:
            if not fal:
                Alt+=","
            fal=False
            Alt+=a
        INFO="SVTYPE=%s,End=%s"%(SV.Type,SV.BreakRight)
        print("%s %s * %s %s 100 PASS %s GT"%(SV.Chrom,SV.BreakLeft,Ref,Alt,INFO),end="",file=OutFile)
        SI=0
        for i in range(len(g.SampleNames)):
            if SI>=len(SV.Samples) or i!=SV.Samples[SI][0]:
                print(" 0/0",end="",file=OutFile)
            else:
                print(" %s/%s"%(SV.Samples[SI][1][0],SV.Samples[SI][1][1]),end="",file=OutFile)
                SI+=1
        print("\n",end="",file=OutFile)

def reportVCFHeader(OutFile):
    Header=DefaultHeader
    for s in g.SampleNames:
        Header+=" "+s
    print(Header,file=OutFile)