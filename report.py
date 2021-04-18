import globals as g
from utils import gettime
from consts import DefaultHeader,CNPriors
def reportVCF(SVs,RefSeq,OutFile,ReportHeader=False,MyGenome=None):
    if ReportHeader:
        reportVCFHeader(OutFile,MyGenome)
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
        SVLEN=SV.BreakRight-SV.BreakLeft
        NS=len(g.SampleNames)
        AN=2*NS
        AC=[0]*len(SV.Alleles)
        for s in SV.Samples:
            if s[1][0]!=0:
                AC[s[1][0]-1]+=1
            if s[1][1]!=0:
                AC[s[1][1]-1]+=1
        AF=[]
        for a in AC:
            AF.append(a/AN)
        first=True
        ACS=""
        AFS=""
        for i in range(len(AC)):
            if not first:
                ACS+=","
                AFS+=","
            first=False
            ACS+="%d"%AC[i]
            AFS+="%s"%AF[i]
        SC="%s"%SV.Score
        INFO="SVTYPE=%s;END=%s;IMPRECISE;SVLEN=%d;AC=%s;AF=%s;NS=%d;AN=%d;SC=%s;VT=SV"%(SV.SVType,SV.BreakRight,SVLEN,ACS,AFS,NS,AN,SC)
        print("%s\t%s\t*\t%s\t%s\t100\tPASS\t%s\tGT:CS"%(SV.Chrom,SV.BreakLeft,Ref,Alt,INFO),end="",file=OutFile)
        SI=0
        for i in range(len(g.SampleNames)):
            if SI>=len(SV.Samples) or i!=SV.Samples[SI][0]:
                print("\t0/0:.",end="",file=OutFile)
            else:
                print("\t%s/%s:%s"%(SV.Samples[SI][1][0],SV.Samples[SI][1][1],SV.Samples[SI][2]),end="",file=OutFile)
                SI+=1
        print("\n",end="",file=OutFile)

def reportVCFHeader(OutFile,MyGenome):
    Header="""##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate="""+gettime()
    Header+="\n##reference=%s"%g.ReferencePath
    Header+="\n##source=jcrd#"
    first=True
    for k in vars(g.Parameters):
        if k[:2]!="__":
            if not first:
                Header+=";"
            first=False
            Header+="%s:%s"%(k,vars(g.Parameters)[k])
    for c in MyGenome:
        Header+="\n##contig=<ID=%s,assembly=b37,length=%d>"%(c.Name,c.Length)
    Header+="""
##ALT=<ID=CNV,Description="Copy Number Polymorphism">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">"""
    for i in range(len(CNPriors)):
        Header+='\n##ALT=<ID=CN%d,Description="Copy number allele: %d copies">'%(i,i)
    Header+="""
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CS,Number=1,Type=Float,Description="Confidence Score">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length. It is only calculated for structural variation MEIs. For other types of SVs; one may calculate the SV length by INFO:END-START+1, or by finding the difference between lengthes of REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=SC,Number=1,Type=Float,Description="Calling score of this CNV">
##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""
    for s in g.SampleNames:
        Header+="\t"+s
    print(Header,file=OutFile)

def reportVcfHeaderGeneral(OutFile,ReferencePath,Parameters,Contigs,SampleNames):
    Header="""##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate="""+gettime()
    Header+="\n##reference=%s"%ReferencePath
    Header+="\n##source=jcrd#"
    first=True
    for k in vars(Parameters):
        if k[:2]!="__":
            if not first:
                Header+=";"
            first=False
            Header+="%s:%s"%(k,vars(Parameters)[k])
    for c in Contigs:
        Header+="\n##contig=<ID=%s,assembly=b37,length=%d>"%(c[0],c[1])
    Header+="""
##ALT=<ID=CNV,Description="Copy Number Polymorphism">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">"""
    for i in range(len(CNPriors)):
        Header+='\n##ALT=<ID=CN%d,Description="Copy number allele: %d copies">'%(i,i)
    Header+="""
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CS,Number=1,Type=Float,Description="Confidence Score">
##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="SV length. It is only calculated for structural variation MEIs. For other types of SVs; one may calculate the SV length by INFO:END-START+1, or by finding the difference between lengthes of REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=SC,Number=1,Type=Float,Description="Calling score of this CNV">
##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""
    for s in SampleNames:
        Header+="\t"+s
    print(Header,file=OutFile)