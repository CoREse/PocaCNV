from benchlib import *
import sys


percentage=0.5
near=-1
Contigs=None
Contigs=None#["22"]
Samples=None
Samples=None#["HG00403"]
SamplesFile=None
PrintResult=True
PrintOuts=False
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/HG00514.BIP-unified.vcf.gz'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/ALL_Illumina_Integrate_20170206.vcf.gz'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/delly.bcf'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/nstd152.GRCh38.variant_call.vcf.gz'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/chr22_indels_HG00403.recode.vcf'
GSF='/mnt/c/Users/CRE/Productive/Programming/data/chr22_indels_CHS.recode.vcf'
#GSF='/mnt/c/Users/CRE/Productive/Programming/Workspace/jc/results/CHS_chr22_hs37d5_nc.vcf'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/delly_hs37d5_HG004xx.bcf'
MyF="test/512-4_chr22_testcall.txt"
Format="VCF"
if len(sys.argv)>1:
    i=1
    while i<len(sys.argv):
        a=sys.argv[i]
        if a=="-G":
            GSF=sys.argv[i+1]
            i+=1
        elif a=='-F':
            Format=sys.argv[i+1]
            i+=1
        elif a=='-C':
            Contigs=[sys.argv[i+1]]
            i+=1
        elif a=='-S':
            Samples=[sys.argv[i+1]]
            i+=1
        elif a=='-SF':
            SamplesFile=sys.argv[i+1]
            i+=1
        elif a=="-PO":
            PrintOuts=True
        else:
            MyF=a
        i+=1

if SamplesFile!=None:
    Samples=[]
    SF=open(SamplesFile,"r")
    for line in SF:
        if line.strip()!="":
            Samples.append(line.strip())
    SF.close()

Gold=VariantRecords("Gold")
Gold.parseVcfCNV(GSF,Contigs,Samples)
Test=VariantRecords("Test")
Test.parseVcfCNV(MyF,Contigs,Samples)

print(Test.interpret(Test.matchAll(Gold),PrintResult))