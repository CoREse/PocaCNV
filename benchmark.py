from benchmark_lib import *
import sys

percentage=0.5
near=-1
Contigs=None
Contigs=["22"]
Samples=None
Samples=["HG00403"]
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
        elif a=="-PO":
            PrintOuts=True
        else:
            MyF=a
        i+=1
myout=parse_my(MyF,Contigs,Samples,format=Format)
goldstandard=parse_vcf(GSF,Contigs,Samples)
if PrintOuts:
    print("goldstandards:")
    for gs in goldstandard:
        print("%s"%gs)
    print("myouts:")
    for m in myout:
        print("%s"%m)
presult=calculate_sensitivity_by_percent(goldstandard,myout,percentage,near)
calculate_fdr_by_percent(goldstandard,myout,percentage,near)
sresult=calculate_sensitivity_by_inclusion(goldstandard,myout)
calculate_fdr_by_inclusion(goldstandard,myout)
if PrintResult:
    print("Results by percent:")
    for r in presult:
        print(r)
    print("Results by inclution:")
    for r in sresult:
        print(r)
