from benchmark_lib import *
import sys

percentage=0.5
near=-1
Contigs=None
#Contigs=["chr22"]
Samples=None
Samples=["HG00514"]
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/HG00514.BIP-unified.vcf.gz'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/ALL_Illumina_Integrate_20170206.vcf.gz'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/delly.bcf'
GSF='/mnt/c/Users/CRE/Productive/Programming/data/nstd152.GRCh38.variant_call.vcf.gz'
goldstandard=parse_vcf(GSF,Contigs,Samples)
if len(sys.argv)>1:
    myout=parse_my(sys.argv[1],Contigs,Samples)
else:
    myout=parse_my("test/512-4_chr22_testcall.txt")
calculate_sensitivity_by_percent(goldstandard,myout,percentage,near)
calculate_fdr_by_percent(goldstandard,myout,percentage,near)
calculate_sensitivity_by_inclusion(goldstandard,myout)
calculate_fdr_by_inclusion(goldstandard,myout)