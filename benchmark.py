from benchmark_lib import *
import sys

percentage=0.5
near=-1
goldstandard=parse_vcf('/mnt/c/Users/CRE/Productive/Programming/data/HG00514.BIP-unified.vcf.gz',["chr22"])
if len(sys.argv)>1:
    myout=parse_my(sys.argv[1])
else:
    myout=parse_my("test/512-4_chr22_testcall.txt")
calculate_sensitivity_by_percent(goldstandard,myout,percentage,near)
calculate_fdr_by_percent(goldstandard,myout,percentage,near)
calculate_sensitivity_by_inclusion(goldstandard,myout)
calculate_fdr_by_inclusion(goldstandard,myout)