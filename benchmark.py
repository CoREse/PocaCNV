from benchmark_lib import *

percentage=0.01
near=100
goldstandard=parse_vcf('/mnt/c/Users/CRE/Productive/Programming/data/HG00514.BIP-unified.vcf.gz',["chr22"])
myout=parse_my("test/512-4_chr22_testcall.txt")
calculate_sensitivity_by_percent(goldstandard,myout,percentage,near)
calculate_fdr_by_percent(goldstandard,myout,percentage,near)