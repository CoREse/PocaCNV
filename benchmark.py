from benchmark_lib import *

goldstandard=parse_vcf('/mnt/c/Users/CRE/Productive/Programming/data/HG00514.BIP-unified.vcf.gz',["chr1"])
myout=parse_my("512-4testcall.txt")
calculate_sensitivity_by_percent(goldstandard,myout,0.5)
calculate_fdr_by_percent(goldstandard,myout,0.5)