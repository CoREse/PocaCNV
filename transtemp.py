from benchmark_lib import *

'''
filename="/home/cre/data/HG00514.BIP-unified.vcf.gz"
vcf=pysam.VariantFile(filename,"r")
goldstandard=[]
included={}
SVTypes={}
for record in vcf.fetch():
    try:
        SVTypes[record.info["SVTYPE"]]=True
    except:
        pass
vcf.close()
print(len(SVTypes),SVTypes)

'''
af=open("bamtraits.txt","r")
bf=open("bamtraits.json","w")
for line in af:
    print(line.replace("'",'"'),end="",file=bf)
af.close()
bf.close()
