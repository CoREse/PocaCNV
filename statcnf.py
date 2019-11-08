import pysam
import sys

NN=11
CNPriors=[0]*NN

for filename in sys.argv[1:]:
    vcf=pysam.VariantFile(filename,"r")
    for record in vcf.fetch():
        SVLen=0
        OccuredCN=[]
        try:
            if record.info["SVTYPE"]=="CNV":
                SVLen=record.rlen
                for v in record.alts:
                    if v[:3]=="<CN":
                        CN=int(v[3:-1])
                        OccuredCN.append(CN)
                for i in range(len(record.info["AF"])):
                    if OccuredCN[i]<NN:
                        CNPriors[OccuredCN[i]]+=SVLen*record.info["AF"][i]
            else:
                continue
        except:
            pass
    vcf.close()

for i in range(NN):
    CNPriors[i]/=3e9
CNPriors[1]=1
for i in range(NN):
    if i==1:
        continue
    CNPriors[1]-=CNPriors[i]
CNPriorsD=[0]*NN
for i in range(NN):
    for j in range(i+1):
        CNPriorsD[i]+=CNPriors[j]*CNPriors[i-j]
print(CNPriorsD)