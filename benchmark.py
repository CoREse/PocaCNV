from benchlib import *
import sys


Percentage=0.8
near=-1
Contigs=None
Contigs=None#["22"]
Samples=None
Samples=None#["HG00403"]
SamplesFile=None#one sample name per line
PrintResult=False
PrintOuts=False
PrintSample=None
MinLength=0#Minimum Variant Length
MinScore=0
MaxAF=1#Max allele frequency, only apply for gold standard
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/HG00514.BIP-unified.vcf.gz'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/ALL_Illumina_Integrate_20170206.vcf.gz'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/delly.bcf'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/nstd152.GRCh38.variant_call.vcf.gz'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/chr22_indels_HG00403.recode.vcf'
GSF='/home/cre/data/chr22_indels_CHS.recode.vcf'
#GSF='/mnt/c/Users/CRE/Productive/Programming/Workspace/jc/results/CHS_chr22_hs37d5_nc.vcf'
#GSF='/mnt/c/Users/CRE/Productive/Programming/data/delly_hs37d5_HG004xx.bcf'
MyF=[]#"test/512-4_chr22_testcall.txt"

Merge=False
MG=""

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
            if Contigs==None:
                Contigs=[sys.argv[i+1]]
            else:
                Contigs.append(sys.argv[i+1])
            i+=1
        elif a=='-S':
            Samples=[sys.argv[i+1]]
            i+=1
        elif a=='-SF':
            SamplesFile=sys.argv[i+1]
            i+=1
        elif a=='-ML':
            MinLength=int(sys.argv[i+1])
            i+=1
        elif a=='-MAF':
            MaxAF=float(sys.argv[i+1])
            i+=1
        elif a=="-PO":
            PrintOuts=True
        elif a=="-PS":
            PrintSample=sys.argv[i+1]
            i+=1
        elif a=="-PD":
            PrintResult=True
        elif a=="-MSC":
            MinScore=float(sys.argv[i+1])
            i+=1
        elif a=="-P":
            Percentage=float(sys.argv[i+1])
            i+=1
        elif a=="-M":
            Merge=True
        elif a=="-MG":#merge groups, 10,1,3 to indicate number of vcfs of each group
            MG=sys.argv[i+1]
            i+=1
        else:
            MyF.append(a)
        i+=1

if SamplesFile!=None:
    Samples=[]
    SF=open(SamplesFile,"r")
    for line in SF:
        if line.strip()!="":
            Samples.append(line.strip())
    SF.close()

Gold=VariantRecords("Gold")
Gold.parseVcfCNV(GSF,Contigs,Samples,MinLength,MaxAF)
if PrintSample!=None:
    Gold.printSample(PrintSample)
    exit(0)

if (not Merge) or len(MyF)==1:
    Test=VariantRecords("Test")
    for i in range(len(MyF)):
        try:
            Test.parseVcfCNV(MyF[i],Contigs,Samples,MinLength,1,MinScore)
        except Exception as e:
            print(e,file=sys.stderr)

    print(Test.interpret(Test.matchAll(Gold,Percentage),PrintResult))
else:
    #Merge rule: same type overlap>=MergeThreshold or match to same gs
    MergeThreshold=0.95
    if MG!="":
        MG=MG.split(",")
        for i in range(len(MG)):
            MG[i]=int(MG[i])
    else:
        for i in range(len(MyF)):
            MG.append(1)  
    Tests=[]
    TResults=[]
    Matched=[]
    MyFL=0
    for i in range(len(MG)):
        Tests.append(VariantRecords("Test%s"%i))
        Test=Tests[-1]
        for j in range(MyFL, MyFL+MG[i]):
            try:
                Test.parseVcfCNV(MyF[j],Contigs,Samples,MinLength,1,MinScore)
            except Exception as e:
                print(e,file=sys.stderr)
        MyFL+=MG[i]
        TResults.append(Test.matchAll(Gold,Percentage))
        Matched.append([])
        for ss in TResults[-1]["Samples"]:
            if TResults[-1]["Samples"][ss]==None:
                continue
            for m in TResults[-1]["Samples"][ss]:
                Matched[-1].append(m[1])
        Tests[i].Variants.sort()
    
    MergedN=[0,0,0,0]#Merged Del, Merged Dup, Merged Matched DEL, merged Matched DUP
    SameMergeN=0
    for i in range(len(MG)):
        Test=Tests[i]
        j=0
        while j< len(Test.Variants):
            v=Test.Variants[j]
            Merged=False
            for k in range(len(MG)):
                if k==i:
                    continue
                m=0
                while m < len(Tests[k].Variants):
                    x=Tests[k].Variants[m]
                    if x.Start>v.End:
                        break
                    if x.Type!=v.Type or x.Sample!=v.Sample:
                        m+=1
                        continue
                    Over=calcOverlap(v.Start,v.End,x.Start,x.End)
                    if (Over/(v.End-v.Start)>=MergeThreshold and Over/(x.End-x.Start)>=MergeThreshold) or (v.MatchTo!=None and v.MatchTo==x.MatchTo):
                        if x.MatchTo!=None:
                            MergedN[2+(0 if v.Type=="DEL" else 1)]+=1
                            if x.MatchTo==v.MatchTo:
                                SameMergeN+=1
                        if v.MatchTo!=None and v.MatchTo==x.MatchTo:
                            try:
                                Matched[i].remove(v.MatchTo)
                                Matched[k].remove(x.MatchTo)
                            except ValueError:
                                pass
                        del Tests[k].Variants[m]
                        Merged=True
                        MergedN[0 if v.Type=="DEL" else 1]+=1
                        continue
                    m+=1
            if False and (not Merged) and v.MatchTo!=None:
                for k in range(len(MG)):
                    if k==i:
                        continue
                    if v.MatchTo in Matched[k]:
                        Matched[k].remove(v.MatchTo)
                        del Test.Variants[j]
                        MergedN[0 if v.Type=="DEL" else 1]+=1
                        MergedN[2+(0 if v.Type=="DEL" else 1)]+=1
                        break
            j+=1
    print(MergedN,SameMergeN,file=sys.stderr)
    Test=VariantRecords("Test")
    for i in range(len(MyF)):
        try:
            Test.parseVcfCNV(MyF[i],Contigs,Samples,MinLength,1,MinScore)
        except Exception as e:
            print(e,file=sys.stderr)

    print(Test.interpret(Test.matchAll(Gold,Percentage),PrintResult,MergedN))
        
    