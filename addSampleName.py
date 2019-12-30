import sys
for a in sys.argv[1:]:
    SampleName=a.split("/")[-1].split("\\")[-1]
    if SampleName[:2]=="rd":
        SampleName=SampleName[2:]
    SampleName=SampleName.split(".")[0].split("_")[0].strip()
    F=open(a,"a")
    print("\n##Sample:%s"%SampleName,end="",file=F)
    F.close()