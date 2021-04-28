import sys
EDDataFileName=sys.argv[1]
BenchDataFileName=sys.argv[2]

EDDataFile=open(EDDataFileName,"r")
BenchDataFile=open(BenchDataFileName,"r")

FineResults=[]
for line in BenchDataFile:
    if line.strip()!="":
        sl=line.split("...")[1]
        if len(sl)<2:
            continue
        for s in sl.split("[")[1:]:
            SampleName=s.split("]")[0]
            BE=s.split("]")[1].split(",")[0].split(":")[1]
            Begin=int(BE.split("-")[0])-1
            End=int(BE.split("-")[1])-1
            FineResults.append((SampleName,Begin,End))

BenchDataFile.close()

def isFine(SampleName,Begin,End):
    for fs in FineResults:
        if SampleName==fs[0] and int(Begin)==int(fs[1]) and int(End)==int(fs[2]):
            return True
    return False

First=True
for line in EDDataFile:
    if line.strip()!="":
        sl=line.strip().split()
        SegNum=sl[0]
        ChromoLength=sl[1]
        SiblingCount=sl[2]
        SiblingRatio=sl[3]
        SampleName=sl[4]
        Begin=sl[5]
        End=sl[6]
        Mu=sl[7]
        MuS=sl[8]
        PassConfidence=sl[9]
        CN=sl[10]
        Confidence=sl[11]
        CScore=sl[12]
        HasSDRP=sl[13]
        HasMultiSDRP=sl[14]
        SDRPRatio=sl[15]
        ACL=sl[16]
        if not First:
            print("")
        First=False
        if isFine(SampleName,Begin,End):
            print("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"%(SegNum,ChromoLength,SiblingCount,SiblingRatio,Begin,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,HasSDRP,HasMultiSDRP,SDRPRatio,ACL,1),end="")
        else:
            print("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"%(SegNum,ChromoLength,SiblingCount,SiblingRatio,Begin,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,HasSDRP,HasMultiSDRP,SDRPRatio,ACL,0),end="")

EDDataFile.close()