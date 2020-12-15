import sys
EDDataFileName=sys.argv[1]
BenchDataFileName=sys.argv[2]

EDDataFile=open(EDDataFileName,"r")
BenchDataFile=open(BenchDataFileName,"r")

FineResults=[]
for line in BenchDataFile:
    if line.strip()!="":
        sl=line.split("...")[1]
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
        SampleName=sl[3]
        Begin=sl[4]
        End=sl[5]
        Mu=sl[6]
        MuS=sl[7]
        PassConfidence=sl[8]
        CN=sl[9]
        Confidence=sl[10]
        CScore=sl[11]
        if not First:
            print("")
        First=False
        if isFine(SampleName,Begin,End):
            print("%s %s %s %s %s %s %s %s %s %s %s %s"%(SegNum,ChromoLength,SiblingCount,Begin,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,1),end="")
        else:
            print("%s %s %s %s %s %s %s %s %s %s %s %s"%(SegNum,ChromoLength,SiblingCount,Begin,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,0),end="")

EDDataFile.close()