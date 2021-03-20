from rdread import *

import pysam
def writeRDData(mygenome,SampleNames,SampleReadCount):
    for i in range(len(SampleNames)):
        writeSampleRDData(mygenome,SampleNames[i],i,SampleReadCount[i],g.RDWindowSize,g.SamplePaths[i])
    return

def getRDFPath(SamPath,Path="data"):
    SamFileName=SamPath.split("/")[-1].split("\\")[-1]
    return "%s/%s.rdf"%(Path,SamFileName)

def writeSampleRDData(mygenome, SampleName, SampleI, SampleReadCount, WindowSize, SamplePath, Path="data"):
    #SamFileName=g.SamplePaths[SampleI].split("/")[-1].split("\\")[-1]
    rdfile=open(getRDFPath(SamplePath,Path),"w")
    print("##Sample:%s"%SampleName,end="",file=rdfile)
    print("\n##SampleReadCount:%s\n##WindowSize:%s"%(SampleReadCount,WindowSize),end="",file=rdfile)
    for c in mygenome.Contigs:
        print("\n#%s %s"%(c.Name,c.Length),end="",file=rdfile)
        for j in range(len(c.RDWindows[SampleI])):
            print("\n%d"%(c.RDWindows[SampleI][j]),end="",file=rdfile)
    rdfile.close()
    return

def writeMixedRDData(mygenome,SampleNames):
    #ReferenceFile:pysam.FastaFile
    for i in range(len(SampleNames)):
        SamFileName=g.SamplePaths[i].split("/")[-1].split("\\")[-1]
        rdfile=open("data/%s.mrd"%(SamFileName),"w")
        print("##Sample:%s"%SampleNames[i],end="",file=rdfile)
        for c in mygenome.Contigs:
            print("\n#%s %s"%(c.Name,c.Length),end="",file=rdfile)
            for j in range(len(c.MixedRDRs[i])):
                print("\n%.8s"%(c.MixedRDRs[i][j]),end="",file=rdfile)
        rdfile.close()
    return

def readRDDataNA0(ContigsWindows, ContigSet, SampleI, FileName):#read, no adding sample
    print(gettime()+"Loading from %s..."%(FileName),file=sys.stderr)
    SampleName=FileName.split("\\")[-1].split("/")[-1][:-4]
    if "rd"==SampleName[:2]:
        SampleNameS=SampleName.split("rd")
        SampleName=""
        for s in SampleNameS[1:]:
            SampleName+=s
    DataFile=open(FileName,"r")
    ContigName=None
    Skip=False
    ConI=0
    ReadCount=0
    SampleReadCount=None
    Windows=None
    ContigSampleReadCount={}
    for line in DataFile:
        if line[0]=='#':
            if line[1]=="#":
                sl=line[2:].split(":",maxsplit=1)
                if sl[0]=="Sample":
                    SampleName=sl[1].strip()
                elif sl[0]=="SampleReadCount":
                    SampleReadCount=int(sl[1].strip())
            else:
                sl=line[1:].split()
                ContigName=sl[0]
                if ContigName not in ContigSet:
                    Skip=True
                    continue
                Length=int(sl[1])
                if ContigName not in ContigSampleReadCount.keys():
                    ContigSampleReadCount[ContigName]=0
                Windows=ContigsWindows[ContigName]
                Skip=False
                ConI=0
            continue
        if Skip:
            continue
        #try:
        #    Windows[ConI]=int(line)
        #except IndexError as e:
        #    if len(Windows)<=ConI:
        #        print("data exceed contig %s's capacity(data no.%s, len of %s:%s)!"%(ContigName,ConI,ContigName,len(Windows)),file=sys.stderr)
        #    else:
        #        raise e
        Windows[ConI]=float(line)
        ReadCount+=Windows[ConI]
        ContigSampleReadCount[ContigName]+=Windows[ConI]
        ConI+=1
    if SampleReadCount==None:
        RC=ReadCount
    else:
        RC=SampleReadCount
    DataFile.close()
    return SampleName, RC, ContigSampleReadCount

def readRDDataNA(SampleI, FileName):#read, no adding sample
    print(gettime()+"Loading from %s..."%(FileName),file=sys.stderr)
    mygenome=g.MyGenome
    SampleName=FileName.split("\\")[-1].split("/")[-1][:-4]
    if "rd"==SampleName[:2]:
        SampleNameS=SampleName.split("rd")
        SampleName=""
        for s in SampleNameS[1:]:
            SampleName+=s
    DataFile=open(FileName,"r")
    ContigName=None
    Skip=False
    Windows=None
    ConI=0
    ReadCount=0
    SampleReadCount=None
    ContigSampleReadCount={}
    for line in DataFile:
        if line[0]=='#':
            if line[1]=="#":
                sl=line[2:].split(":",maxsplit=1)
                if sl[0]=="Sample":
                    SampleName=sl[1].strip()
                elif sl[0]=="SampleReadCount":
                    SampleReadCount=int(sl[1].strip())
            else:
                sl=line[1:].split()
                ContigName=sl[0]
                if not mygenome.hasContig(ContigName):
                    Skip=True
                    continue
                Length=int(sl[1])
                TheContig=mygenome.get(ContigName)
                Windows=TheContig.RDWindows[SampleI]
                if ContigName not in ContigSampleReadCount.keys():
                    ContigSampleReadCount[ContigName]=0
                Skip=False
                ConI=0
            continue
        if Skip:
            continue
        #try:
        #    Windows[ConI]=int(line)
        #except IndexError as e:
        #    if len(Windows)<=ConI:
        #        print("data exceed contig %s's capacity(data no.%s, len of %s:%s)!"%(ContigName,ConI,ContigName,len(Windows)),file=sys.stderr)
        #    else:
        #        raise e
        Windows[ConI]=float(line)
        ReadCount+=Windows[ConI]
        ContigSampleReadCount[ContigName]+=Windows[ConI]
        ConI+=1
    if SampleReadCount==None:
        RC=ReadCount
    else:
        RC=SampleReadCount
    DataFile.close()
    return SampleName, RC, ContigSampleReadCount

def readRDDataAll(SampleNames, FileNames):
    mygenome=g.MyGenome
    for i in FileNames:
        mygenome.addSample("")
    #ContigsWindows={}
    #ContigSet=set()
    #for i in range(len(mygenome.Contigs)):
    #    ContigSet.add(mygenome.Contigs[i].Name)
    if g.ThreadN==1:
        for i in range(len(FileNames)):
            #for j in range(len(mygenome.Contigs)):
            #    ContigsWindows[mygenome.Contigs[j].Name]=mygenome.Contigs[j].RDWindows[i]
            SampleName, RC, CRC=readRDDataNA(i,FileNames[i])
            mygenome.changeSampleName(i,SampleName)
            SampleNames.append(SampleName)
            g.SampleReadCount.append(RC)
            for j in range(len(mygenome.Contigs)):
                mygenome.Contigs[j].ContigSampleReadCounts[i]=CRC[mygenome.Contigs[j].Name]
        print(gettime()+"RDFs read. "+getMemUsage(),file=sys.stderr)
    else:
        ctx=mp.get_context("fork")
        pool=ctx.Pool(g.ThreadN)
        args=[]
        for i in range(len(FileNames)):
            #for j in range(len(mygenome.Contigs)):
            #    ContigsWindows[mygenome.Contigs[j].Name]=mygenome.Contigs[j].RDWindows[i]
            args.append((i,FileNames[i]))
        addPool(pool)
        Results=pool.starmap(readRDDataNA,args)
        print(gettime()+"RDFs read. "+getMemUsage(),file=sys.stderr)
        delPool()
        pool.terminate()
        for i in range(len(Results)):
            SampleNames.append(Results[i][0])
            g.SampleReadCount.append(Results[i][1])
            for j in range(len(mygenome.Contigs)):
                mygenome.Contigs[j].ContigSampleReadCounts[i]=Results[i][2][mygenome.Contigs[j].Name]

def readRDData(mygenome, SampleNames, FileName):
    SampleName=FileName.split("\\")[-1].split("/")[-1][:-4]
    if "rd"==SampleName[:2]:
        SampleNameS=SampleName.split("rd")
        SampleName=""
        for s in SampleNameS[1:]:
            SampleName+=s
    DataFile=open(FileName,"r")
    ContigName=None
    mygenome.addSample(SampleName)
    Skip=False
    Windows=None
    ConI=0
    ReadCount=0
    SampleReadCount=None
    for line in DataFile:
        if line[0]=='#':
            if line[1]=="#":
                sl=line[2:].split(":",maxsplit=1)
                if sl[0]=="Sample":
                    SampleName=sl[1].strip()
                elif sl[0]=="SampleReadCount":
                    SampleReadCount=int(sl[1].strip())
            else:
                sl=line[1:].split()
                ContigName=sl[0]
                if not mygenome.hasContig(ContigName):
                    Skip=True
                    continue
                Length=int(sl[1])
                TheContig=mygenome.get(ContigName)
                Windows=TheContig.RDWindows[-1]
                ContigSampleReadCounts=TheContig.ContigSampleReadCounts
                Skip=False
                ConI=0
            continue
        if Skip:
            continue
        #try:
        #    Windows[ConI]=int(line)
        #except IndexError as e:
        #    if len(Windows)<=ConI:
        #        print("data exceed contig %s's capacity(data no.%s, len of %s:%s)!"%(ContigName,ConI,ContigName,len(Windows)),file=sys.stderr)
        #    else:
        #        raise e
        Windows[ConI]=float(line)
        ReadCount+=Windows[ConI]
        ContigSampleReadCounts[-1]+=Windows[ConI]
        ConI+=1
    SampleNames.append(SampleName)
    mygenome.changeSampleName(-1,SampleName)
    if SampleReadCount==None:
        g.SampleReadCount.append(ReadCount)
    else:
        g.SampleReadCount.append(SampleReadCount)
    DataFile.close()
    return