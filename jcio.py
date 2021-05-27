from rdread import *

import pysam
import pickle
from DRP import *
import h5py
import array

def readSamAndSaveSD(thegenome, SamplePaths):
    RDFPaths=SamplePaths.copy()
    SAMPaths=[]
    RtS=[]
    for i in range(len(SamplePaths)):
        #print(gettime()+"Reading %s..."%(g.SamplePaths[i]),file=sys.stderr)
        RtS.append(-1)
        if SamplePaths[i].split(".")[-1]=="sd" or SamplePaths[i].split(".")[-1]=="hdf5":
            continue
        else:
            RDFPaths[i]=getSDPath(SamplePaths[i])
            SAMPaths.append(SamplePaths[i])
            RtS[-1]=len(SAMPaths)-1
    if len(SAMPaths)>0:
        print(gettime()+"Start reading SAM files and save them to SD files for following analysis...",file=sys.stderr)
        if g.ThreadNR>1:
            ctx=mp.get_context()
            pool=ctx.Pool(g.ThreadNR)
            args=[]
            for i in range(len(SAMPaths)):
                args.append((thegenome.genVacant(),SAMPaths[i],g.ReferencePath,g.RDWindowSize))
            addPool(pool)
            Results=pool.starmap(readSamToSD,args)
            delPool()
            pool.terminate()
        else:
            Results=[]
            for i in range(len(SAMPaths)):
                Results.append(readSamToSD(thegenome.genVacant(),SAMPaths[i],g.ReferencePath,g.RDWindowSize))
        NoProblem=True
        RDFPathsO=RDFPaths
        RDFPaths=[]
        for i in range(len(Results)):
            if Results[i]==False:
                NoProblem=False
        for i in range(len(RDFPathsO)):
            if RtS[i]!=-1:
                if Results[RtS[i]]:
                    RDFPaths.append(RDFPathsO[i])
            else:
                RDFPaths.append(RDFPathsO[i])
        if NoProblem:
            print(gettime()+"All SAMFile(s) read and SDs made. "+getMemUsage(),file=sys.stderr)
        else:
            print(gettime()+"A part of SAMFile(s) read and SDs made. "+getMemUsage(),file=sys.stderr)
    return RDFPaths

def getClipLength(Cigar):
    Length=0
    for i in range(len(Cigar)):
        if Cigar[i][0]==4 or Cigar[i][0]==5:
            Length+=Cigar[i][1]
    return Length

def readSamToSD(thegenome, FilePath, ReferencePath, WindowSize):
    print(gettime()+"Reading reads from %s..."%(FilePath),file=sys.stderr)
    SamFile=pysam.AlignmentFile(FilePath,"rb",reference_filename=ReferencePath)
    Name=SamFile.header["RG"][0]["SM"]
    thegenome.addSample(Name)
    ReadCount=0
    SampleIndex=0
    UnmappedCount=0
    LastRefId=-1
    try:
        for read in SamFile:
            ReadCount+=1
            #read: pysam.AlignedSegment
            if read.is_unmapped:
                UnmappedCount+=1
            else:
                if LastRefId!=read.reference_id:
                    TheContig=thegenome.get(read.reference_name)#This is slightly faster than using getcontigid, and is correct at all time
                    LastRefId=read.reference_id
                #TheContig=thegenome.get(read.reference_name)#This is slightly faster than using getcontigid, and is correct at all time
                if TheContig==None:
                    continue
                if read.is_paired and read.is_read1 and (not read.mate_is_unmapped):
                    isize=abs(read.template_length)
                    if isize > g.MedianInsertionSize+3*g.ISSD:
                        TheContig.DRPs[SampleIndex].append(DRP(read.reference_start,read.reference_start+read.template_length,read.reference_end,read.next_reference_start,0))
                    elif isize < g.MedianInsertionSize-3*g.ISSD:
                        TheContig.DRPs[SampleIndex].append(DRP(read.reference_start,read.reference_start+read.template_length,read.reference_end,read.next_reference_start,1))
                    if g.AdaptiveIS:
                        g.MedianInsertionSize=(g.MedianInsertionSize*g.AdaptiveCount+isize)/(g.AdaptiveCount+1)
                        g.ISSD=(g.ISSD*g.AdaptiveCount+abs(isize-g.MedianInsertionSize))/(g.AdaptiveCount+1)
                        g.AdaptiveCount+=1
                TheContig.RDWindows[SampleIndex][int((int((read.reference_start+read.reference_end)/2))/WindowSize)]+=1
                ClipLength=getClipLength(read.cigar)
                TheContig.AverageClipLengths[SampleIndex][int((int((read.reference_start+read.reference_end)/2))/WindowSize)]+=ClipLength
                TheContig.ContigSampleReadCounts[SampleIndex]+=1
                #ReadCount+=1
    except Exception as e:
        print("WARN:[Sample:%s] %s"%(Name,e),file=sys.stderr)
        return False
    SamFile.close()
    for TheContig in thegenome.Contigs:
        TheContig.DRPs[SampleIndex].sort(key=lambda d: d.Start)
        for i in range(len(TheContig.AverageClipLengths[SampleIndex])):
            if TheContig.RDWindows[SampleIndex][i]!=0:
                TheContig.AverageClipLengths[SampleIndex][i]/=TheContig.RDWindows[SampleIndex][i]
    print(gettime()+"Sample %s read."%(Name),file=sys.stderr)
    print(gettime()+"Storing sd data for %s..."%(Name),file=sys.stderr)
    writeSampleSDData(thegenome,Name,SampleIndex,ReadCount,WindowSize,FilePath)
    print(gettime()+"Stored sd data for %s."%(Name),file=sys.stderr)
    return True

def writeSampleSDData(mygenome, SampleName, SampleI, SampleReadCount, WindowSize, SamplePath, Path="data"):
    sdfile=open(getSDPath(SamplePath,Path),"wb")
    SampleData={}
    SampleData["SampleName"]=SampleName
    SampleData["SampleReadCount"]=SampleReadCount
    SampleData["RDWindowSize"]=WindowSize
    SampleData["Contigs"]=[]
    for c in mygenome.Contigs:
        SampleData["Contigs"].append({"Name":c.Name,"Length":c.Length,"RDWindows":c.RDWindows[SampleI],"AverageClipLengths":c.AverageClipLengths[SampleI],"ContigSampleReadCounts":c.ContigSampleReadCounts[SampleI],"DRPs":c.DRPs[SampleI]})
    pickle.dump(SampleData,sdfile)
    sdfile.close()
    return

def readSampleSDData(mygenome, FileName):
    print(gettime()+"Reading reads from %s..."%(FileName),file=sys.stderr)
    SDFile=open(FileName,"rb")
    try:
        SampleData=pickle.load(SDFile)
    except Exception as e:
        print("[ERROR]:%s, FileName:%s"%(e,FileName),file=sys.stderr)
        raise(e)
    SDFile.close()

    SampleName=SampleData["SampleName"]
    mygenome.addSample(SampleName)
    if g.RDWindowSize!=SampleData["RDWindowSize"]:
        print(gettime()+"Wrong RDWindowSize of SD file %s, quitting..."%FileName,file=sys.stderr)
        exit(0)
    for c in SampleData["Contigs"]:
        if not mygenome.hasContig(c["Name"]):
            continue
        TheContig=mygenome.get(c["Name"])
        TheContig.RDWindows[-1]=c["RDWindows"]
        TheContig.AverageClipLengths[-1]=c["AverageClipLengths"]
        if TheContig.Length!=c["Length"]:
            print(gettime()+"WARN: contig %s length doesn't matchup."%c["Name"],file=sys.stderr)
        TheContig.DRPs[-1]=c["DRPs"]
        TheContig.ContigSampleReadCounts[-1]=c["ContigSampleReadCounts"]
    return (SampleName,SampleData["SampleReadCount"],mygenome)

def readSDDataAll(TheGenome, SampleNames, FileNames):
    mygenome=TheGenome
    if g.ThreadNR==1:
        for i in range(len(FileNames)):
            Result=readSampleSDData(mygenome,FileNames[i])
            g.SampleReadCount.append(Result[1])
            SampleNames.append(Result[0])
            g.SampleNameIndexes[Result[0]]=i
        print(gettime()+"SDs read. "+getMemUsage(),file=sys.stderr)
    else:
        ctx=mp.get_context("fork")
        pool=ctx.Pool(g.ThreadNR)
        args=[]
        for i in range(len(FileNames)):
            args.append((mygenome.genVacant(),FileNames[i]))
        addPool(pool)
        Results=pool.starmap(readSampleSDData,args)
        delPool()
        pool.terminate()
        for j in range(len(Results)):
            r=Results[j]
            SampleName=r[0]
            g.SampleReadCount.append(r[1])
            SampleNames.append(SampleName)
            g.SampleNameIndexes[SampleName]=j
            r=r[2]
            for i in range(len(mygenome.Contigs)):
                c=mygenome.Contigs[i]
                SampleRDWindows=r.Contigs[i].RDWindows[0]
                SampleAverageClipLengths=r.Contigs[i].AverageClipLengths[0]
                ContigSampleReadCount=r.Contigs[i].ContigSampleReadCounts[0]
                SampleDRPs=r.Contigs[i].DRPs[0]
                Ploidy=r.Contigs[i].Ploidies[0]
                c.addSampleWithData(SampleName,SampleRDWindows,SampleAverageClipLengths,ContigSampleReadCount,SampleDRPs,Ploidy)
            mygenome.SampleNames.append(SampleName)
            mygenome.SampleN+=1

def readSampleHDF5Data(mygenome, FileName):
    print(gettime()+"Reading reads from %s..."%(FileName),file=sys.stderr)
    try:
        HDF5File=h5py.File(FileName,"r")
    except Exception as e:
        print("[ERROR]:%s, FileName:%s"%(e,FileName),file=sys.stderr)
        raise(e)
    
    SampleName=HDF5File.attrs["SampleName"]
    if type(SampleName)!=str:
        SampleName=str(SampleName.decode())
    SampleReadCount=HDF5File.attrs["SampleReadCount"]
    mygenome.addSample(SampleName)
    if g.RDWindowSize!=HDF5File.attrs["RDWindowSize"]:
        print(gettime()+"Wrong RDWindowSize of Sample HDF5 file %s, quitting..."%FileName,file=sys.stderr)
        exit(0)
    for c in HDF5File.keys():
        if not mygenome.hasContig(c):
            continue
        TheContig=mygenome.get(c)
        ContigG=HDF5File[c]
        TheContig.RDWindows[-1]=array.array("f",ContigG["RDWindows"][...])
        TheContig.AverageClipLengths[-1]=array.array("f",ContigG["AverageClipLengths"][...])
        if TheContig.Length!=ContigG.attrs["Length"]:
            print(gettime()+"WARN: contig %s length doesn't matchup."%c,file=sys.stderr)
        DRPTs=ContigG["DRPs"][...]
        DRPs=[0]*len(DRPTs)
        for i in range(len(DRPTs)):
            D=DRPTs[i]
            DRPs[i]=DRP(D[0],D[1],D[2],D[3],D[4])
        TheContig.DRPs[-1]=DRPs
        TheContig.ContigSampleReadCounts[-1]=ContigG.attrs["ContigSampleReadCounts"]
    HDF5File.close()
    return (SampleName,SampleReadCount,mygenome)

def readHDF5DataAll(TheGenome, SampleNames, FileNames):
    mygenome=TheGenome
    if True or g.ThreadNR==1:
        for i in range(len(FileNames)):
            Result=readSampleHDF5Data(mygenome,FileNames[i])
            g.SampleReadCount.append(Result[1])
            SampleNames.append(Result[0])
            g.SampleNameIndexes[Result[0]]=i
        print(gettime()+"SDs read. "+getMemUsage(),file=sys.stderr)
    else:
        ctx=mp.get_context("fork")
        pool=ctx.Pool(g.ThreadNR)
        args=[]
        for i in range(len(FileNames)):
            args.append((mygenome.genVacant(),FileNames[i]))
        addPool(pool)
        Results=pool.starmap(readSampleHDF5Data,args)
        delPool()
        pool.terminate()
        for j in range(len(Results)):
            r=Results[j]
            SampleName=r[0]
            g.SampleReadCount.append(r[1])
            SampleNames.append(SampleName)
            g.SampleNameIndexes[SampleName]=j
            r=r[2]
            for i in range(len(mygenome.Contigs)):
                c=mygenome.Contigs[i]
                SampleRDWindows=r.Contigs[i].RDWindows[0]
                SampleAverageClipLengths=r.Contigs[i].AverageClipLengths[0]
                ContigSampleReadCount=r.Contigs[i].ContigSampleReadCounts[0]
                SampleDRPs=r.Contigs[i].DRPs[0]
                Ploidy=r.Contigs[i].Ploidies[0]
                c.addSampleWithData(SampleName,SampleRDWindows,SampleAverageClipLengths,ContigSampleReadCount,SampleDRPs,Ploidy)
            mygenome.SampleNames.append(SampleName)
            mygenome.SampleN+=1

def readSampleSDDataToHDF5(FileName):
    print(gettime()+"Reading reads from %s..."%(FileName),file=sys.stderr)
    SDFile=open(FileName,"rb")
    try:
        SampleData=pickle.load(SDFile)
    except Exception as e:
        print("[ERROR]:%s, FileName:%s"%(e,FileName),file=sys.stderr)
        raise(e)
    SDFile.close()

    HD5FileName=FileName+".hdf5" if FileName[-3:]!=".sd" else FileName[:-3]+".hdf5"
    HD5File=h5py.File(HD5FileName,'w')

    HD5File.attrs["SampleName"]=SampleData["SampleName"]
    HD5File.attrs["SampleReadCount"]=SampleData["SampleReadCount"]
    HD5File.attrs["RDWindowSize"]=SampleData["RDWindowSize"]
    for c in SampleData["Contigs"]:
        ContigG=HD5File.create_group(c["Name"])
        ContigG.attrs["Length"]=c["Length"]
        ContigG.attrs["ContigSampleReadCounts"]=c["ContigSampleReadCounts"]
        ContigG["RDWindows"]=c["RDWindows"]
        ContigG["AverageClipLengths"]=c["AverageClipLengths"]
        DRPs=[]
        for D in c["DRPs"]:
            DRPs.append((D.Start,D.End,D.InnerStart,D.InnerEnd,D.SupportedVariantType))
        ContigG["DRPs"]=DRPs

    HD5File.close()

def readSDDataToHDF5All(FileNames):
    for i in range(len(FileNames)):
        Result=readSampleSDDataToHDF5(FileNames[i])
    print(gettime()+"SDs read and HDF5s made. "+getMemUsage(),file=sys.stderr)

def readSamToRDF(thegenome,FilePath,ReferencePath,WindowSize):
    print(gettime()+"Reading reads from %s..."%(FilePath),file=sys.stderr)
    SamFile=pysam.AlignmentFile(FilePath,"rb",reference_filename=ReferencePath)
    Name=SamFile.header["RG"][0]["SM"]
    thegenome.addSample(Name)
    ReadCount=0
    SampleIndex=0
    UnmappedCount=0
    #Reads=SamFile.fetch(until_eof=True)
    #Reads:pysam.libcalignmentfile.IteratorRowAll
    try:
#        while 1:
        for read in SamFile:
#            try:
#                read=next(SamFile)
#            except Exception as e:
#                if type(e)==StopIteration:
#                    break
#                print("WARN:(%s)%s"%(type(e),e),file=sys.stderr)
#                continue
            ReadCount+=1
            read: pysam.AlignedSegment
            if read.is_unmapped:
                UnmappedCount+=1
                #ReadCount+=1#Unmapped reads also count(to accurately measure the coverage)
            else:
                TheContig=thegenome.get(read.reference_name)#This is slightly faster than using getcontigid, and is correct at all time
                if TheContig==None:
                    continue
                TheContig.RDWindows[SampleIndex][int((int((read.reference_start+read.reference_end)/2))/WindowSize)]+=1
                #ReadCount+=1
    except Exception as e:
        print("WARN:[Sample:%s] %s"%(Name,e),file=sys.stderr)
        return False
    SamFile.close()
    print(gettime()+"Sample %s read."%(Name),file=sys.stderr)
    print(gettime()+"Storing rd data for %s..."%(Name),file=sys.stderr)
    writeSampleRDData(thegenome,Name,SampleIndex,ReadCount,WindowSize,FilePath)
    print(gettime()+"Stored rd data for %s."%(Name),file=sys.stderr)
    return True

def writeRDData(mygenome,SampleNames,SampleReadCount):
    for i in range(len(SampleNames)):
        writeSampleRDData(mygenome,SampleNames[i],i,SampleReadCount[i],g.RDWindowSize,g.SamplePaths[i])
    return

def getRDFPath(SamPath,Path="data"):
    SamFileName=SamPath.split("/")[-1].split("\\")[-1]
    return "%s/%s.rdf"%(Path,SamFileName)

def getSDPath(SamPath,Path="data"):
    SamFileName=SamPath.split("/")[-1].split("\\")[-1]
    return "%s/%s.sd"%(Path,SamFileName)

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