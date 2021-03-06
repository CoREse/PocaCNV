from contig import *
from readdepth import *
import globals as g

def readRDDataAndSaveRDF(thegenome, SamplePaths):
    RDFPaths=SamplePaths.copy()
    SAMPaths=[]
    RtS=[]
    for i in range(len(SamplePaths)):
        #print(gettime()+"Reading %s..."%(g.SamplePaths[i]),file=sys.stderr)
        RtS.append(-1)
        if SamplePaths[i].split(".")[-1]=="rdf":
            continue
        else:
            RDFPaths[i]=getRDFPath(SamplePaths[i])
            SAMPaths.append(SamplePaths[i])
            RtS[-1]=len(SAMPaths)-1
    if len(SAMPaths)>0:
        print(gettime()+"Start reading SAM files and save them to RDF files for following analysis...",file=sys.stderr)
        if g.ThreadNR>1:
            ctx=mp.get_context()
            pool=ctx.Pool(g.ThreadNR)
            args=[]
            for i in range(len(SAMPaths)):
                args.append((thegenome.genVacant(),SAMPaths[i],g.ReferencePath,g.RDWindowSize))
            addPool(pool)
            Results=pool.starmap(readSamToRDF,args)
            delPool()
            pool.close()
        else:
            Results=[]
            for i in range(len(SAMPaths)):
                Results.append(readSamToRDF(thegenome.genVacant(),SAMPaths[i],g.ReferencePath,g.RDWindowSize))
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
            print(gettime()+"All SAMFile(s) read and RDFs made. "+getMemUsage(),file=sys.stderr)
        else:
            print(gettime()+"A part of SAMFile(s) read and RDFs made. "+getMemUsage(),file=sys.stderr)
    return RDFPaths

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