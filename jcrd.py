import pysam
import sys
import time
import globals
from calling import *
from readpair import *
from readdepth import *
from filters import *
from contig import *
from array import array
from report import *
from joint import uniformlyCombine
import json
import multiprocessing as mp
from jcio import *

import os
import psutil

if __name__ == "__main__":
    process = psutil.Process(os.getpid())

    mp.set_start_method("spawn")
    g=globals
    #parameters
    #g.WriteRDData=False
    g.WriteRDDataOnly=False
    g.WriteMixedRDData=False
    g.Contigs=set()#if not vacant, contain only those contigs
    g.ExcludeRegionsFileName=None
    g.DumpCandidates=None
    g.LoadCandidates=None
    #g.ExcludeRegionsFileName="data/mdust-v28-p1.bed"
    #g.ExcludeRegionsFileName="data/exclusion_regions.txt"

    g.WriteRDDataTo=None

    g.ReferencePath=None
    g.SamplePaths=[]
    g.ThreadNR=None
    def usage():
        print("""python3 jcrd [ARGs] -T ReferenceFile SampleFile1.sam/bam/cram SampleFile2.sam/bam/cram ...
    Arguments:
        -T,--Reference FILENAME   give referencefile(fasta)(str)
        -WO                       Write RD data only(no calling)
        -C             ContigName Contain (only) contig, can be used multiple times(str)
        -J,-j          ThreadNum  Default:8. Thread number(int)
        -JR            ThreadNum  Thread number for reading sam files, default same as threadnum
        -WS            Size       Default:100. WindowSize(int)
        -WM                       Write Mixed RD data and other data
        -LC            FILENAME   Load candidates(include all data) from dumped file(str)
        -DC            FILENAME   Dump candidates to json FILENAME(str)
        -EN            EDataName  EDataName if you want to write EData(for training)
        -EP            EDataPath  Default:data. Path that stores EData
        """,file=sys.stderr)
        exit(0)
    def getParas():
        i=1
        while i< len(sys.argv):
            try:
                a=sys.argv[i]
                if a=='-T' or a=="--Reference":
                    g.ReferencePath=sys.argv[i+1]
                    i+=1
                #elif a=='-WT':
                #    g.WriteRDDataTo=sys.argv[i+1]
                #    i+=1
                elif a=='-W':
                    print("WARNING: -W is obsolete(all rdf is writen by default), try -WT",file=sys.stderr)
                elif a=='-WO':
                    g.WriteRDData=True
                    g.WriteRDDataOnly=True
                elif a=="-WM":
                    g.WriteMixedRDData=True
                elif a=='-C':
                    g.Contigs.add(sys.argv[i+1])
                    i+=1
                elif a=='-J' or a=='-j':
                    g.ThreadN=int(sys.argv[i+1])
                    i+=1
                elif a=='-JR':
                    g.ThreadNR=int(sys.argv[i+1])
                    i+=1
                elif a=='-WS':
                    g.RDWindowSize=int(sys.argv[i+1])
                    i+=1
                elif a=="-DC":
                    g.DumpCandidates=sys.argv[i+1]
                    i+=1
                elif a=="-LC":
                    g.LoadCandidates=sys.argv[i+1]
                    i+=1
                elif a=="-EN":
                    g.EDataName=sys.argv[i+1]
                    i+=1
                elif a=="-EP":
                    g.EDataPath=sys.argv[i+1]
                    i+=1
                else:
                    g.SamplePaths.append(sys.argv[i])
            except Exception as e:
                print(e,file=sys.stderr)
                usage()
            i+=1
        #if g.WriteRDDataOnly:
        #    g.ThreadN=1
        if g.ThreadNR==None:
            g.ThreadNR=g.ThreadN
        if g.ReferencePath==None or len(g.SamplePaths)==0:
            usage()

    getParas()

    g.Processes.append(process)

    print(gettime()+"Run command:",end="",file=sys.stderr)
    for a in sys.argv:
        print(" %s"%a,end="",file=sys.stderr)
    print(gettime()+"Joint calling started...", file=sys.stderr)
    print(gettime()+"Reading reference...",file=sys.stderr)
    ReferenceFile=pysam.FastaFile(g.ReferencePath)
    PosCount=0
    mygenome=Genome(ReferenceFile.filename)
    for tid in range(ReferenceFile.nreferences):
        g.RefInd[ReferenceFile.references[tid]]=tid
        g.RefStartPos.append(PosCount)
        PosCount+=ReferenceFile.lengths[tid]
        if len(g.Contigs)!=0:
            if ReferenceFile.references[tid] not in g.Contigs:
                continue
        c=Contig(ReferenceFile.references[tid],ReferenceFile.lengths[tid],globals.RDWindowSize,mygenome)
        mygenome.RefID.append(tid)
        mygenome.append(c)
    ReferenceFile.close()
    RefLength=PosCount
    mygenome.setGenomeLength(RefLength)
    g.MyGenome=mygenome
    print(gettime()+"Reference %s read. Length:%s, Contigs:%s."%(g.ReferencePath,PosCount,len(mygenome)),file=sys.stderr)
    print(gettime()+"Reading samples...",file=sys.stderr)
    SDPaths=readSamAndSaveSD(mygenome,g.SamplePaths)#return SD SamplePaths
    if g.WriteRDDataOnly:
        exit(0)

    HeadReported=False
    for ci in range(len(mygenome.Contigs)):
        print(gettime()+"Starting call for contig %s..."%(mygenome.Contigs[ci].Name),file=sys.stderr)
        ContigGenome=Genome(mygenome.Name)
        TheContig=mygenome.Contigs[ci].getNew(ContigGenome)
        ContigGenome.RefID.append(mygenome.RefID[ci])
        ContigGenome.append(TheContig)


        SampleIndex=0

        if g.LoadCandidates==None:
            PairCount=0
            ReadCount=0
            LCount=0
            RCount=0
            UnmappedCount=0
            SampleNames=[]

            #RDFPaths=readRDDataAndSaveRDF(mygenome,g.SamplePaths)#return RDF SamplePaths
            if g.WriteRDDataOnly:
                exit(0)
            readSDDataAll(ContigGenome, SampleNames,SDPaths)
            ContigGenome.calcContigReadCounts()
            #print(ReadCount,PairCount,LCount,RCount,UnmappedCount,file=sys.stderr)
            #exit(0)
            globals.SampleNames=SampleNames
            if not HeadReported:
                reportVCFHeader(sys.stdout,mygenome)
                HeadReported=True
            #calclulateSequenceDepthRatio()

            if g.ExcludeRegionsFileName!=None:
                EAFile=open(g.ExcludeRegionsFileName,"r")
                ExcludedAreasByContig=readExcludedAreas(EAFile,ReferenceFile)
                EAFile.close()

            print("%s"%(getMemUsage()),file=sys.stderr)
            print(gettime()+"Samples read, calculating RD data...", file=sys.stderr)

            #for c in mygenome:
            #    analyzeRD(c.RDWindows,c.Length,c, True)
            #    print(gettime()+"%s analyzed. Memory usage: %.6sgb"%(c.Name, process.memory_info().vms/1024/1024/1024),file=sys.stderr)

            #print("Memory usage: %.6sgb"%(process.memory_info().vms/1024/1024/1024),file=sys.stderr)
            #writeMixedRDData(mygenome,ReferenceFile,SampleNames)
            #exit(0)

            RDICandidates=[]
            for c in ContigGenome:
                ploidyProcessing(c)
                c.RDICandidates=analyzeRD(c.RDWindows,c.Length,c)
                #c.RDICandidates=filtExcludedAreas(c.RDICandidates)
                RDICandidates.append(c.RDICandidates)
                print(gettime()+"%s analyzed. %s"%(c.Name, getMemUsage()),file=sys.stderr)

            '''if g.WriteMixedRDData:
                print(gettime()+"Writing mixed rdrs data...",file=sys.stderr)
                writeMixedRDData(mygenome,g.SampleNames)
                OtherData={}
                OtherData["RDWindowStandards"]={}
                OtherData["Segmentation"]={}
                OtherData["SDSegments"]={}
                for C in mygenome:
                    OtherData["RDWindowStandards"][C.Name]=C.RDWindowStandards
                    OtherData["Segmentation"][C.Name]={}
                    try:
                        OtherData["SDSegments"][C.Name]=C.SDSegments
                    except:
                        pass
                    for s in range(len(C.Intervals)):
                        SName=g.SampleNames[s]
                        OtherData["Segmentation"][C.Name][SName]=[]
                        for I in C.Intervals[s]:
                            OtherData["Segmentation"][C.Name][SName].append((I.WEnd,I.AverageRD))
                OtherData["RDWindowSize"]=g.RDWindowSize
                ODF=open("data/OtherData.json","w")
                json.dump(OtherData,ODF)
                ODF.close()'''

            print("Number of RDI candidates:%d"%(len(RDICandidates)),file=sys.stderr)
            print("%s"%(getMemUsage()),file=sys.stderr)

            Candidates=RDICandidates
            """
            for c in mygenome:
                Candidates.append(combineCandidateSets(DRCandidates,RDICandidates))
                Candidates[-1]=filtExcludedAreas(Candidates[-1])
            """
            #analyzeRD
            #analyzeDR
            #combineEvidence
            CCount=0
            for cs in Candidates:
                CCount+=len(cs)
            print(gettime()+"Number of candidates:%d. Filtering candidates in LPR..."%CCount,file=sys.stderr)
            if g.ExcludeRegionsFileName!=None:
                for i in range(len(Candidates)):
                    Candidates[i]=filtExcludedAreas(Candidates[i],ExcludedAreasByContig,ContigGenome[i])
            CCount=0
            for cs in Candidates:
                CCount+=len(cs)
            print(gettime()+"Number of filtered candidates:%d. Uniformly combining..."%(CCount),file=sys.stderr)
            #for i in range(len(Candidates)):
            #    Candidates[i]=uniformlyCombine(Candidates[i],mygenome[i])
            CCount=0
            for cs in Candidates:
                CCount+=len(cs)
            print(gettime()+"Number of uniformly combined candidates:%d. CNV calling..."%(CCount),file=sys.stderr)
        '''else:
            import pickle
            AllDataFile=open(g.LoadCandidates,"rb")
            AllData=pickle.load(AllDataFile)
            AllDataFile.close()
            Candidates=AllData["Candidates"]
            g.SampleNames=AllData["SampleNames"]
            mygenome=AllData["mygenome"]
        if g.DumpCandidates!=None:
            import pickle
            AllData={}
            AllData["SampleNames"]=g.SampleNames
            AllData["Candidates"]=Candidates
            AllData["mygenome"]=mygenome
            AllDataFile=open(g.DumpCandidates,"wb")
            AllData=pickle.dump(AllData,AllDataFile)
            AllDataFile.close()'''
        ReferenceFile=pysam.FastaFile(g.ReferencePath)
        for i in range(len(ContigGenome)):
            SVs=[]
            print(gettime()+"Calling CNV for %s"%ContigGenome[i].Name,file=sys.stderr)
            if False and g.ThreadN==1:
                for C in Candidates[i]:
                    SV=callSV(ReferenceFile,C,ContigGenome[i])
                    if SV!="":
                        SV.Chrom=ContigGenome[i].Name
                        SVs.append(SV)
            else:
                #pool=mp.Pool(g.ThreadN)
                #addPool(pool)
                #Args=[]
                #for C in Candidates[i]:
                #    Args.append((C,mygenome[i]))
                #Scores=pool.starmap(getRDScore,Args)
                Scores=getRDScores(Candidates[i],ContigGenome[i],g.ThreadN)
                print(gettime()+"Scores get. "+getMemUsage(),file=sys.stderr)
                #delPool()
                #pool.terminate()
                for j in range(len(Candidates[i])):
                    C=Candidates[i][j]
                    SV=callSV(ReferenceFile,C,ContigGenome[i],Scores[j])
                    if SV!="":
                        SV.Chrom=ContigGenome[i].Name
                        SVs.append(SV)
            SVs.sort(key=lambda s:s.BreakLeft)
            reportVCF(SVs,ReferenceFile.fetch(ContigGenome[i].Name),sys.stdout)
        ReferenceFile.close()
        print(gettime()+"Contig %s called."%(mygenome.Contigs[ci].Name),file=sys.stderr)
    print(gettime()+"All done.",file=sys.stderr)