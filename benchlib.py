#Written by CRE. 20201013.
def calcOverlap(B1,E1,B2,E2):
    if B1<=E2 and B1>=B2:
        Overlap=min(E2-B1,E1-B1)
    elif E1>=B2 and E1<=E2:
        Overlap=E1-B2
    elif B2>=B1 and B2<=E1:#1 covers 2
        Overlap=E2-B2
    else:
        Overlap=-1
    return Overlap

def inclusion(In1,In2):#return 0: not overlapped, 1: overlapped, 2: In1 include In2, 3: In2 inlcude In1, 4: identical
    if In1==In2:
        return 4
    elif In1[0]>=In2[0] and In1[1]<=In2[1]:
        return 3
    elif In1[0]<=In2[0] and In1[1]>=In2[1]:
        return 2
    elif In1[0]>=In2[0] and In1[0]<=In2[1] or In1[1]>=In2[0] and In1[1]<=In2[1]:
        return 1
    return 0


class Variant:
    def __init__(self,Sample="",Type="",Chrom="",Start=0,End=0):
        self.Type=Type
        self.Chrom=Chrom
        self.Start=Start
        self.End=End
        self.Sample=Sample
        self.Length=self.End-self.Start
    
    def __str__(self):
        return "[%s]%s:%s-%s,%s,%s"%(self.Sample,self.Chrom,self.Start,self.End,self.End-self.Start,self.Type)
    
    def match(self,Other,Mode=0.8):#if Mode="inclusion" use inclusion, else use Mode percentage overlap
        if self.Type!=Other.Type or self.Chrom!=Other.Chrom:
            return False
        if Mode=="inclusion":
            if inclusion((self.Start,self.End),(Other.Start,Other.End))>1:
                return True
            return False
        else:
            Overlap=calcOverlap(self.Start,self.End,Other.Start,Other.End)
            if Overlap/self.Length>=Mode and Overlap/Other.Length>=Mode:
                return True
            return False

    def __lt__(self,other):
        if self.Chrom<other.Chrom:
            return True
        if self.Chrom>other.Chrom:
            return False
        if self.Start<other.Start:
            return True
        if self.Start>other.Start:
            return False
        if self.End<other.End:
            return True
        return False
    def __gt__(self,other):
        if self<other:
            return False
        if self.Chrom==other.Chrom and self.Start==other.Start and self.End==other.End:
            return False
        return True

class Sample:
    def __init__(self,Name):
        self.Name=Name
        self.Variants=[]
        self.IndexParameter=10000
        #self.IndexedLength=0#list length last time indexing taken place
        self.Index={}#by chrom
        self.Chroms=[]
        self.ChromIndex=None
    
    def addVariant(self,NewVariant):
        self.Variants.append(NewVariant)
        AddChrom=True
        for c in self.Chroms:
            if c==NewVariant.Chrom:
                AddChrom=False
                break
        if AddChrom:
            self.Chroms.append(NewVariant.Chrom)
        self.Chroms.sort()
        self.Variants.sort()
    
    def addVariantRaw(self,NewVariant):#without sorting, need call force sort afterwards
        self.Variants.append(NewVariant)
        AddChrom=True
        for c in self.Chroms:
            if c==NewVariant.Chrom:
                AddChrom=False
                break
        if AddChrom:
            self.Chroms.append(NewVariant.Chrom)
    
    def forceSort(self):
        self.Variants.sort()
        self.Chroms.sort()

    def forceIndex(self):#must call after sorted
        self.ChromIndex=set()
        for c in self.Chroms:
            self.ChromIndex.add(c)
        if len(self.Variants)==0:
            self.Index={}
        else:
            self.Index={}
            IndexLength=int(self.Variants[-1].Start/self.IndexParameter)
            FirstStartOccurredAt=0
            if len(self.Chroms)!=0:
                Chrom=self.Chroms[0]
                self.Index[Chrom]=[0]
                i=0
                while FirstStartOccurredAt<len(self.Variants):
                    if self.Variants[FirstStartOccurredAt].Start>=i*self.IndexParameter:
                        self.Index[Chrom].append(FirstStartOccurredAt)
                        i+=1
                    if self.Variants[FirstStartOccurredAt].Chrom!=Chrom:
                        Chrom=self.Variants[FirstStartOccurredAt].Chrom
                        self.Index[Chrom]=[0]
                        i=0
                    FirstStartOccurredAt+=1
    
    def getStartIndex(self,Chrom,Start):
        if self.ChromIndex!=None:
            if Chrom not in self.ChromIndex:
                return len(self.Variants)
        else:
            NoChrom=True
            for c in self.Chroms:
                if c==Chrom:
                    NoChrom=False
                    break
            if NoChrom:
                return len(self.Variants)
        i=0
        if len(self.Index)!=0:
            if Chrom in self.Index.keys():
                try:
                    i=self.Index[Chrom][int(Start/self.IndexParameter)]
                except:#beyond index
                    i=self.Index[Chrom][-1]
        while i<len(self.Variants) and self.Variants[i].Start<Start:
            i+=1
        return i
    
    def match(self,other,Mode=0.8):# returns ([All matched self variants],per other matched variant)
        Result=[]
        if Mode!="inclusion":
            ratio=Mode
        for v in other.Variants:
            start=0 if v.Start-v.Length/ratio-1<0 else int(v.Start-v.Length/ratio-1)
            end=v.End+1
            #SI=self.getStartIndex(v.Chrom,start)
            #EI=self.getStartIndex(v.Chrom,end)
            SI=0
            EI=len(self.Variants)
            TempPair=([],v)
            for i in range(SI,EI):
                if self.Variants[i].match(v,Mode):
                    TempPair[0].append(self.Variants[i])
            if len(TempPair[0])!=0:
                Result.append(TempPair)
        return Result

class Record:
    def __init__(self,Chrom="",Start=0,End=0):
        self.Chrom=Chrom
        self.Start=Start
        self.End=End
        self.Length=self.End-self.Start
        self.Variants=[]
        self.Samples={}
    def addVariant(self,NewVariant):
        self.Variants.append(NewVariant)
        if NewVariant.Sample in self.Samples.keys():
            self.Samples[NewVariant.Sample].append(NewVariant)
        else:
            self.Samples[NewVariant.Sample]=[NewVariant]

import pysam

class VariantRecords:
    def __init__(self,Name=""):
        self.Variants=[]
        self.Name=Name
        self.Records=[]
        self.Samples={}
    
    def getSample(self, SampleName):
        if SampleName not in self.Samples.keys():
            self.Samples[SampleName]=Sample(SampleName)
        return self.Samples[SampleName]
    
    def matchSample(self,Other,SampleName,Mode=0.8):#other as standard
        try:
            return self.Samples[SampleName].match(Other.Samples[SampleName],Mode)
        except Exception:
            return None
    
    def matchAll(self,Other,Mode=0.8):#other as standard
        Result={}
        Result["Samples"]={}
        RSamples=Result["Samples"]
        Result["Self"]=self
        Result["Other"]=Other
        for SN in Other.Samples.keys():
            if SN in self.Samples.keys():
                RSamples[SN]=self.Samples[SN].match(Other.Samples[SN],Mode)
            else:
                RSamples[SN]=None
        return Result
    
    def interpret(self,maresult,PR=False):
        OtherN=0
        SelfN=0
        Matched=0
        SelfMatched=0
        OtherDupN=0
        OtherDelN=0
        MatchedDup=0
        MatchedDel=0
        SelfDupN=0
        SelfDelN=0
        SelfMatchedDup=0
        SelfMatchedDel=0

        RSamples=maresult["Samples"]
        Other=maresult["Other"]

        SelfN=len(self.Variants)
        OtherN=len(Other.Variants)
        for v in self.Variants:
            if v.Type.upper()=="DEL":
                SelfDelN+=1
            elif v.Type.upper()=="DUP":
                SelfDupN+=1
        for v in Other.Variants:
            if v.Type.upper()=="DEL":
                OtherDelN+=1
            elif v.Type.upper()=="DUP":
                OtherDupN+=1
        
        for SN in Other.Samples.keys():
            if RSamples[SN]!=None:
                Matched+=len(RSamples[SN])
                for i in range(len(RSamples[SN])):
                    SelfMatched+=len(RSamples[SN][i][0])
                    if RSamples[SN][i][1].Type.upper()=="DEL":
                        MatchedDel+=1
                        SelfMatchedDel+=len(RSamples[SN][i][0])
                    elif RSamples[SN][i][1].Type.upper()=="DUP":
                        MatchedDup+=1
                        SelfMatchedDup+=len(RSamples[SN][i][0])
        Out= "All: Sensitivity: %s (%s/%s), PPV: %s (%s/%s)\nDel: Sensitivity: %s (%s/%s), PPV: %s (%s/%s)\nDup: Sensitivity: %s (%s/%s), PPV: %s (%s/%s)"%(\
            Matched/OtherN if OtherN!=0 else 0,Matched,OtherN,SelfMatched/SelfN if SelfN!=0 else 0,SelfMatched,SelfN\
                ,MatchedDel/OtherDelN if OtherDelN!=0 else 0,MatchedDel,OtherDelN,SelfMatchedDel/SelfDelN if SelfDelN!=0 else 0,SelfMatchedDel,SelfDelN\
                    ,MatchedDup/OtherDupN if OtherDupN!=0 else 0,MatchedDup,OtherDupN,SelfMatchedDup/SelfDupN if SelfDupN!=0 else 0,SelfMatchedDup,SelfDupN
                )#PPV:positive predictive value,1-FDR
        if PR:
            for SN in Other.Samples.keys():
                if RSamples[SN]!=None:
                    for i in range(len(RSamples[SN])):
                        Out+="\nGS:%s..."%RSamples[SN][i][1]
                        for j in range(len(RSamples[SN][i][0])):
                            if j!=0:
                                Out+=","
                            Out+="%s"%RSamples[SN][i][0][j]
        return Out
    
    def parseVcfCNV(self,filename,contigs=None,samples=None,MinLength=0):
        vcf=pysam.VariantFile(filename,"r")
        result=[]
        included={}
        includedSamples={}
        if contigs!=None:
            for c in contigs:
                included[c]=True
        if samples!=None:
            for s in samples:
                includedSamples[s.upper()]=True
        for record in vcf.fetch():
            End=record.pos
            try:
                End=record.pos+abs(record.info["SVLEN"])
            except:
                try:
                    End=record.pos+abs(record.info["SVLEN"][0])
                except:
                    try:
                        End=record.info["END"]
                    except:
                        try:
                            End=record.stop
                        except:
                            pass
            try:
                Length=End-record.pos
                if Length<MinLength:
                    continue
                chrom=record.chrom
                #if len(record.chrom)<3 and (int(record.chrom[:2])<23 or record.chrom.upper()=="X" or record.chrom.upper()=="Y"):
                #        chrom="chr"+record.chrom
                if contigs!=None:
                    if not chrom in included.keys():
                        continue
                CNAlts=[1]
                tempRecord=Record(chrom,record.pos,End)
                for a in record.alts:
                    if "<CN" in a:
                        CNAlts.append(int(a[a.find("<CN")+3:-1]))
                    elif "DEL" in a.upper():
                        CNAlts.append(0)
                    elif "DUP" in a.upper():
                        CNAlts.append(2)
                    else:
                        CNAlts.append(1)
                try:
                    for s in record.samples:
                        SampleName=record.samples[s].name
                        if samples==None or SampleName in includedSamples.keys():
                            if record.samples[s].allele_indices!=(0,0):
                                CN=2
                                if record.info["SVTYPE"]=="DEL" or record.info["SVTYPE"]=="DUP" or record.info["SVTYPE"]=="CNV":
                                    CN=CNAlts[record.samples[s].allele_indices[0]]+CNAlts[record.samples[s].allele_indices[1]]
                                if CN==2:
                                    continue
                                else:
                                    if CN<2:
                                        SVType="DEL"
                                    else:
                                        SVType="DUP"
                                    self.Variants.append(Variant(SampleName,SVType,chrom,record.pos,End))
                                    self.getSample(SampleName).addVariantRaw(self.Variants[-1])
                                    tempRecord.addVariant(self.Variants[-1])
                except Exception as e:
                    continue
                self.Records.append(tempRecord)
            except:
                continue
        for SN in self.Samples.keys():
            self.Samples[SN].forceSort()
            self.Samples[SN].forceIndex()
        vcf.close()
        return self.Records