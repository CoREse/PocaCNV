import pysam

import sys

def calcOverlap(B1,E1,B2,E2):
    if B1<=E2 and B1>=B2:
        Overlap=min(E2-B1,E1-B1)
    elif E1>=B2 and E1<=B2:
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

class interval:
    
    def __init__(self, chrom, start, end, type):  
        self.chrom = chrom 
        self.start = start  
        self.end = end
        self.type = type
    
    def __str__(self):
        return "%s:%s-%s,%s,%s"%(self.chrom,self.start,self.end,self.end-self.start,self.type)

    def get_chrom(self):
        return self.chrom
    
    def get_start(self):
        return self.start
    
    def get_end(self):
        return self.end
    
    def get_type(self):
        return self.type

    def get_length(self):
        return self.end - self.start + 1
    
    def is_overlap_pure(self, other, percent=0):
        if self.start > other.get_end() or self.end < other.get_start():
            return False
        if self.start >= other.get_start(): 
            if self.end <= other.get_end():
                overlap = self.get_length()
            else:
                overlap = other.get_end() - self.start + 1
            if float(overlap) / float(self.get_length()) >= percent and float(overlap) / float(other.get_length()) >= percent:
                return True
            else:
                return False
        if self.start < other.get_start():
            if other.get_end() <= self.end:
                overlap = other.get_length()
            else:
                overlap = self.end - other.get_start() + 1
            if float(overlap) / float(self.get_length()) >= percent and float(overlap) / float(other.get_length()) >= percent:               
                return True
            else:
                return False
    
    def is_overlap(self, other):
        if self.chrom != other.get_chrom() or self.type != other.get_type():
            return False
        elif self.start > other.get_end() or self.end < other.get_start():
            return False
        else:
            return True
    
    def is_overlap_by_percent(self, other, percent):
        if self.chrom != other.get_chrom() or self.type != other.get_type():
            return False
        if self.start > other.get_end() or self.end < other.get_start():
            return False
        if self.start >= other.get_start(): 
            if self.end <= other.get_end():
                overlap = self.get_length()
            else:
                overlap = other.get_end() - self.start + 1
            if float(overlap) / float(self.get_length()) >= percent and float(overlap) / float(other.get_length()) >= percent:
                return True
            else:
                return False
        if self.start < other.get_start():
            if other.get_end() <= self.end:
                overlap = other.get_length()
            else:
                overlap = self.end - other.get_start() + 1
            if float(overlap) / float(self.get_length()) >= percent and float(overlap) / float(other.get_length()) >= percent:               
                return True
            else:
                return False
    def is_overlap_by_percent_or_near(self, other, percent, near):#for insertion, if break point is near within near considered overlapped
        if self.chrom != other.get_chrom() or self.type != other.get_type():
            return False
        if self.type=="INS":
            spoint=int((self.start+self.end)/2)
            opoint=int((other.start+other.end)/2)
            if abs(spoint-opoint)<=near:
                return True
        if self.start > other.get_end() or self.end < other.get_start():
            return False
        if self.start >= other.get_start(): 
            if self.end <= other.get_end():
                overlap = self.get_length()
            else:
                overlap = other.get_end() - self.start + 1
            if float(overlap) / float(self.get_length()) >= percent and float(overlap) / float(other.get_length()) >= percent:
                return True
            else:
                return False
        if self.start < other.get_start():
            if other.get_end() <= self.end:
                overlap = other.get_length()
            else:
                overlap = self.end - other.get_start() + 1
            if float(overlap) / float(self.get_length()) >= percent and float(overlap) / float(other.get_length()) >= percent:               
                return True
            else:
                return False
    def is_included_by(self,other):
        if self.start>=other.get_start() and self.end<=other.get_end():
            return True
        return False

def parse_cnvnator(input_file, min_length=0, max_length=sys.maxsize):
    records = []
    f = open(input_file, 'r')
    for line in f:
        record = line.split()
        type = record[0]
        chrom = record[1].split(':')[0]
        region = record[1].split(':')[1].split('-')
        start = int(region[0])
        end = int(region[1])
        if type == 'deletion':
            tmp = interval(chrom, start, end, 'DEL')
            if tmp.get_length() >= min_length and tmp.get_length() < max_length:
                    records.append(tmp)
        if type == 'duplication':
            tmp = interval(chrom, start, end, 'DUP')
            if tmp.get_length() >= min_length and tmp.get_length() < max_length:
                    records.append(tmp)           
    f.close()
    return records

class Match:
    def __init__(self,Type=0,gs=None,Target=None):#Type:0:percent or near; 1:inclusion
        self.Type=Type
        self.gs=gs
        self.Target=Target
        self.gsp=0
        self.tgp=0
        self.analyze()

    def addGS(self,gs):
        self.gs=gs
        self.analyze()
    def addTarget(self,Target):
        self.Target=Target
        self.analyze()
    def analyze(self):
        if self.gs==None or self.Target==None:
            return
        Over=calcOverlap(self.gs.start,self.gs.end,self.Target.start,self.Target.end)
        self.gsp=Over/(self.gs.end-self.gs.start)
        self.tgp=Over/(self.Target.end-self.Target.start)
    def __str__(self):
        return "gs:%s,tg:%s,gsp:%d%%,tgp:%d%%"%(self.gs,self.Target,self.gsp*100,self.tgp*100)+(", gs included by tg"if self.gsp>=1 else "")+(", tg included by gs"if self.tgp>=1 else"")

            
def calculate_sensitivity_by_percent(goldstandard, callset, percent,near=-1,printmatch=False):
    count = 0
    gaincount=0
    losscount=0
    inscount=0
    gaintotal=0
    losstotal=0
    instotal=0
    results=[]
    for interval1 in goldstandard:
        if interval1.type=="DUP":
            gaintotal+=1
        elif interval1.type=="DEL":
            losstotal+=1
        elif interval1.type=="INS":
            instotal+=1
        tag = False    
        for interval2 in callset:
            if interval1.type!=interval2.type:
                continue
            if interval1.is_overlap_by_percent_or_near(interval2, percent,near):
                tag = True
                results.append(Match(0,interval1,interval2))
        if tag == True:
            if printmatch:
                print(interval1)
            count = count + 1
            if interval1.type=="DUP":
                gaincount+=1
            elif interval1.type=="DEL":
                losscount+=1
            elif interval1.type=="INS":
                inscount+=1
    print('All=%s, Found=%s, Recalled=%s(percentage=%s)' % (len(goldstandard),str(len(callset)),count,percent))
    print('Sensitivity: All:%s, Dup:%s(%s/%s), Del:%s(%s/%s), Ins:%s(%s/%s)' %(str(float(count) / float(len(goldstandard))) if len(goldstandard)!=0 else "N/A",\
        gaincount/gaintotal if gaintotal !=0 else "N/A",gaincount,gaintotal,\
            losscount/losstotal if losstotal!=0 else "N/A",losscount,losstotal,\
                inscount/instotal if instotal!=0 else "N/A",inscount,instotal))      
    return results

def calculate_fdr_by_percent(goldstandard, callset, percent, near=-1, printmatch=False):
    count = 0
    gaincount=0
    losscount=0
    inscount=0
    gaintotal=0
    losstotal=0
    instotal=0
    results=[]
    for interval1 in callset:
        if interval1.type=="DUP":
            gaintotal+=1
        elif interval1.type=="DEL":
            losstotal+=1
        elif interval1.type=="INS":
            instotal+=1
        tag = False    
        for interval2 in goldstandard:
            if interval1.type!=interval2.type:
                continue
            if interval1.is_overlap_by_percent_or_near(interval2, percent,near):
                tag = True
                results.append(Match(0,interval2,interval1))
        if tag == True:
            if printmatch:
                print(interval1)
            count = count + 1
            if interval1.type=="DUP":
                gaincount+=1
            elif interval1.type=="DEL":
                losscount+=1
            elif interval1.type=="INS":
                inscount+=1
    print('All=%s, Found=%s, Correct=%s(percentage=%s)' % (len(goldstandard),str(len(callset)),count,percent))
    print('FDR(Right/All): All:%s, Dup:%s(%s/%s), Del:%s(%s/%s), Ins:%s(%s/%s)' \
        %(str(1-float(count) / float(len(callset))) if len(callset)!=0 else "N/A",1-gaincount/gaintotal if gaintotal !=0 else "N/A",gaincount,gaintotal,\
            1-losscount/losstotal if losstotal!=0 else "N/A",losscount,losstotal,\
                1-inscount/instotal if instotal!=0 else "N/A",inscount,instotal))
    return results

##consider caught if callset fully include goldstandard
def calculate_sensitivity_by_inclusion(goldstandard, callset, printmatch=False):
    count = 0
    gaincount=0
    losscount=0
    inscount=0
    gaintotal=0
    losstotal=0
    instotal=0
    results=[]
    for interval1 in goldstandard:
        if interval1.type=="DUP":
            gaintotal+=1
        elif interval1.type=="DEL":
            losstotal+=1
        elif interval1.type=="INS":
            instotal+=1
        tag = False    
        for interval2 in callset:
            if interval1.type!=interval2.type:
                continue
            if interval1.is_included_by(interval2) or interval2.is_included_by(interval1):
                results.append(Match(1,interval1,interval2))
                tag = True
        if tag == True:
            if printmatch:
                print(interval1)
            count = count + 1
            if interval1.type=="DUP":
                gaincount+=1
            elif interval1.type=="DEL":
                losscount+=1
            elif interval1.type=="INS":
                inscount+=1
    print('All=%s, Found=%s, Recalled=%s(by inclusion)' % (len(goldstandard),str(len(callset)),count))
    print('Sensitivity: All:%s, Dup:%s(%s/%s), Del:%s(%s/%s), Ins:%s(%s/%s)' %(str(float(count) / float(len(goldstandard)))if len(goldstandard)!=0 else "N/A",\
        gaincount/gaintotal if gaintotal !=0 else "N/A",gaincount,gaintotal,\
            losscount/losstotal if losstotal!=0 else "N/A",losscount,losstotal,\
                inscount/instotal if instotal!=0 else "N/A",inscount,instotal))   
    return results

def calculate_fdr_by_inclusion(goldstandard, callset, printmatch=False):
    count = 0
    gaincount=0
    losscount=0
    inscount=0
    gaintotal=0
    losstotal=0
    instotal=0
    results=[]
    for interval1 in callset:
        if interval1.type=="DUP":
            gaintotal+=1
        elif interval1.type=="DEL":
            losstotal+=1
        elif interval1.type=="INS":
            instotal+=1
        tag = False    
        for interval2 in goldstandard:
            if interval1.type!=interval2.type:
                continue
            if interval2.is_included_by(interval1) or interval1.is_included_by(interval2):
                tag = True
                results.append(Match(1,interval2,interval1))
        if tag == True:
            if printmatch:
                print(interval1)
            count = count + 1
            if interval1.type=="DUP":
                gaincount+=1
            elif interval1.type=="DEL":
                losscount+=1
            elif interval1.type=="INS":
                inscount+=1
    print('All=%s, Found=%s, Correct=%s(by inclusion)' % (len(goldstandard),str(len(callset)),count))
    print('FDR(Right/All): All:%s, Dup:%s(%s/%s), Del:%s(%s/%s), Ins:%s(%s/%s)' \
        %(str(1-float(count) / float(len(callset))) if len(callset)!=0 else "N/A",1-gaincount/gaintotal if gaintotal !=0 else "N/A",gaincount,gaintotal,\
            1-losscount/losstotal if losstotal!=0 else "N/A",losscount,losstotal,\
                1-inscount/instotal if instotal!=0 else "N/A",inscount,instotal))
    return results

def parse_vcf(filename,contigs=None,samples=None):
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
            chrom=record.chrom
            if samples!=None and len(samples)==1:
                sa=None
            #if len(record.chrom)<3 and (int(record.chrom[:2])<23 or record.chrom.upper()=="X" or record.chrom.upper()=="Y"):
            #        chrom="chr"+record.chrom
            if contigs!=None:
                if not chrom in included.keys():
                    continue
            if samples!=None:
                try:
                    skip=True
                    if len(record.samples)!=0:
                        for s in record.samples:
                            if record.samples[s].name in includedSamples.keys() or (len(samples)==1 and samples[0] in record.samples[s].name):
                                if record.samples[s].allele_indices!=(0,0):
                                    sa=record.samples[s].allele_indices
                                    skip=False
                        if skip:
                            continue
                    else:
                        if not (record.info["SAMPLE"].upper() in includedSamples.keys() or (len(samples)==1 and samples[0] in record.info["SAMPLE"])):
                            continue
                except Exception as e:
                    continue
            temp=None
            if record.info["SVTYPE"]=="INS":
                temp=interval(chrom,record.pos,record.stop+1,"INS")
            elif record.info["SVTYPE"]=="DEL":
                temp=interval(chrom,record.pos,End,"DEL")
            elif record.info["SVTYPE"]=="DUP":
                temp=interval(chrom,record.pos,End,"DUP")
            elif record.info["SVTYPE"]=="CNV":
                CNs=[]
                for a in record.alts:
                    if "<CN" in a:
                        CNs.append(int(a[a.find("<CN")+3:-1]))
                    else:
                        CNs.append(1)
                if len(CNs)==0:
                    continue
                if len(CNs)>1:
                    Type="CNV"
                    if samples!=None and len(samples)==1:
                        Al=1
                        if sa[0]==0:
                            Al=CNs[sa[1]]
                        if sa[1]==0:
                            Al=CNs[sa[0]]
                        if sa[0]==sa[1]:
                            Al=CNs[sa[0]]
                        if Al==0:
                            Type="DEL"
                        elif Al>1:
                            Type="DUP"
                elif CNs[0]==0:
                    Type="DEL"
                else:
                    Type="DUP"
                temp=interval(chrom,record.pos,End,Type)
            else:
                continue
            result.append(temp)
        except:
            pass
    vcf.close()
    return result

def parse_my(filename,contigs=None, samples=None, format=None):
    if format.upper()=="VCF" or filename[-4:]==".vcf" or filename[-4:]==".bcf":
        return parse_vcf(filename,contigs,samples)
    myout=open(filename,"r")
    records=[]
    for line in myout:
        try:
            line=line.replace(" ","")
            sl=line.split(",")
            blockstart=int(sl[0].split("-")[0].split(":")[1])
            blockend=int(sl[0].split("-")[1].split(":")[1])
            breakstart=int(sl[2].split(":")[-1])
            breakend=int(sl[3].split(":")[1].split("]")[0])
            temp=interval(sl[0].split(":")[0],breakstart,breakend,"DEL" if sl[1]=="DEL" else ("INS" if sl[1]=="INS" else "DUP"))
            vsamples=sl[-2].split()
            skip=True
            for s in samples:
                for v in vsamples:
                    if s in v:
                        skip=False
            if skip:
                continue
            records.append(temp)
        except:
            continue
    myout.close()
    return records