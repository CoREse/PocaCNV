import pysam

import sys

class interval:
    
    def __init__(self, chrom, start, end, type):  
        self.chrom = chrom 
        self.start = start  
        self.end = end
        self.type = type
    
    def __str__(self):
        return "%s:%s-%s,%s"%(self.chrom,self.start,self.end,self.type)

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

            
def calculate_sensitivity_by_percent(goldstandard, callset, percent,near=-1,printmatch=False):
    count = 0
    gaincount=0
    losscount=0
    inscount=0
    gaintotal=0
    losstotal=0
    instotal=0
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

def calculate_fdr_by_percent(goldstandard, callset, percent, near=-1, printmatch=False):
    count = 0
    gaincount=0
    losscount=0
    inscount=0
    gaintotal=0
    losstotal=0
    instotal=0
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

##consider caught if callset fully include goldstandard
def calculate_sensitivity_by_inclusion(goldstandard, callset, printmatch=False):
    count = 0
    gaincount=0
    losscount=0
    inscount=0
    gaintotal=0
    losstotal=0
    instotal=0
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
            if interval1.is_included_by(interval2):
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

def calculate_fdr_by_inclusion(goldstandard, callset, printmatch=False):
    count = 0
    gaincount=0
    losscount=0
    inscount=0
    gaintotal=0
    losstotal=0
    instotal=0
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
            if interval2.is_included_by(interval1):
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
    print('All=%s, Found=%s, Correct=%s(by inclusion)' % (len(goldstandard),str(len(callset)),count))
    print('FDR(Right/All): All:%s, Dup:%s(%s/%s), Del:%s(%s/%s), Ins:%s(%s/%s)' \
        %(str(1-float(count) / float(len(callset))) if len(callset)!=0 else "N/A",1-gaincount/gaintotal if gaintotal !=0 else "N/A",gaincount,gaintotal,\
            1-losscount/losstotal if losstotal!=0 else "N/A",losscount,losstotal,\
                1-inscount/instotal if instotal!=0 else "N/A",inscount,instotal))

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
                End=record.info["END"]
            except:
                try:
                    End=record.stop
                except:
                    pass
        try:
            chrom=record.chrom
            #if len(record.chrom)<3 and (int(record.chrom[:2])<23 or record.chrom.upper()=="X" or record.chrom.upper()=="Y"):
            #        chrom="chr"+record.chrom
            if contigs!=None:
                if not chrom in included.keys():
                    continue
            if samples!=None:
                try:
                    skip=True
                    for s in record.samples:
                        if record.samples[s].name in includedSamples.keys():
                            if record.samples[s].allele_indices!=(0,0):
                                skip=False
                    if skip:
                        continue
                    #if not record.info["SAMPLE"].upper() in includedSamples.keys():
                    #    continue
                except Exception as e:
                    pass
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
                    CNs.append(int(a[3:-1]))
                if len(CNs)==0:
                    continue
                if len(CNs)>1:
                    Type="CNV"
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

def parse_my(filename,contigs=None, samples=None):
    if filename[-4:]==".vcf" or filename[-4:]==".bcf":
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