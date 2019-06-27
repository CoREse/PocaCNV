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
            tmp = interval(chrom, start, end, 'LOSS')
            if tmp.get_length() >= min_length and tmp.get_length() < max_length:
                    records.append(tmp)
        if type == 'duplication':
            tmp = interval(chrom, start, end, 'GAIN')
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
        if interval1.type=="GAIN":
            gaintotal+=1
        elif interval1.type=="LOSS":
            losstotal+=1
        elif interval1.type=="INS":
            instotal+=1
        tag = False    
        for interval2 in callset:
            if interval1.is_overlap_by_percent_or_near(interval2, percent,near):
                tag = True
        if tag == True:
            if printmatch:
                print(interval1)
            count = count + 1
            if interval1.type=="GAIN":
                gaincount+=1
            elif interval1.type=="LOSS":
                losscount+=1
            elif interval1.type=="INS":
                inscount+=1
    print('All=%s, Found=%s, Recalled=%s(percentage=%s)' % (len(goldstandard),str(len(callset)),count,percent))
    print('Sensitivity: All:%s, Gain:%s(%s/%s), Loss:%s(%s/%s), Ins:%s(%s/%s)' %(str(float(count) / float(len(goldstandard))),\
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
        if interval1.type=="GAIN":
            gaintotal+=1
        elif interval1.type=="LOSS":
            losstotal+=1
        elif interval1.type=="INS":
            instotal+=1
        tag = False    
        for interval2 in goldstandard:
            if interval1.is_overlap_by_percent_or_near(interval2, percent,near):
                tag = True
        if tag == True:
            if printmatch:
                print(interval1)
            count = count + 1
            if interval1.type=="GAIN":
                gaincount+=1
            elif interval1.type=="LOSS":
                losscount+=1
            elif interval1.type=="INS":
                inscount+=1
    print('All=%s, Found=%s, Correct=%s(percentage=%s)' % (len(goldstandard),str(len(callset)),count,percent))
    print('FDR(Right/All): All:%s, Gain:%s(%s/%s), Loss:%s(%s/%s), Ins:%s(%s/%s)' \
        %(str(1-float(count) / float(len(callset))),1-gaincount/gaintotal if gaintotal !=0 else "N/A",gaincount,gaintotal,\
            1-losscount/losstotal if losstotal!=0 else "N/A",losscount,losstotal,\
                1-inscount/instotal if instotal!=0 else "N/A",inscount,instotal))

def parse_vcf(filename,contigs=None):
    vcf=pysam.VariantFile(filename,"r")
    result=[]
    included={}
    if contigs!=None:
        for c in contigs:
            included[c]=True
    for record in vcf.fetch():
        try:
            if contigs!=None:
                if not included[record.chrom]:
                    continue
            temp=None
            if record.info["SVTYPE"]=="INS":
                temp=interval(record.chrom,record.pos,record.stop+1,"INS")
            elif record.info["SVTYPE"]=="DEL":
                temp=interval(record.chrom,record.pos,record.pos+abs(record.info["SVLEN"]),"LOSS")
            elif record.info["SVTYPE"]=="DUP":
                temp=interval(record.chrom,record.pos,record.pos+abs(record.info["SVLEN"]),"GAIN")
            else:
                continue
            result.append(temp)
        except:
            pass
    vcf.close()
    return result

def parse_my(filename):
    myout=open(filename,"r")
    records=[]
    for line in myout:
        line=line.replace(" ","")
        sl=line.split(",")
        blockstart=int(sl[0].split("-")[0].split(":")[1])
        blockend=int(sl[0].split("-")[1].split(":")[1])
        breakstart=int(sl[-2].split(":")[-1])
        breakend=int(sl[-1].split(":")[1].split("]")[0])
        temp=interval(sl[0].split(":")[0],breakstart,breakend,"LOSS" if sl[1]=="DEL" else "INS")
        records.append(temp)
    myout.close()
    return records