'''
Created on Nov 23, 2015

@author: Liu
'''
import sys

class interval:
    
    def __init__(self, chrom, start, end, type):  
        self.chrom = chrom 
        self.start = start  
        self.end = end
        self.type = type
        
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

            
def calculate_sensitivity_by_percent(goldstandard, callset, percent):
    count = 0
    for interval1 in goldstandard:
        tag = False    
        for interval2 in callset:
            if interval1.is_overlap_by_percent(interval2, percent):
                tag = True
        if tag == True:
            count = count + 1
    print('Count=' + str(len(callset)))        
    print('Sensitivity=' + str(float(count) / float(len(goldstandard))))

def calculate_fdr_by_percent(goldstandard, callset, percent):
    count = 0
    for interval1 in callset:
        tag = False    
        for interval2 in goldstandard:
            if interval1.is_overlap_by_percent(interval2, percent):
                tag = True
        if tag == True:
            count = count + 1
    print('Count=' + str(len(callset)))        
    print('FDR=' + str(1-float(count) / float(len(callset))))
    

