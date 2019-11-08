import benchmark
import sys

def parse_goldstandard(input_file, min_length=0):
    records = []
    f = open(input_file, 'r')
    for line in f:
        record = line.split()
        chrom='chr'+record[0]
        tmp = benchmark.interval(chrom, int(record[1]), int(record[2]), record[3])
        if record[3]=='LOSS' and (int(record[2])-int(record[1])+1) >= min_length:
            records.append(tmp)
    f.close()
    return records

def parse_triocnv_rd(sample_index, input_file, min_length=0, max_length=sys.maxsize):
    records = []
    f = open(input_file, 'r')
    f.readline()
    for line in f:
        record = line.split()
        cn = int(record[sample_index+3])
        if cn > 2:
            tmp = benchmark.interval(record[0], int(record[1]), int(record[2]), 'GAIN')
            if int(record[3]) >= min_length and int(record[3]) < max_length:
                records.append(tmp)
        if cn < 2:
            tmp = benchmark.interval(record[0], int(record[1]), int(record[2]), 'LOSS')
            if int(record[3]) >= min_length and int(record[3]) < max_length:
                records.append(tmp)
    f.close()
    return records


def parse_triocnv_drp(sample_index, input_file):
    records = []
    f = open(input_file, 'r')
    f.readline()
    for line in f:
        record = line.split()
        num_drps = int(record[sample_index + 5])
        if record[4] == 'DELETION':
            tmp = benchmark.interval(record[0], int(record[1]), int(record[2]), 'LOSS')
            if num_drps > 0 :
                records.append(tmp)
        if record[4] == 'TANDEM_DUPLICATION':
            tmp = benchmark.interval(record[0], int(record[1]), int(record[2]), 'GAIN')
            if num_drps > 0 :
                records.append(tmp)
    f.close()
    return records

def parse_integrated_sv_map(sv_file,sample,output_file):
    records = []
    sample_index=-1
    f = open(sv_file, 'r')
    output = open(output_file, 'w')
    for line in f:
        if(line.startswith('#CHROM')):
            header_record=line.split()
            sample_index=header_record.index(sample)
        if(not line.startswith('#')):
            record=line.split()
            if(record[0]=='chrX'):
                break;
            if(record[sample_index]!='0|0'):
                chrom=record[0]
                start=int(record[1])
                end=0
                type=''
                info = record[7].split(";");
                for item in info:
                    if item.startswith('SVTYPE'):
                        type=item[7:]
                    if item.startswith('END'):
                        end=int(item[4:])
                if(type=='DEL'):
                    tmp = benchmark.interval(chrom, start, end, 'LOSS')
                    output.write(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+'LOSS'+'\n')
                    records.append(tmp)
                if(type=='DUP'):
                    tmp = benchmark.interval(chrom, start, end, 'GAIN')
                    output.write(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+'GAIN'+'\n')
                    records.append(tmp)
    f.close()
    output.close()
    return records
            
if __name__ == '__main__':
#     parse_integrated_sv_map('D:/My_Project/TrioCNV2/ceutrio/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf','NA12878'
#                             ,'D:/My_Project/TrioCNV2/ceutrio/NA12878.txt')
    goldstandard = parse_goldstandard('D:/My_Project/TrioCNV2/ceutrio/NA12878.txt',100)
    triocnv_rd=parse_triocnv_rd(2, 'D:/My_Project/TrioCNV2/ceutrio/TrioCNV_RD_Chr1.txt')
    triocnv_drp=parse_triocnv_drp(2, 'D:/My_Project/TrioCNV2/ceutrio/TrioCNV_DRP.txt')
    cnvnator=benchmark.parse_cnvnator('D:/My_Project/TrioCNV2/ceutrio/NA12878.cnvnator.txt')
    benchmark.calculate_sensitivity_by_percent(goldstandard,cnvnator,0.5)

