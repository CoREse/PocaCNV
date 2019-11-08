import pysam
import sys

vcf=pysam.VariantFile(sys.argv[1],"r")
svn=0
#print(vcf.header)
i=0
#print(vcf.header)
vcf.subset_samples(["UW_PacBio_YRI_NA19240"])
for record in vcf.fetch():
    record: pysam.libcbcf.VariantRecord
    print(record)
    break
    maxdiffn=0
    try:
        if record.info["SVTYPE"]!="DEL" and record.info["SVTYPE"]!="INS" and record.info["SVTYPE"]!="NONE" and record.info["SVTYPE"]=="DUP":
            print(record)
            break
        else:
            continue
    except:
        continue
    if i<10:
        print(record)
        try:
            print(type(record.info["SVLEN"]))
        except:
            print("no 1234")
        i+=1
        continue
    break
    for a in record.alts:
        maxdiffn=max(maxdiffn,abs(record.rlen-len(a)))
    if record.alts[0]=='<INS>':
        if svn<10:
            print(record)
            print(record.info["SVCLASS"])
        #print(record)
        svn+=1
print(svn)