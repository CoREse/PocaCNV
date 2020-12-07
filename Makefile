CC=g++
AR=ar
CPPFLAGS= -Wall -O3 -Ihtslib
LDFLAGS=-lz -lm -lbz2 -llzma -lpthread
LIBS=
PYTHON=python3.8
INCLUDE=/usr/include/$(PYTHON)

PROJECT_OBJS=jc.o
PROJECT_HEADERS=
EXAMPLE_OBJS=
HTSLIB=htslib/libhts.a

.PHONY: all test clean

all: CGetRDScores.so

CGetRDScores.so: getRDScores.cpp
	$(CC) $^ -o $@ -fPIC -fopenmp -lgsl -lgslcblas -shared -o CGetRDScores.so -I$(INCLUDE) -l$(PYTHON)

$(DNASEQ):$(DNASEQ_OBJS)
	$(AR) -rc $@ $(DNASEQ_OBJS)

example: $(EXAMPLE_OBJS) $(DNASEQ)
	$(LINK.cpp) $^ -o $@

otest: $(PROJECT_OBJS) $(HTSLIB)
	$(LINK.cpp) -o $@ $^ $(LDFLAGS)

test:
	python3 test.py
oldtest:
	time python3 jc.py ~/data/Homo_sapiens_assembly38.fasta.gz ~/data/HG0051*chr22.cram > test/512-4_chr22_testcall.txt
	python3 benchmark.py

fasttest:
	time python3 jc.py ~/data/hg38.fa ~/data/51*_ch1_2M.bam > test/fasttest.txt

qtest:
	time python3 jcrd.py -C chr22 -T ~/data/GRCh38_full_analysis_set_plus_decoy_hla.fa data/rd*chr22.cram.rdf >data/qtest.vcf
qtest100:
	time python3 jcrd.py -C chr22 -T ~/data/GRCh38_full_analysis_set_plus_decoy_hla.fa data/100rd*chr22.cram.rdf -WS 10000 >data/qtest100.vcf
qtest100t1:
	time python3 jcrd.py -J 1 -C chr22 -T ~/data/GRCh38_full_analysis_set_plus_decoy_hla.fa data/100rd*chr22.cram.rdf -WS 10000 >data/qtest100.vcf
qtestco:
	time python3 jcrd.py -C chr22 -T ~/data/GRCh38_full_analysis_set_plus_decoy_hla.fa -LC data/AllData.data data/rd*chr22.cram.rdf >data/qtest.vcf
qtestc:
	time python3 jcrd.py -T ~/data/GRCh38_full_analysis_set_plus_decoy_hla.fa data/test*.cram -W >data/qtestc.vcf
wotest:
	time python3 jcrd.py -T ~/data/GRCh38_full_analysis_set_plus_decoy_hla.fa data/test*.cram -WO
qtest108:
	time python3 -u jcrd.py -T ~/data/hs37d5.fa.gz -WS 2000 ~/data/CHS/2krd/*.rdf -C 22 > data/qtest108.vcf
	python3 benchmark.py -C 22 data/qtest108.vcf
test22:
	time python3 -u jcrd.py -T ~/data/hs37d5.fa.gz -WS 100 ~/data/CHS/*.rdf -C 22 > data/test22.vcf
	python3 benchmark.py -C 22 data/test22.vcf
test1:
	time python3 -u jcrd.py -T ~/data/hs37d5.fa.gz -WS 100 ~/data/CHS/*.rdf -C 1 > data/test1.vcf
	python3 benchmark.py -G ~/data/1000gp/chr1_indel_sv_chs.vcf -C 22 data/test1.vcf
debug:
	bash debugs/debug.sh
bench:
	python3 benchmark.py

$(HTSLIB): htslib/*
	cd htslib && make

clean:
	rm *.o test
