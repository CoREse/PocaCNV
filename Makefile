CC=g++
AR=ar
CPPFLAGS= -Wall -O3 -Ihtslib
LDFLAGS=-lz -lm -lbz2 -llzma -lpthread
LIBS=
PYTHON=python3.6m
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
qtestc:
	time python3 jcrd.py -T ~/data/GRCh38_full_analysis_set_plus_decoy_hla.fa data/test*.cram -W >data/qtestc.vcf
bench:
	python3 benchmark.py

$(HTSLIB): htslib/*
	cd htslib && make

clean:
	rm *.o test
