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
LAUNCHER=./launcher

.PHONY: all test clean

all: build

CGetRDScores.so: getRDScores.cpp
	$(CC) $^ -o $@ -O3 -fPIC -fopenmp -lgsl -lgslcblas -shared -o CGetRDScores.so -I$(INCLUDE) -l$(PYTHON)

jcinputraw: jcinput.cpp
	$(CC) $^ -o $@ -g -pg -fopenmp -Lhdf5/lib -lhdf5 -Lhtslib -lhts -lz --std=c++11

$(DNASEQ):$(DNASEQ_OBJS)
	$(AR) -rc $@ $(DNASEQ_OBJS)

example: $(EXAMPLE_OBJS) $(DNASEQ)
	$(LINK.cpp) $^ -o $@

otest: $(PROJECT_OBJS) $(HTSLIB)
	$(LINK.cpp) -o $@ $^ $(LDFLAGS)

cython: rdprocessing.pyx
	python3 cythonsetup.py build_ext --inplace
build: cython CGetRDScores.so jcinputraw

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
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 2000 ~/data/0/CHS/2krd/*.rdf -C 22 > data/qtest108.vcf
	python3 benchmark.py -G ~/data/0/1000gp/chr22_indel_sv_chs.vcf -C 22 data/qtest108.vcf
test2122: build
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CHS*.cram.sd -C 22 -C 21 > data/test2122.vcf
	python3 benchmark.py -G ~/data/0/1000gp/chr22_indel_sv_chs.vcf -C 22 data/test2122.vcf
test22: build
	#time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CHS*.cram.hdf5 -C 22 > data/test22.vcf
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 ~/data/CHS/crams/*CHS*.cram.hdf5 -C 22 > data/test22.vcf
	python3 benchmark.py -G ~/data/0/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf -C 22 -SF /data/0/cre/CHS/samples.txt data/test22.vcf
testCDX22: build
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CDX*.cram.sd -C 22 > data/testCDX22.vcf
	python3 benchmark.py -G ~/data/0/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf -C 22 -SF /home/cre/data/CDX/samples.txt data/testCDX22.vcf
testCDX10: build
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CDX*.cram.sd -C 10 > data/testCDX10.vcf
	python3 benchmark.py -G ~/data/0/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf -C 10 -SF /home/cre/data/CDX/samples.txt data/testCDX10.vcf
test22light: build
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/HG0040*CHS*.cram.hdf5 -C 22 > data/test22light.vcf
	python3 benchmark.py -G ~/data/0/1000gp/chr22_indel_sv_chs.vcf -C 22 data/test22light.vcf
test22tg: build
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*.cram.sd -C 22 > data/test22tg.vcf
	python3 benchmark.py -G ~/data/0/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf -C 22 -SF /data/0/cre/CHS/samples.txt data/test22tg.vcf
test1: cython
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CHS*.cram.hdf5 -C 1 > data/test1.vcf
	python3 benchmark.py -G ~/data/0/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf -C 1 -SF /data/0/cre/CHS/samples.txt data/test1.vcf
testX:
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CHS*.cram.sd -J 16 -C X > data/testX.vcf
	python3 benchmark.py -G ~/data/0/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -C X data/testX.vcf
testY:
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 ~/data/0/CHS/*.rdf -C Y > data/testY.vcf
	python3 benchmark.py -G ~/data/0/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz -C Y data/testY.vcf
chsedata: cython
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CHS*.cram.sd -EN CHS > data/edata.vcf
clmedata: cython
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CLM*.cram.sd -EN CLM > data/clmedata.vcf 2> data/clmedata.log
cdxedata: cython
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CDX*.cram.sd -EN CDX > data/cdxedata.vcf 2> data/cdxedata.log
clmwhole: cython
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 ~/data/CLM/*CLM*.cram -EN CLM > data/clmwhole.vcf 2> data/clmwhole.log
cdxwhole: cython
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 ~/data/CDX/*CDX*.cram -EN CDX > data/cdxwhole.vcf 2> data/cdxwhole.log
3g: build
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CHS*.cram.sd -EN CHS > data/test3gCHS.vcf 2> logs/test3g.log
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CDX*.cram.sd -EN CDX > data/test3gCDX.vcf 2>> logs/test3g.log
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CLM*.cram.sd -EN CLM > data/test3gCLM.vcf 2>> logs/test3g.log
2g: build
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CDX*.cram.sd -EN CDX > data/test3gCDX.vcf 2>> logs/test3g.log
	time $(LAUNCHER) -T ~/data/0/hs37d5.fa.gz -WS 100 data/*CLM*.cram.sd -EN CLM > data/test3gCLM.vcf 2>> logs/test3g.log
debug:
	bash debugs/debug.sh
bench:
	python3 benchmark.py

$(HTSLIB): htslib/*
	cd htslib && make

clean:
	rm *.o test
