CC=g++
AR=ar
CPPFLAGS= -Wall -O3 -Ihtslib
LDFLAGS=-lz -lm -lbz2 -llzma -lpthread
LIBS=

PROJECT_OBJS=jc.o
PROJECT_HEADERS=
EXAMPLE_OBJS=
HTSLIB=htslib/libhts.a

.PHONY: all test clean

all: test

$(DNASEQ):$(DNASEQ_OBJS)
	$(AR) -rc $@ $(DNASEQ_OBJS)

example: $(EXAMPLE_OBJS) $(DNASEQ)
	$(LINK.cpp) $^ -o $@

otest: $(PROJECT_OBJS) $(HTSLIB)
	$(LINK.cpp) -o $@ $^ $(LDFLAGS)

test:
	time python3 jc.py ~/data/Homo_sapiens_assembly38.fasta.gz ~/data/HG0051*chr22.cram > test/512-4_chr22_testcall.txt
	python3 benchmark.py

$(HTSLIB): htslib/*
	cd htslib && make

clean:
	rm *.o test
