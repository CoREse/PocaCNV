#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <vector>
#include <string>

using namespace std;

typedef struct Bam_file{
	bool _is_record_set;//set to 1 when _brec has data
	htsFile* _hfp;//the file of BAM/CRAM
	bam_hdr_t* _hdr;//header for BAM/CRAM file
	hts_idx_t* _hidx;//index for bam/cram file
	hts_itr_t* _hitr;//Iterator for bam/cram file
	bam1_t _brec;//current BAM record

	// track for debug only:
	unsigned _record_no;
	char * _stream_name;
	bool _is_region;
	char* _region;
}Bam_file;

void bam_file_close(Bam_file * bf){
	hts_close(bf->_hfp);
	bam_hdr_destroy(bf->_hdr);
	//hts_idx_destroy(bf->_hidx);
	//hts_itr_destroy(bf->_hitr);
	//if(bf->_stream_name) free(bf->_stream_name);
	//if(bf->_region) free(bf->_region);
	memset(bf, 0, sizeof(Bam_file));
}

struct DRP
{
	int Start,End,InnerStart,InnerEnd,SVT;
	DRP(int start, int end, int innerstart, int innerend, int svt)
	:Start(start),End(end),InnerStart(innerstart),InnerEnd(innerend),SVT(svt)
	{}
};

struct Contig
{
	const unsigned int Size;
	float *RDWindows;
	float *AverageClipLengths;
	unsigned int ContigReadCount;
	list<DRP> DRPs;
	string Name;
	void init()
	{
		for (int i=0;i<Size;++i)
		{
			RDWindows[i]=0;
			AverageClipLengths[i]=0;
		}
	}
	Contig(int Size, const char * SName)
	:Size(Size),ContigReadCount(0),DRPs(),Name(SName)
	{
		RDWindows=(float*) malloc(sizeof(float)*Size);
		AverageClipLengths=(float*) malloc(sizeof(float)*Size);
	}
	~Contig()
	{
		free(RDWindows);
		free(AverageClipLengths);
		RDWindows=NULL;
		AverageClipLengths=NULL;
	}
	void clear()
	{
		DRPs.clear();
		init();
	}
	Contig* genNew()
	{
		return new Contig(Size,Name.c_str());
	}
	private:
	//prohibit copy
	Contig(const Contig&);
	Contig& operator=(const Contig&);
};

Contig * getContigs(const char * ReferenceFN, int * NSeq, int WindowSize=100)
{
	faidx_t * Ref=fai_load(ReferenceFN);
	*NSeq=faidx_nseq(Ref);
	Contig * Contigs=(Contig*) malloc(sizeof(Contig)*(*NSeq));
	for (int i=0;i<(*NSeq);++i)
	{
		const char * ContigName=faidx_iseq(Ref,i);
		int SeqLen=faidx_seq_len(Ref,ContigName);
		int Size=SeqLen/WindowSize;
		if (SeqLen%WindowSize!=0) ++Size;
		new (Contigs+i) Contig(Size,ContigName);
	}
	fai_destroy(Ref);
	return Contigs;
}

bool compareDRP(const DRP& a, const DRP& b)
{
	return a.Start<b.Start;
}

#define read_is_unmapped(b) (((b)->core.flag&BAM_FUNMAP) != 0)
#define read_mate_is_unmapped(b) (((b)->core.flag&BAM_FMUNMAP) != 0)
#define read_is_paired(b) (((b)->core.flag&BAM_FPAIRED) != 0)
#define read_is_read1(b) (((b)->core.flag&BAM_FREAD1) != 0)
#define read_is_read2(b) (((b)->core.flag&BAM_FREAD2) != 0)

float getClipLength(bam1_t * br)
{
	uint32_t* cigars=bam_get_cigar(br);
	float ClipLength=0;
	for (int i=0;i<br->core.n_cigar;++i)
	{
		if ((cigars[i]&0xf)==4 || (cigars[i]&0xf)==5) ClipLength+=cigars[i]>>4;
	}
	return ClipLength;
}

int MedianInsertionSize=550;
int ISSD=150;
int RDWindowSize=100;

int main(int argc, char* argv[])
{
    char * referenceFilename=argv[1];
    char * filename=argv[2];
	int NSeq;
	Contig * Contigs=getContigs(referenceFilename,&NSeq,RDWindowSize);
    Bam_file *bf=new Bam_file;
    bf->_hfp = hts_open(filename, "rb");

	//set reference file
	if (NULL != referenceFilename)
	{
		char referenceFilenameIndex[128];
		strcpy(referenceFilenameIndex, referenceFilename);
		strcat(referenceFilenameIndex, ".fai");
		int ret = hts_set_fai_filename(bf->_hfp, referenceFilenameIndex);
	}
	bf->_hdr = sam_hdr_read(bf->_hfp);
	int CordinTrans[bf->_hdr->n_targets];
	for (int i=0;i<bf->_hdr->n_targets;++i)
	{
		bool NoFound=true;
		for (int j=0;j<NSeq;++j)
		{
			if (strcmp(Contigs[j].Name.c_str(),bf->_hdr->target_name[i])==0)
			{
				CordinTrans[i]=j;
				NoFound=false;
				break;
			}
		}
		if (NoFound) fprintf(stderr,"[WARN] There's no contig in reference named %s.\n",bf->_hdr->target_name[i]);
	}

    bam1_t *br=bam_init1();
	int ReadCount=0, UnmappedCount=0;
	while(sam_read1(bf->_hfp, bf->_hdr, br) >=0)//read record
	{
		++ReadCount;
		if (read_is_unmapped(br))
		{
			++UnmappedCount;
		}
		else
		{
			Contig &TheContig=Contigs[CordinTrans[br->core.tid]];
			int End=bam_endpos(br);
			if (read_is_paired(br) && read_is_read1(br) && (!read_mate_is_unmapped(br)))
			{
				//TODO: change DRP calc way.
        		int isize=abs(br->core.isize);
				if (isize> MedianInsertionSize+3*ISSD)
				{
					TheContig.DRPs.push_back(DRP(br->core.pos,br->core.pos+br->core.isize,End,br->core.mpos,0));
				}
				else if (isize< MedianInsertionSize-3*ISSD)
				{
					TheContig.DRPs.push_back(DRP(br->core.pos,br->core.pos+br->core.isize,End,br->core.mpos,1));
				}
			}
			int Index=(br->core.pos+End)/2/RDWindowSize;
			++TheContig.RDWindows[Index];
			float ClipLength=getClipLength(br);
			TheContig.AverageClipLengths[Index]+=ClipLength;
			++TheContig.ContigReadCount;
		}
	}
    bam_file_close(bf);

	FILE * rdf=fopen("data/jcinput.rd.txt","w");
	FILE * avf=fopen("data/jcinput.av.txt","w");
	FILE * drf=fopen("data/jcinput.dr.txt","w");
	for (int i=0;i<NSeq;++i)
	{
		Contigs[i].DRPs.sort(compareDRP);
		for (auto j=Contigs[i].DRPs.begin();j!=Contigs[i].DRPs.end();++j) fprintf(drf,"%d %d %d %d %d\n",j->Start,j->End,j->InnerStart,j->InnerEnd,j->SVT);
		for (int j=0;j<Contigs[i].Size;++j)
		{
			if (Contigs[i].RDWindows[j]!=0) Contigs[i].AverageClipLengths[j]/=Contigs[i].RDWindows[j];
			fprintf(rdf,"%.6f\n",Contigs[i].RDWindows[j]);
			fprintf(avf,"%.6f\n",Contigs[i].AverageClipLengths[j]);
		}
		break;
	}
	fclose(rdf);
	fclose(avf);
	fclose(drf);

	for (int i=0;i<NSeq;++i)
	{
		Contigs[i].~Contig();
	}
	free(Contigs);
	delete bf;
    return 0;
}