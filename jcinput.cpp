#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <vector>
#include <string>
#include "hdf5-1.12.0/src/hdf5.h"
#include <mutex>
#include <unistd.h>
//#include <omp.h>
#include "htslib/htslib/bgzf.h"
#include <thread>

using namespace std;

const int ReadThreadN=1;

float MedianInsertionSize=481.2;
float ISSD=36.4;
//bool DynamicDRP=false;
//int DynamicCount=1000;
int RDWindowSize=100;

mutex TopLock, WriteLock, StdinLock;

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
		init();
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

inline Contig * genCopyOfContigs(const Contig * ContigModels, const int NSeq)
{
	Contig * Contigs=(Contig*) malloc(sizeof(Contig)*(NSeq));
	for (int i=0;i<NSeq;++i)
	{
		new (Contigs+i) Contig(ContigModels[i].Size,ContigModels[i].Name.c_str());
	}
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

void saveHDF5(const char * FileName, const char * SampleName, const Contig * Contigs, int NSeq ,int ReadCount, int RDWindowSize)
{
	hid_t       file_id, group_id;  // identifiers
	hid_t       aid, atype, attr;
	herr_t      status;

	
	aid = H5Screate(H5S_SCALAR);
	atype = H5Tcopy (H5T_C_S1);
	//status = H5Tset_cset(atype,H5T_CSET_UTF8);
	status = H5Tset_size (atype, strlen(SampleName)+1);

	// Create a new file using default properties. 
	file_id = H5Fcreate(FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	group_id = H5Gopen (file_id, "/", H5P_DEFAULT);
	attr = H5Acreate (group_id, "SampleName", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite (attr, atype, SampleName);
	status = H5Aclose (attr);
	
	atype = H5Tcopy (H5T_NATIVE_INT32);
	attr = H5Acreate (group_id, "SampleReadCount", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite (attr, atype, &ReadCount);
	status = H5Aclose (attr);

	atype = H5Tcopy (H5T_NATIVE_INT32);
	attr = H5Acreate (group_id, "RDWindowSize", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite (attr, atype, &RDWindowSize);
	status = H5Aclose (attr);

	atype = H5Tcopy (H5T_NATIVE_FLOAT);
	attr = H5Acreate (group_id, "ISAV", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite (attr, atype, &MedianInsertionSize);
	status = H5Aclose (attr);

	atype = H5Tcopy (H5T_NATIVE_FLOAT);
	attr = H5Acreate (group_id, "ISSD", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Awrite (attr, atype, &ISSD);
	status = H5Aclose (attr);

	status = H5Gclose(group_id);
	status = H5Sclose (aid);

	for (int i=0;i<NSeq;++i)
	{
		group_id = H5Gcreate (file_id,Contigs[i].Name.c_str(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		
		hsize_t dim[1];
		dim[0]=Contigs[i].Size;
		aid = H5Screate_simple(1,dim,NULL);
		atype = H5Tcopy (H5T_IEEE_F32BE);

		attr = H5Dcreate (group_id, "RDWindows", atype, aid, H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
		status = H5Dwrite (attr, H5T_NATIVE_FLOAT, H5S_ALL,H5S_ALL,H5P_DEFAULT, Contigs[i].RDWindows);

		status = H5Sclose (aid);
		status = H5Dclose (attr);

		aid = H5Screate_simple(1,dim,NULL);
		atype = H5Tcopy (H5T_IEEE_F32BE);

		attr = H5Dcreate (group_id, "AverageClipLengths", atype, aid, H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
		status = H5Dwrite (attr, H5T_NATIVE_FLOAT, H5S_ALL,H5S_ALL,H5P_DEFAULT, Contigs[i].AverageClipLengths);

		status = H5Sclose (aid);
		status = H5Dclose (attr);

		hsize_t dim2[2];
		dim2[0]=Contigs[i].DRPs.size();
		dim2[1]=5;

		//int **data=(int**)malloc(dim2[0]*sizeof(int*));
		//for (int j=0;j<dim2[0];++j) data[j]=(int*)malloc(dim2[1]*sizeof(int));
		//int data[dim2[0]][dim2[1]];
		//int *[5]data=(int *[5])malloc(sizeof(int)*dim2[0]*dim2[1]);
		DRP* data=(DRP *)malloc(sizeof(DRP)*dim2[0]);
		_List_const_iterator<DRP> it=Contigs[i].DRPs.begin();
		for (int j=0;j<Contigs[i].DRPs.size();++j)
		{
			data[j]=*it;
			/*
			data[j][0]=it->Start;
			data[j][1]=it->End;
			data[j][2]=it->InnerStart;
			data[j][3]=it->InnerEnd;
			data[j][4]=it->SVT;*/
			it=next(it);
		}

		aid = H5Screate_simple(2,dim2,NULL);
		atype = H5Tcopy (H5T_NATIVE_INT32);

		attr = H5Dcreate (group_id, "DRPs", atype, aid, H5P_DEFAULT, H5P_DEFAULT,H5P_DEFAULT);
		status = H5Dwrite (attr, atype, H5S_ALL,H5S_ALL,H5P_DEFAULT, data);

		status = H5Sclose (aid);
		status = H5Dclose (attr);

		aid = H5Screate(H5S_SCALAR);
		atype = H5Tcopy (H5T_NATIVE_INT32);

		attr = H5Acreate (group_id, "ContigSampleReadCounts", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite (attr, atype, &(Contigs[i].ContigReadCount));
		status = H5Aclose (attr);

		atype = H5Tcopy (H5T_NATIVE_INT32);

		attr = H5Acreate (group_id, "Length", atype, aid, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite (attr, atype, &(Contigs[i].Size));
		status = H5Aclose (attr);

		status = H5Sclose (aid);

		status = H5Gclose(group_id);
		//for (int j=0;j<dim2[0];++j) free(data[j]);
		free(data);
	}
		
	// Close the file. 
	status = H5Fclose(file_id);
}

struct BrBlock
{
	bam1_t * Block;
	int Size;
	int Count;
	bool Mature;
	bool Over;//read function will create a Over==true block as a last sentinel
	BrBlock* Next;
	BrBlock(int Size=10000)
	:Size(Size),Count(0),Over(false),Next(NULL),Mature(false)
	{
		if (Size!=0) Block=(bam1_t*)calloc(Size,sizeof(bam1_t));
	}
	~BrBlock()
	{
		freeBlock();
	}
	void freeBlock()
	{
		if (Block!=NULL) free(Block);
		Block=NULL;
	}
	private:
	//prohibit copy
	BrBlock(const BrBlock&);
	BrBlock& operator=(const BrBlock&);
};

void handlebr(bam1_t *br, Contig * Contigs, int * CordinTrans, int& ReadCount, int & UnmappedCount)
{
		WriteLock.lock();
		++ReadCount;
		WriteLock.unlock();
		if (read_is_unmapped(br))
		{
			WriteLock.lock();
			++UnmappedCount;
			WriteLock.unlock();
		}
		else
		{
			Contig &TheContig=Contigs[CordinTrans[br->core.tid]];
			int End=bam_endpos(br);
			if (read_is_paired(br) && (!read_mate_is_unmapped(br)) && br->core.isize>0)
			{
				//TODO: change DRP calc way.
        		int isize=br->core.isize;
				if (isize> MedianInsertionSize+3*ISSD)
				{
					WriteLock.lock();
					TheContig.DRPs.push_back(DRP(br->core.pos,br->core.pos+br->core.isize,End,br->core.mpos,0));
					WriteLock.unlock();
				}
				else if (isize< MedianInsertionSize-3*ISSD)
				{
					WriteLock.lock();
					TheContig.DRPs.push_back(DRP(br->core.pos,br->core.pos+br->core.isize,End,br->core.mpos,1));
					WriteLock.unlock();
				}
			}
			int Index=(br->core.pos+End)/2/RDWindowSize;
			float ClipLength=getClipLength(br);
			WriteLock.lock();
			++TheContig.RDWindows[Index];
			TheContig.AverageClipLengths[Index]+=ClipLength;
			++TheContig.ContigReadCount;
			WriteLock.unlock();
		}
}

void handleBrBlock(BrBlock * TheBlock, Contig * Contigs, int * CordinTrans, int& ReadCount, int & UnmappedCount)
{
	int BlockReadCount=0, BlockUnmappedCount=0;
	int BrReadCount=0,BrUnmappedCount=0;
	for (int i=0;i<TheBlock->Count;++i)
	{
		BrReadCount=0;
		BrUnmappedCount=0;
		handlebr(TheBlock->Block+i,Contigs,CordinTrans,BrReadCount,BrUnmappedCount);
		BlockReadCount+=BrReadCount;
		BlockUnmappedCount+=BrUnmappedCount;
	}
	WriteLock.lock();
	ReadCount+=BlockReadCount;
	UnmappedCount+=BlockUnmappedCount;
	WriteLock.unlock();
}

void getBlockAndProcess(BrBlock ** UpMostBlock, Contig * Contigs, int * CordinTrans, int& ReadCount, int & UnmappedCount)
{
	while (true)
	{
		if ((*UpMostBlock)->Over) break;
		TopLock.lock();
		if (*UpMostBlock==NULL || !((*UpMostBlock)->Mature))
		{
			TopLock.unlock();
			sleep(1);
		}
		else
		{
			BrBlock * TheBlock=*UpMostBlock;
			if (TheBlock->Over)
			{
				TopLock.unlock();
				break;
			}
			else *UpMostBlock=(*UpMostBlock)->Next;
			TopLock.unlock();
			handleBrBlock(TheBlock,Contigs,CordinTrans,ReadCount,UnmappedCount);
			delete TheBlock;
		}
	}
}

void readBamToBrBlock(htsFile * SamFile,bam_hdr_t *Header, BrBlock** Top)
{
	BrBlock * TheBlock=new BrBlock();
	TheBlock->Next=new BrBlock();
	while (sam_read1(SamFile, Header, TheBlock->Block+TheBlock->Count) >=0)
	{
		++TheBlock->Count;
		if (TheBlock->Count>=TheBlock->Size)
		{
			TheBlock->Mature=true;
			if (*Top==NULL)//The only one who could fill NULL top
			{
				*Top=TheBlock;
			}
			TheBlock=TheBlock->Next;
			TheBlock->Next=new BrBlock();
		}
	}
	TheBlock->Mature=true;
	TheBlock->Next->Over=true;
}

void takeStdinAndHandleBr(int &ReadCount, int &UnmappedCount, Contig * Contigs, int * CordinTrans, bam_hdr_t * Header)
{
    bam1_t *br=bam_init1();
	size_t linebuffersize=1024*0124;
	char * linebuffer=(char*) malloc(linebuffersize);
	StdinLock.lock();
	int length=getline(&linebuffer,&linebuffersize,stdin);
	StdinLock.unlock();
	while (length>=0)
	{
		kstring_t ks={length,linebuffersize,linebuffer};
		sam_parse1(&ks,Header,br);
		handlebr(br,Contigs,CordinTrans,ReadCount,UnmappedCount);
		StdinLock.lock();
		length=getline(&linebuffer,&linebuffersize,stdin);
		StdinLock.unlock();
	}
	free(linebuffer);
	bam_destroy1(br);
}

void readSam(const Contig * ContigModels, const int NSeq, const char * ReferenceFileName, const char * SampleFileName, bool UseStdin=false, const char * datadir="data")
{
	fprintf(stderr,"Reading %s...\n",SampleFileName);
	char HDF5FileName[strlen(SampleFileName)+strlen(datadir)+10];
	strcpy(HDF5FileName,SampleFileName);
	int pathloc=strlen(SampleFileName)-1;
	for (;pathloc>=0;--pathloc)
	{
		if (SampleFileName[pathloc]=='/') break;
	}
	strcpy(HDF5FileName,datadir);
	if (HDF5FileName[strlen(HDF5FileName)-1]!='/') strcat(HDF5FileName,"/");
	strcat(HDF5FileName,SampleFileName+pathloc+1);
	strcat(HDF5FileName,".hdf5");
	htsFile* SamFile;//the file of BAM/CRAM
	bam_hdr_t* Header;//header for BAM/CRAM file
    SamFile = hts_open(SampleFileName, "rb");

	//set reference file
	if (NULL != ReferenceFileName)
	{
		char referenceFilenameIndex[128];
		strcpy(referenceFilenameIndex, ReferenceFileName);
		strcat(referenceFilenameIndex, ".fai");
		int ret = hts_set_fai_filename(SamFile, referenceFilenameIndex);
	}
	Header = sam_hdr_read(SamFile);
	char SampleName[1024];
	kstring_t ks=KS_INITIALIZE;
	ks_resize(&ks,1000);
	sam_hdr_find_line_pos(Header,"RG",0,&ks);
	int splitn;
	int * splitoffsets=ksplit(&ks,'\t',&splitn);
	for (int i=0;i<splitn;++i)
	{
		if (memcmp(ks.s+splitoffsets[i],"SM",2)==0)
		{
			int End=ks_len(&ks);
			if (i<splitn-1) End=splitoffsets[i+1];
			memcpy(SampleName,ks.s+splitoffsets[i]+3,End-splitoffsets[i]-3);
			SampleName[End-splitoffsets[i]-3]='\0';
		}
	}
	ks_free(&ks);
	Contig * Contigs=genCopyOfContigs(ContigModels,NSeq);
	int CordinTrans[Header->n_targets];
	for (int i=0;i<Header->n_targets;++i)
	{
		bool NoFound=true;
		for (int j=0;j<NSeq;++j)
		{
			if (strcmp(Contigs[j].Name.c_str(),Header->target_name[i])==0)
			{
				CordinTrans[i]=j;
				NoFound=false;
				break;
			}
		}
		if (NoFound) fprintf(stderr,"[WARN] There's no contig in reference named %s.\n",Header->target_name[i]);
	}

	int ReadCount=0, UnmappedCount=0;
	if (UseStdin)
	{
		if (ReadThreadN==1) takeStdinAndHandleBr(ReadCount,UnmappedCount,Contigs,CordinTrans,Header);
		else
		{
			thread ** Threads=(thread **)malloc(sizeof(thread)*ReadThreadN);
			for (int i =0;i<ReadThreadN;++i)
			{
				Threads[i]=new thread(takeStdinAndHandleBr,ref(ReadCount),ref(UnmappedCount),Contigs,(int*)CordinTrans,Header);
			}
			for (int i=0;i<ReadThreadN;++i)
			{
				Threads[i]->join();
			}
			
			for (int i=0;i<ReadThreadN;++i)
			{
				delete Threads[i];
			}
			free(Threads);
		}
	}
	else
	{
    	bam1_t *br=bam_init1();
		while(sam_read1(SamFile, Header, br) >=0)//read record
		{
			handlebr(br,Contigs,CordinTrans,ReadCount,UnmappedCount);
		}
		bam_destroy1(br);
	}
	/*
	BrBlock * Top[1];
	*Top=NULL;
	readBamToBrBlock(SamFile,Header,Top);
	getBlockAndProcess(Top,Contigs,CordinTrans,ReadCount,UnmappedCount);
	delete *Top;*/

	for (int i=0;i<NSeq;++i)
	{
		Contigs[i].DRPs.sort(compareDRP);
		for (int j=0;j<Contigs[i].Size;++j)
		{
			if (Contigs[i].RDWindows[j]!=0) Contigs[i].AverageClipLengths[j]/=Contigs[i].RDWindows[j];
		}
	}

	saveHDF5(HDF5FileName, SampleName, Contigs, NSeq, ReadCount,RDWindowSize);

	for (int i=0;i<NSeq;++i)
	{
		Contigs[i].~Contig();
	}
	free(Contigs);
	hts_close(SamFile);
	bam_hdr_destroy(Header);
}

int main(int argc, char* argv[])
{
	float AISAV=-1, AISSD=-1;
	AISAV=atoi(argv[1]);
	AISSD=atoi(argv[2]);
	if (AISAV>0) MedianInsertionSize=AISAV;
	if (AISSD>0) ISSD=AISSD;
    char * ReferenceFilename=argv[3];
	int NSeq;
	Contig * Contigs=getContigs(ReferenceFilename,&NSeq,RDWindowSize);
    //omp_set_num_threads(2);
    //#pragma omp parallel for
	for (int i=4;i<argc;++i)
	{
		readSam(Contigs,NSeq,ReferenceFilename,argv[i],1);
	}
	
	free(Contigs);
    return 0;
}