#include <Python.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "include/stats.hpp"
//#include <boost/math/distributions/poisson.hpp>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_cdf.h>

struct Cand
{
    int Size;
    int *EBs;
    int *EEs;
    int *SampleIs;
    int *ECNs;
    double *PassCs;
    double *Confs;
    Cand():Size(0),EBs(NULL),EEs(NULL),SampleIs(NULL),ECNs(NULL),PassCs(NULL),Confs(NULL){}
};

void allocCand(Cand *a)
{
    if (a->Size!=0)
    {
        a->EBs=(int*)malloc(sizeof(int)*a->Size);
        a->EEs=(int*)malloc(sizeof(int)*a->Size);
        a->SampleIs=(int*)malloc(sizeof(int)*a->Size);
        a->ECNs=(int*)malloc(sizeof(int)*a->Size);
        a->PassCs=(double*)malloc(sizeof(double)*a->Size);
        a->Confs=(double*)malloc(sizeof(double)*a->Size);
    }
}
void freeCand(Cand *a)
{
    if (a->EBs!=NULL) free(a->EBs);
    if (a->EEs!=NULL) free(a->EEs);
    if (a->SampleIs!=NULL) free(a->SampleIs);
    if (a->ECNs!=NULL) free(a->ECNs);
    if (a->PassCs!=NULL) free(a->PassCs);
    if (a->Confs!=NULL) free(a->Confs);
}
#define min(x,y) (x>y?y:x)
double pcdf(int mu, int mus);
double mulikely(int mu, int mus)
{
    //double cd=stats::ppois(mus,mu);
    double cd=pcdf(mus,mu);
    return min(cd,1.0-cd);
}
inline double poisson(int x,int lambda)
{
    if (lambda==0) return 0.0;
    //return gsl_ran_poisson_pdf(x,lambda);
    return stats::dpois(x,lambda);
    //return boost::math::pdf(x,lambda);
    //return pow(lambda,x)/tgamma(x+1)*exp(-lambda);
}
inline double pcdf(int mus, int mu)
{
    if (mu==0) return 0.0;
    //return gsl_cdf_poisson_P(mus,mu);
    return stats::ppois(mus,mu);
    //double Result=0.0;
    //for (int i=0;i<=mus;++i) Result+=poisson(i,mu);
    //return Result;
}

double getScore(Cand *TheCand, double **RDWsAcc, double* StandardsAcc,int SampleN,int WindowN,double* CNPriors,int CNPN)
{
    double Score=0,P=1,CN2L=1;
    for (int i=0;i<TheCand->Size;++i)
    {
        int mu=0;
        int mus=0;
        int SampleI=TheCand->SampleIs[i];
        mus=RDWsAcc[SampleI][TheCand->EEs[i]]-RDWsAcc[SampleI][TheCand->EBs[i]]+0.5;
        mu=StandardsAcc[TheCand->EEs[i]]-StandardsAcc[TheCand->EBs[i]]+0.5;
        
        double PSum=0.0;
        int SumN=0;
        for (int i=double(mu)*0.75;i<double(mu)*1.25;++i)
        {
            PSum+=mulikely(i,mus);
            ++SumN;
        }
        PSum/=SumN;
        CN2L*=PSum;
        TheCand->PassCs[i]=1.0-PSum;

        int eCN=TheCand->ECNs[i];
        double MP=0,MCN=0,Pmus=0;
        //CNPN=len(CNPriors)-1
        
        for (int j=0;j<CNPN;++j)
            Pmus+=CNPriors[j]*poisson(mus,mu*j/2==0?1:mu*j/2);
        if (Pmus==0)
        {
            MCN=eCN;
            MP=1;
        }
        else{
            for (int CN=min(0,eCN-1);CN<min(CNPN+1,eCN+2);++CN){
                if (CN>CNPN)
                    CN=CNPN;
                double Pmuscn=poisson(mus,mu*CN/2==0?1:mu*CN/2)*CNPriors[CN];
                double Pd=Pmuscn/Pmus;
                if (Pd>MP)
                {
                    MP=Pd;
                    MCN=CN;
                }
            }
        }
        TheCand->ECNs[i]=MCN;
        TheCand->Confs[i]=MP;
    }
    return 1-CN2L;
}

double* getScores(Cand *Cands,int Size, double **RDWsAcc, double* StandardsAcc,int SampleN,int WindowN,double* CNPriors,int CNPN,int ThreadN)
{
    double * Scores=(double*)malloc(sizeof(double)*Size);
    omp_set_num_threads(ThreadN);
    #pragma omp parallel for
    for (int i=0;i<Size;++i)
    {
        Scores[i]=getScore(Cands+i,RDWsAcc,StandardsAcc,SampleN,WindowN,CNPriors,CNPN);
    }
    return Scores;
}


PyObject* getRDScores(PyObject *self, PyObject *args)
{
    //Py_Initialize();
    fprintf(stderr,"getting vars...");
    fflush(stderr);
    PyObject *pModule;
    PyObject *Candidates, *TheContig, *PyThreadN;
    PyArg_ParseTuple(args,"OOO",&Candidates,&TheContig,&PyThreadN);
    int ThreadN=PyLong_AsLong(PyThreadN);

    
    PyObject * PyCNPriors= PyObject_GetAttrString(PyImport_Import(PyUnicode_FromString("consts")),"CNPriors");
    int CNPN=int(PyList_Size(PyCNPriors))-1;
    
    double *CNPriors=(double*)malloc(sizeof(double)*(CNPN+1));
    for (int i=0;i<=CNPN;++i) CNPriors[i]=PyFloat_AsDouble(PyList_GetItem(PyCNPriors,i));
    
    PyObject* Candidate;

    int i=0;
    int Size=int(PyList_Size(Candidates));

    PyObject * RDWindowsAcc=PyObject_GetAttrString(TheContig,"RDWindowsAcc");
    PyObject * RDStandardsAcc=PyObject_GetAttrString(TheContig,"RDWindowStandardsAcc");
    
    int WindowN=int(PyLong_AsLong(PyObject_CallMethod(RDStandardsAcc,"__len__",NULL)))-1;
    int SampleN=int(PyList_Size(RDWindowsAcc));

    double *StandardsAcc=(double*)malloc((WindowN+1)*sizeof(double));
    double* *WindowsAcc=(double * *)malloc(sizeof(double*)*SampleN);
    for (int j=0;j<SampleN;++j)
    {
        WindowsAcc[j]=(double*)malloc(sizeof(double)*(WindowN+1));
    }
    for (i=0;i<SampleN;++i)
    for (int j=0;j<=WindowN;++j)
    {
        WindowsAcc[i][j]=PyFloat_AsDouble(PyObject_CallMethod(PyList_GetItem(RDWindowsAcc,i),"__getitem__","O",PyLong_FromLong(j)));
    }
    
    for (i=0;i<WindowN;++i)
    {
        StandardsAcc[i]=PyFloat_AsDouble(PyObject_CallMethod(RDStandardsAcc,"__getitem__","O",PyLong_FromLong(i)));
    }

    Cand * Cands=(Cand*)malloc(sizeof(Cand)*Size);

    for (i=0;i<Size;++i)
    {
        Candidate=PyList_GetItem(Candidates,i);
        PyObject * Es=PyObject_GetAttrString(Candidate,"Evidences");
        Cands[i].Size=int(PyList_Size(Es));
        allocCand(Cands+i);
        for (int j=0;j<Cands[i].Size;++j)
        {
            PyObject * E=PyList_GetItem(Es,j);
            PyObject * EData=PyObject_GetAttrString(E,"Data");
            Cands[i].EBs[j]=PyLong_AsLong(PyObject_GetAttrString(EData,"WBegin"));//Es[j].Data.WBegin
            Cands[i].EEs[j]=PyLong_AsLong(PyObject_GetAttrString(EData,"WEnd"));
            Cands[i].SampleIs[j]=PyLong_AsLong(PyObject_GetAttrString(EData,"Sample"));
            Cands[i].ECNs[j]=PyLong_AsLong(PyObject_GetAttrString(EData,"CN"));
        }
    }

    fprintf(stderr,"calculating scores...");
    double * Scores=getScores(Cands, Size, WindowsAcc, StandardsAcc, SampleN, WindowN, CNPriors,CNPN,ThreadN);
    
    fprintf(stderr,"writing back...");

    PyObject *Results=PyList_New(0);
    
    for (i=0;i<Size;++i)
    {
        Candidate=PyList_GetItem(Candidates,i);
        PyObject * Es=PyObject_GetAttrString(Candidate,"Evidences");
        PyList_Append(Results,PyFloat_FromDouble(Scores[i]));
        for (int j=0;j<Cands[i].Size;++j)
        {
            PyObject_SetAttrString(PyObject_GetAttrString(PyList_GetItem(Es,j),"Data"),"CN",PyLong_FromLong(Cands[i].ECNs[j]));
            PyObject_SetAttrString(PyList_GetItem(Es,j),"PassConfidence",PyFloat_FromDouble(Cands[i].PassCs[j]));
            PyObject_SetAttrString(PyList_GetItem(Es,j),"Confidence",PyFloat_FromDouble(Cands[i].Confs[j]));
        }
    }

    fprintf(stderr,"cleaning up...");
    free(StandardsAcc);
    free(Scores);
    free(CNPriors);
    for (i=0;i<SampleN;++i) free(WindowsAcc[i]);
    for (i=0;i<Size;++i) freeCand(Cands+i);
    free(Cands);
    fprintf(stderr,"before return fine.Size:%d",i);
    return Results;
}

/*
PyObject* getRDScores(PyObject *self, PyObject *args)
{
    Py_Initialize();
    fprintf(stderr,"enter fine.");
    PyObject *pModule,*pGetRDScore;
    PyObject *pArgs, *pValue, *Candidates, *TheContig;
    
    PyArg_ParseTuple(args,"OO",&Candidates,&TheContig);

    // import 
    pModule = PyImport_Import(PyUnicode_FromString("calling"));

    // great_module.great_function 
    pGetRDScore = PyObject_GetAttrString(pModule, "getRDScore"); 
    //fprintf(stderr,"name:%s",_PyUnicode_AsString(PyObject_CallObject(PyObject_GetAttrString(pGetRDScore,"__str__"),NULL)));
    
    // build args
    pArgs = PyTuple_New(2);
    PyObject* Candidate;
    PyObject *Results=PyList_New(PyList_Size(Candidates));

    fprintf(stderr,"here fine.");
    //long Size=PyLong_AsLong(PyList_Size(Candidates))
    int i=0;
    int Size=int(PyList_Size(Candidates));

    //omp_set_num_threads(8);
    //#pragma omp parallel for private(Candidate,i,pValue)
    for (i=0;i<Size;++i)
    {
        Candidate=PyList_GetItem(Candidates,i);
    //fprintf(stderr,"name:%s",_PyUnicode_AsString(PyObject_CallObject(PyObject_GetAttrString(Candidate,"__str__"),NULL)));
        //PyTuple_SetItem(pArgs,0, Candidate);
        //PyTuple_SetItem(pArgs,1, TheContig);
        // call 
        //pValue = PyObject_CallObject(pGetRDScore, pArgs);
        //pValue=PyObject_CallFunction(pGetRDScore,"OO",Candidate,TheContig);
        pValue=PyObject_CallMethod(pModule,"getRDScore","OO",Candidate,TheContig);
        PyList_SetItem(Results,i,pValue);
    }
    
        fprintf(stderr,"before return fine.Size:%d",i);
    return Results;
}*/

static PyMethodDef GRDMethods[] = {
    {"CGetRDScores",  getRDScores, METH_VARARGS,
     "get RD scores"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef CGetRDScores =
{
    PyModuleDef_HEAD_INIT,
    "CGetRDScores", /* name of module */
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    GRDMethods
};

PyMODINIT_FUNC PyInit_CGetRDScores(void)
{
    return PyModule_Create(&CGetRDScores);
}