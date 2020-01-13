#include <Python.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "include/stats.hpp"

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
double mulikely(int mu, int mus)
{
    double cd=stats::ppois(mus,mu);
    return min(cd,1.0-cd);
}

double getScore(Cand *TheCand, double **RDWsAcc, double* StandardsAcc,int SampleN,int WindowN,double* CNPriors,int CNPN)
{
    double Score=0,P=1,CN2L=1;
    for (int i=0;i<TheCand->Size;++i)
    {
        int mu=0;
        int mus=0;
        int SampleI=TheCand->SampleIs[i];
        fprintf(stderr,"fine2.");
        mus=RDWsAcc[SampleI][TheCand->EEs[i]]-RDWsAcc[SampleI][TheCand->EBs[i]]+0.5;
        mu=StandardsAcc[TheCand->EEs[i]]-StandardsAcc[TheCand->EBs[i]]+0.5;
        fprintf(stderr,"fine3.");
        
        TheCand->PassCs[i]=1.0-mulikely(mu,mus);
        CN2L*=1.0-TheCand->PassCs[i];

        fprintf(stderr,"fine4.");
        int eCN=TheCand->ECNs[i];
        double MP=0,MCN=0,Pmus=0;
        //CNPN=len(CNPriors)-1
        
        for (int j=0;j<CNPN;++j)
            Pmus+=CNPriors[j]*stats::dpois(mus,int(mu*j/2));
        if (Pmus==0)
        {
            MCN=eCN;
            MP=1;
        }
        else{
            for (int CN=min(0,eCN-1);CN<eCN+2;++CN){
                if (CN>CNPN)
                    CN=CNPN;
                double Pmuscn=stats::dpois(mus,int(mu*CN/2))*CNPriors[CN];
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

double* getScores(Cand *Cands,int Size, double **RDWsAcc, double* StandardsAcc,int SampleN,int WindowN,double* CNPriors,int CNPN)
{
    double * Scores=(double*)malloc(sizeof(double)*Size);
    for (int i=0;i<Size;++i)
    {
        fprintf(stderr,"fine1.--%d--",i);
        Scores[i]=getScore(Cands+i,RDWsAcc,StandardsAcc,SampleN,WindowN,CNPriors,CNPN);
    }
    return Scores;
}


PyObject* getRDScores(PyObject *self, PyObject *args)
{
    Py_Initialize();
    fprintf(stderr,"getting vars...");
    PyObject *pModule,*pGetRDScore;
    PyObject *pArgs, *pValue, *Candidates, *TheContig;
    
    PyArg_ParseTuple(args,"OO",&Candidates,&TheContig);

    // import
    pModule = PyImport_Import(PyUnicode_FromString("calling"));
    PyObject * PyCNPriors= PyObject_GetAttrString(PyImport_Import(PyUnicode_FromString("consts")),"CNPriors");
    int CNPN=int(PyList_Size(PyCNPriors))-1;
    double *CNPriors=(double*)malloc(sizeof(double)*(CNPN+1));
    for (int i=0;i<=CNPN;++i) CNPriors[i]=PyFloat_AsDouble(PyList_GetItem(PyCNPriors,i));

    // great_module.great_function
    pGetRDScore = PyObject_GetAttrString(pModule, "getRDScore"); 
    //fprintf(stderr,"name:%s",_PyUnicode_AsString(PyObject_CallObject(PyObject_GetAttrString(pGetRDScore,"__str__"),NULL)));
    
    // build args 
    pArgs = PyTuple_New(2);
    PyObject* Candidate;

    //long Size=PyLong_AsLong(PyList_Size(Candidates))
    int i=0;
    int Size=int(PyList_Size(Candidates));
    //double *Scores=(double *)malloc(sizeof(double)*Size);

    PyObject * RDWindowsAcc=PyObject_GetAttrString(TheContig,"RDWindowsAcc");
    PyObject * RDStandardsAcc=PyObject_GetAttrString(TheContig,"RDWindowStandardsAcc");
    
    int WindowN=int(PyList_Size(RDStandardsAcc));
    int SampleN=int(PyList_Size(RDWindowsAcc));

    double *StandardsAcc=(double*)malloc((WindowN+1)*sizeof(double));
    double* *WindowsAcc=(double * *)malloc(sizeof(double*)*SampleN);
    for (int j=0;j<SampleN;++j)
    {
        WindowsAcc[j]=(double*)malloc(sizeof(double)*(WindowN+1));
    }
    for (i=0;i<SampleN;++i)
    for (int j=0;j<WindowN;++j)
    {
        WindowsAcc[i][j]=PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(RDWindowsAcc,i),j));
    }
    for (i=0;i<WindowN;++i)
    {
        StandardsAcc[i]=PyFloat_AsDouble(PyList_GetItem(RDStandardsAcc,i));
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
    double * Scores=getScores(Cands, Size, WindowsAcc, StandardsAcc, SampleN, WindowN, CNPriors,CNPN);
    
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