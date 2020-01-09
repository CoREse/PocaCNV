#include <Python.h>
#include <stdio.h>

PyObject* getRDScores(PyObject *self, PyObject *args)
{
    fprintf(stderr,"enter fine.");
    PyObject *pModule,*pGetRDScore;
    PyObject *pArgs, *pValue, *Candidates, *TheContig;
    
    PyArg_ParseTuple(args,"OO",&Candidates,&TheContig);

    /* import */
    pModule = PyImport_Import(PyUnicode_FromString("calling"));

    /* great_module.great_function */
    pGetRDScore = PyObject_GetAttrString(pModule, "getRDScore"); 
    //fprintf(stderr,"name:%s",_PyUnicode_AsString(PyObject_CallObject(PyObject_GetAttrString(pGetRDScore,"__str__"),NULL)));
    
    /* build args */
    pArgs = PyTuple_New(2);
    PyObject* Candidate;
    PyObject *Results=PyList_New(PyList_Size(Candidates));

    fprintf(stderr,"here fine.");
    //long Size=PyLong_AsLong(PyList_Size(Candidates))
    int i=0;
    for (i=0;i<PyList_Size(Candidates);++i)
    {
        Candidate=PyList_GetItem(Candidates,i);
    //fprintf(stderr,"name:%s",_PyUnicode_AsString(PyObject_CallObject(PyObject_GetAttrString(Candidate,"__str__"),NULL)));
        //PyTuple_SetItem(pArgs,0, Candidate);
        //PyTuple_SetItem(pArgs,1, TheContig);
        /* call */
        //pValue = PyObject_CallObject(pGetRDScore, pArgs);
        //pValue=PyObject_CallFunction(pGetRDScore,"OO",Candidate,TheContig);
        pValue=PyObject_CallMethod(pModule,"getRDScore","OO",Candidate,TheContig);
        PyList_SetItem(Results,i,pValue);
    }
    
        fprintf(stderr,"before return fine.Size:%d",i);
    return Results;
}

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