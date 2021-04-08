#!python
#cython: language_level=3, std=c++11
cimport cython
from cython.parallel import parallel,prange
import ctypes
from libc.stdlib cimport malloc,free
from libc.math cimport pow,fabs
from cpython cimport array
from array import array
import sys
from utils import *

import globals as g

def getSampleSum(TheContig, SampleI, WBegin, WEnd):
    SampleRD=0
    for j in range(WBegin,WEnd):
        SampleRD+=TheContig.RDWindows[SampleI][j]
    return SampleRD

cdef double CgetSampleSum(float ** CRDWindowsAcc, unsigned long SampleI, unsigned long WBegin, unsigned long WEnd) nogil:
    #cdef double SampleRD=0
    #cdef unsigned long j=0
    #SampleRD=CRDWindowsAcc[SampleI][WEnd]-CRDWindowsAcc[SampleI][WBegin]
    #for j from WBegin<=j<WEnd:
    #    SampleRD+=CRDWindows[SampleI][j]
    return CRDWindowsAcc[SampleI][WEnd]-CRDWindowsAcc[SampleI][WBegin]

from libcpp.unordered_set cimport unordered_set
cdef void CgetSP(double *SP, float ** CRDwindowsAcc, unsigned long SampleN, double * SampleReadCount, unsigned long WBegin, unsigned long WEnd, double NSD=3, double MinimumTake=0.8) nogil:
    cdef double SRS=0,SRC=0
    if WEnd<=WBegin:
        SP[0]=0
        SP[1]=0
    cdef double *SampleRDs=<double*>malloc(SampleN*sizeof(double))
    cdef unsigned long i,j
    for i from 0<=i<SampleN:
        SampleRDs[i]=CgetSampleSum(CRDwindowsAcc, i, WBegin, WEnd)
        SRS+=SampleRDs[i]
        SRC+=SampleReadCount[i]
    cdef double EstimatedP=SRS/SRC
    cdef unordered_set[unsigned long] RemovedSet
    cdef double LastP=0
    cdef double *SampleSTDs=<double*>malloc(SampleN*sizeof(double))
    while LastP!=EstimatedP:
        RemovedSet.clear()
        for i from 0<=i<SampleN:
            SampleSTDs[i]=EstimatedP*SampleReadCount[i]
            if fabs(SampleRDs[i]-SampleSTDs[i])>NSD*pow(SampleSTDs[i],0.5):
                RemovedSet.insert(i)
        if RemovedSet.size()>SampleN*MinimumTake:
            break
        SRS=0
        SRC=0
        for i from 0<=i<SampleN:
            if RemovedSet.count(i)==0:
                SRS+=SampleRDs[i]
                SRC+=SampleReadCount[i]
        LastP=EstimatedP
        EstimatedP=SRS/SRC
    SP[0]=SRS
    SP[1]=SRC
    free(SampleRDs)
    free(SampleSTDs)
    return

#should be irrelevant to ploidy if we use local, since all windows in the same sample have the same ploidies.
def getSP(TheContig, WBegin, WEnd, NSD=3, MinimumTake=0.8, local=g.StatLocal):#get rd sum and sample read count sum
    SRS=0
    SRC=0
    if WEnd<=WBegin:
        return (0,0)
    SampleN=len(TheContig.SampleNames)
    SampleRDs=[0]*SampleN
    Stat=g
    if local:
        Stat=TheContig
    #start with middle p
    #this seems better
    '''SamplePs=[0]*SampleN
    for i in range(SampleN):
        SampleRDs[i]=getSampleSum(TheContig,i,WBegin,WEnd)
        SamplePs[i]=SampleRDs[i]/Stat.SampleReadCount[i]
    SamplePs.sort()
    EstimatedP=SamplePs[int(SampleN/2)]
    SRS=EstimatedP
    SRC=1
    if SampleN%2==0:
        EstimatedP+=SamplePs[int(SampleN/2)-1]
        EstimatedP/=2'''
    #start with whole
    for i in range(SampleN):
        SampleRDs[i]=getSampleSum(TheContig,i,WBegin,WEnd)
        SRS+=SampleRDs[i]
        SRC+=Stat.SampleReadCount[i]
    EstimatedP=SRS/SRC
    RemovedSet=set()
    LastP=0
    SampleSTDs=[0]*SampleN
    while LastP!=EstimatedP:
        RemovedSet=set()
        for i in range(SampleN):
            SampleSTDs[i]=EstimatedP*Stat.SampleReadCount[i]
            if abs(SampleRDs[i]-SampleSTDs[i])>NSD*SampleSTDs[i]**0.5:
                RemovedSet.add(i)
        if len(RemovedSet)>SampleN*MinimumTake:
            break
        SRS=0
        SRC=0
        for i in range(SampleN):
            if i not in RemovedSet:
                SRS+=SampleRDs[i]
                SRC+=Stat.SampleReadCount[i]
        LastP=EstimatedP
        EstimatedP=SRS/SRC
    #the efficient way
    '''for i in range(WBegin,WEnd):
        SRS+=g.StatisticalRDWindowSums[i]
        SRC+=g.StatisticalReadCounts[i]
    SRC/=WEnd-WBegin'''
    #the efficient way 2
    '''for i in range(WBegin,WEnd):
        SRS+=g.RDWindowSums[i]
    SRC=g.AllReadCount'''
    return (SRS,SRC)

cdef double CgetNormalRDDirect(double * StatisticalReadCounts, double * StatisticalRDWindowSums, double * SampleReadCount, unsigned long SampleI, unsigned long WindowI) nogil:
    return 0 if StatisticalReadCounts[WindowI]==0 else StatisticalRDWindowSums[WindowI]*(SampleReadCount[SampleI]/StatisticalReadCounts[WindowI])

def processingRD(RDWindows, SampleN, WindowsN, MixedRDRs, RDWindowSums, RDWindowsAcc, StatisticalReadCounts, StatisticalRDWindowSums, SampleReadCount, Ploidies):
    Mean2Adjust=False
    Smooth=1
    RDWindowAverages=array("d",[0]*WindowsN)

    
    cdef double * SampleSums=<double *>malloc(SampleN*sizeof(double))
    cdef double * SampleAverages=<double *>malloc(SampleN*sizeof(double))
    cdef double * CPloidies=<double *>malloc(SampleN*sizeof(double))
    for i in range(SampleN):
        SampleSums[i]=0
        SampleAverages[i]=0
        CPloidies[i]=Ploidies[i]
    SampleReadCountsArray=array("d",SampleReadCount)
    #cdef size_t SampleReadCountsArrayAddress=ctypes.addressof(SampleReadCountsArray)
    #cdef unsigned long CSampleReadCount=TheContig.SampleReadCount
    cdef unsigned long CSampleN=SampleN
    cdef unsigned long CWindowsN=WindowsN
    #cdef size_t CSRCA=ctypes.addressof(StatisticalReadCounts)
    #cdef size_t CSRDWA=ctypes.addressof(StatisticalRDWindowSums)
    cdef array.array TempArray
    TempArray=SampleReadCountsArray
    cdef double* CSampleReadCount=<double*>TempArray.data.as_voidptr
    TempArray=StatisticalReadCounts
    cdef double* CStatisticalReadCounts=<double*>TempArray.data.as_voidptr
    TempArray=StatisticalRDWindowSums
    cdef double* CStatisticalRDWindowSums=<double*>TempArray.data.as_voidptr
    TempArray=RDWindowSums
    cdef double* CRDWindowSums=<double*>TempArray.data.as_voidptr
    TempArray=RDWindowAverages
    cdef double* CRDWindowAverages=<double*>TempArray.data.as_voidptr
    cdef double Standard
    cdef double WR
    cdef unsigned long Cj=0,Ci=0
    cdef float ** CRDWindows=<float **>malloc(SampleN*sizeof(float*))
    cdef float ** CMixedRDRs=<float **>malloc(SampleN*sizeof(float*))
    cdef float ** CRDWindowsAcc=<float **>malloc(SampleN*sizeof(float*))
    #cdef size_t CTempAddress=ctypes.addressof(RDWindowAverages)
    cdef size_t TempAddr
    cdef int CThreadN=g.ThreadN
    #RDWindowsArrays=[]
    #for i in range(SampleN):
    #    RDWindowsArrays.append(array("d",RDWindows[i]))

    for i in range(SampleN):
        #CTempAddress=ctypes.addressof(RDWinodws[i])
        #TempArray=RDWindows[i]
        #TempAddr=ctypes.addressof(RDWindows[i])
        #TempAddr=TempArray.data.as_voidptr
        #TempArray=RDWindowsArrays[i]
        TempArray=RDWindows[i]
        #CRDWindows[i]=<double*>(TempAddr)
        CRDWindows[i]=<float*>(TempArray.data.as_voidptr)
        #CTempAddress=ctypes.addressof(MixedRDRs[i])
        TempArray=MixedRDRs[i]
        CMixedRDRs[i]=<float*>(TempArray.data.as_voidptr)
        TempArray=RDWindowsAcc[i]
        CRDWindowsAcc[i]=<float*>(TempArray.data.as_voidptr)
    TempArray=None

    cdef float Sum
    for Ci in prange(CSampleN,nogil=True,num_threads=CThreadN):
        Sum=0
        for Cj from 0<=Cj<CWindowsN:
            CRDWindowsAcc[Ci][Cj]=Sum
            Sum=CRDWindows[Ci][Cj]+Sum
        CRDWindowsAcc[Ci][CWindowsN]=Sum
    
    cdef double * SP
    '''
    spf=open("data/cSPs.txt","w")
    mrdf=open("data/c1stMRDs.txt","w")
    stf=open("data/cStandards.txt","w")
    avef=open("data/cAverages.txt","w")
    wrf=open("data/cWRs.txt","w")
    wr2f=open("data/cWR2s.txt","w")'''
    #calc MRDs
    with nogil, parallel(num_threads=CThreadN):
        SP=<double *>malloc(2*sizeof(double))
        for Ci in prange(CWindowsN, schedule="guided"):
            for Cj from 0<=Cj<CSampleN:
                CRDWindowSums[Ci]+=CRDWindows[Cj][Ci]
                SampleSums[Cj]+=CRDWindows[Cj][Ci]
            #for Cj from 0<=Cj<100000000:
            #    CRDWindows[0][0]+=0.1*pow(Cj/100000000,0.5)
            CRDWindowAverages[Ci]=CRDWindowSums[Ci]/<double>CSampleN
            CgetSP(SP, CRDWindowsAcc, CSampleN, CSampleReadCount, Ci, Ci+1)
            #oSP=getSP(TheContig,Ci,Ci+1)
            #print("(%s, %d)"%(SP[0],int(SP[1])),file=spf)
            CStatisticalRDWindowSums[Ci]=SP[0]
            CStatisticalReadCounts[Ci]=SP[1]
        free(SP)
    
    AllZeroLeft=-1
    AllZeroRight=WindowsN
    Left1=False
    Right1=False
    for i in range(WindowsN):
        if not Right1:
            j=WindowsN-i
            if j>=0:
                if RDWindowSums[i]!=0:
                    Right1=True
                else:
                    AllZeroRight=j
        if not Left1:
            if RDWindowSums[i]!=0:
                Left1=True
            else:
                AllZeroLeft=i
        if Right1 and Left1:
            break
    
    cdef unsigned long S0Value
    cdef long CAllZeroLeft=AllZeroLeft
    cdef unsigned long CAllZeroRight=AllZeroRight

    #SampleSumAverage=0
    #for j in range(SampleN):
    #    SampleAverages[j]=SampleSums[j]/WindowsN
    #    SampleSumAverage+=SampleSums[j]
    #SampleSumAverage/=SampleN

    '''srdf=open("data/cSRDWindowsSums.txt","w")
    srcf=open("data/cSRCounts.txt","w")
    samrcf=open("data/cSampleReadCounts.txt","w")
    rdwf=open("data/c1stRDWins.txt","w")
    for i in range(WindowsN):
        print(CStatisticalRDWindowSums[i],file=srdf)
        print(int(CStatisticalReadCounts[i]),file=srcf)
        print(CRDWindows[0][i],file=rdwf)
        print(CRDWindowAverages[i],file=avef)
    for i in range(SampleN):
        print(int(SampleReadCount[i]),file=samrcf)'''
    print("Cython for start.",file=sys.stderr)
    for Ci in prange(CSampleN,nogil=True, num_threads=CThreadN):
        for Cj from 0<=Cj<CWindowsN:
        #for j in prange(WindowsN):
            S0Value=0
            if Cj<=CAllZeroLeft or Cj>=CAllZeroRight:
                S0Value=1
            Standard=CgetNormalRDDirect(CStatisticalReadCounts,CStatisticalRDWindowSums,CSampleReadCount,Ci,Cj)
            #Standard2=getNormalRD(TheContig,Ci,Cj)
            #print(Standard,file=stf)
            #if Standard!=Standard2:
            #    print("not identical!",Ci,Cj,g.StatisticalReadCounts[Cj],g.StatisticalRDWindowSums[Cj],TheContig.SampleReadCount[Ci],Standard2,CStatisticalReadCounts[Cj],CStatisticalRDWindowSums[Cj],CSampleReadCount[Ci],Standard, WR, WR2,file=sys.stderr)
            #    exit(0)
            WR=CRDWindows[Ci][Cj]/Standard if Standard!=0 else S0Value
            #WR2=RDWindows[Ci][Cj]/Standard2 if Standard2!=0 else S0Value
            #if WR!=WR2:
            #    print("WR not identical!",Ci,Cj,g.StatisticalReadCounts[Cj],g.StatisticalRDWindowSums[Cj],TheContig.SampleReadCount[Ci],Standard2,CStatisticalReadCounts[Cj],CStatisticalRDWindowSums[Cj],CSampleReadCount[Ci],Standard,file=sys.stderr)
            #    exit(0)
            #print(WR,file=wrf,end=" ")
            if Standard==0 and CRDWindows[Ci][Cj]!=0:
                WR=CRDWindows[Ci][Cj]/CRDWindowAverages[Cj]
            #print(WR,file=wrf)
            #print(WR*2,file=wr2f)
            CMixedRDRs[Ci][Cj]=WR
            CMixedRDRs[Ci][Cj]*=CPloidies[Ci]
    #for i in range(WindowsN):
    #    print(CMixedRDRs[0][i],file=mrdf)
    free(CRDWindows)
    free(CMixedRDRs)
    #for i in range(SampleN):
    #    for j in range(WindowsN):
    #        print(MixedRDRs[i][j])
    print("Cython test completed.",file=sys.stderr)

    if Smooth>1:
        if Smooth%2==0:
            Smooth+=1
        for i in range(SampleN):
            Smoother=array("f",[0]*WindowsN)
            for j in range(WindowsN):
                SH=int(Smooth/2)
                SB=0 if j<SH else j-SH
                SE=WindowsN if j+SH>WindowsN else j+SH
                for k in range(SB,SE):
                    Smoother[j]+=MixedRDRs[i][k]
                Smoother[j]/=Smooth
            MixedRDRs[i]=Smoother
    if Mean2Adjust:
        #sf=open("data/SampleCountMeanData.txt","w")
        #first=True
        for i in range(SampleN):
            #if not first:
            #    print("\n",end="",file=sf)
            #first=False
            Mean=0
            for j in range(WindowsN):
                Mean+=MixedRDRs[i][j]
            Mean/=WindowsN
            ARatio=2.0/Mean
            for j in range(WindowsN):
                MixedRDRs[i][j]*=ARatio

    free(SampleSums)
    free(SampleAverages)
    free(CPloidies)