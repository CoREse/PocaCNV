import subprocess
import os
import sys
from utils import gettime
from multiprocessing import Pool,Manager

Samples=[]
i=1
ThreadN=1
while i < len(sys.argv):
    if sys.argv[i]=="-j" or sys.argv[i]=="-J":
        ThreadN=int(sys.argv[i+1])
        i+=1
    else:
        Samples.append(sys.argv[i])
    i+=1
Reference=Samples[0]
Samples=Samples[1:]
print(Reference,Samples)
SamplesArray=Samples
SamplesArray.sort()

PythonName="python3"
ProgramName="jcrd.py"
print(gettime()+"Starting making RD data...",file=sys.stderr)
def run1(Ref,Samp):
    print(gettime()+"Making RD data for %s..."%(Samp),file=sys.stderr)    
    run=subprocess.Popen([PythonName,"-u",ProgramName,'-T',Reference,Samp,"-WO"],universal_newlines=True)
    run.wait()
    print(gettime()+"Made RD data for %s."%(Samp),file=sys.stderr) 

if ThreadN!=1:
    MyPool=Pool(ThreadN)
    MyManager=Manager()
    Args=[]
    for S in Samples:
        Args.append((Reference,S))
    MyPool.starmap(run1,Args)
else:
    for S in Samples:
        run1(Reference,S)
print(gettime()+"Making RD data done.",file=sys.stderr)