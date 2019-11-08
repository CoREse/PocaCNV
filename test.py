import subprocess
import time
import os
import glob
import sys

Note="""
change reference
"""
Reference="/home/cre/data/GRCh38_full_analysis_set_plus_decoy_hla.fa"
Samples="/home/cre/data/HG0051*chr22.cram"
PythonName="python3"
ProgramName="jc.py"
TestDirRoot="test"


RunDir=TestDirRoot+"/"+"TestRun"+time.strftime("%Y%m%d%H%M%S",time.localtime())
SamplesArray=glob.glob(Samples)

run=subprocess.Popen(["ls",TestDirRoot],bufsize=1,stdout=subprocess.PIPE,universal_newlines=True)
Runs=[]
while True:
    line=run.stdout.readline()
    if line=="":
        break
    if len(line)<7:
        continue
    if line[:7]=="TestRun":
        Runs.append(line.strip())
run.wait()
Runs.sort()
LastDir=TestDirRoot+"/"+Runs[-1]
LastSource=LastDir+"/src"

while os.path.isdir(RunDir):
    RunDir+="+"
os.mkdir(RunDir)
SourceDir=RunDir+"/src"
os.mkdir(SourceDir)

OutFile=open(RunDir+"/"+"out.txt","w")
ErrFile=open(RunDir+"/"+"err.txt","w")
LogFile=open(RunDir+"/"+"log.txt","w")
PatchFile=open(RunDir+"/"+"diff.patch","w")

subprocess.Popen(["cp * %s"%(SourceDir)],shell=True).wait()
subprocess.Popen(["diff","-N",SourceDir,LastSource],bufsize=1,stdout=PatchFile,stderr=subprocess.PIPE,universal_newlines=True).wait()
PatchFile.close()

print("Run %s\nStart at:%s\nProgram:%s\nReference:%s\nSamples:%s\nNote:%s\n"%(RunDir.split("/")[-1],ProgramName,time.ctime(),Reference,Samples,Note),file=LogFile)
LogFile.flush()

Start=time.time()
run=subprocess.Popen([PythonName,"-u",ProgramName,Reference]+SamplesArray,bufsize=1,stdout=OutFile,stderr=subprocess.PIPE,universal_newlines=True)
while True:
    line=run.stderr.readline()
    if line=="":
        break
    ErrFile.write(line)
    print(line, end="", file=sys.stderr)
run.wait()
print("Time elapsed:%s"%(time.time()-Start),file=LogFile)
subprocess.Popen([PythonName,"benchmark.py",RunDir+"/"+"out.txt"],stdout=LogFile).wait()

OutFile.close()
ErrFile.close()
LogFile.close()

subprocess.call(["cat", RunDir+"/"+"log.txt"])