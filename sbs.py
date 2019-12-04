from scipy.stats import norm
import sys
from utils import gettime
vconfidence=0.999
steplength=10
verbose=False
sums=None

def bsegment(data,begin,end):
    if end<=begin:
        return begin
    tv=getv(data,begin,end)
    sep=begin
    for i in range(begin+4,end-3):#reduce coincidences
        lv=getv(data,begin,i)
        rv=getv(data,i,end)
        if tv>lv+rv:
            tv=lv+rv
            sep=i
    if sep==begin:
        return begin
    if validate(data,begin,end,sep):
        return sep
    return begin

def getv(data,begin,end):
    if begin==end:
        return 0
    v=0
    a=0
    if sums==None:
        for i in range(begin,end):
            a+=data[i]
    else:
        a=sums[end]-sums[begin]
    a/=end-begin
    for i in range(begin,end):
        v+=(data[i]-a)**2
    v/=end-begin
    return v

def getav(data,begin,end):
    if begin==end:
        return 0,0
    v=0
    a=0
    if sums==None:
        for i in range(begin,end):
            a+=data[i]
    else:
        a=sums[end]-sums[begin]
    a/=end-begin
    for i in range(begin,end):
        v+=(data[i]-a)**2
    v/=end-begin
    return a,v

def geta(data,begin,end):
    if begin==end:
        return 0,0
    a=0
    if sums==None:
        for i in range(begin,end):
            a+=data[i]
    else:
        a=sums[end]-sums[begin]
    a/=end-begin
    return a

def validate(data,begin,end,i):
    la,lv=getav(data,begin,i)
    ra,rv=getav(data,i,end)
    if lv==0:
        lv=0.1
    if rv==0:
        rv=0.1
    lint=norm.interval(vconfidence,la,lv/(end-i))
    rint=norm.interval(vconfidence,ra,rv/(i-begin))
    if lint[0]<ra<lint[1] or rint[0]<la<rint[1]:
        return False
    return True

def segment(data,begin=None,end=None):
    global sums
    print(gettime()+"calculating sums...",file=sys.stderr)
    if sums==None:
        sums=[]
        s=0
        for d in data:
            sums.append(s)
            s+=d
        sums.append(s)
    if begin==None:
        begin=0
    if end==None:
        end=len(data)
    i=begin+steplength
    separators=[]
    ls=begin
    countstep=100000
    count=countstep
    print(gettime()+"step segmenting...",file=sys.stderr)
    while i<end:
        if verbose:
            if i>count:
                print("%s/%s analyzed."%(count,end),file=sys.stderr)
                count+=countstep
        if validate(data,ls,(i+steplength) if i+steplength<end else end,i):
            separators.append(i)
            ls=i
            i+=steplength
        else:
            i+=1
    separators.append(end)
    '''
    print(separators)
    print(gettime()+"2 bsegmenting...",file=sys.stderr)
    nseps=[]
    ls=begin
    for s in separators:
        sep=bsegment(data,ls,s)
        if sep!=ls:
            nseps.append(sep)
        ls=s
        nseps.append(s)
    ls=begin
    nnseps=[]
    print(gettime()+"validating...",file=sys.stderr)
    for i in range(len(nseps)-1):
        s=nseps[i]
        e=nseps[i+1]
        if validate(data,ls,e,s):
            nnseps.append(s)
        ls=s'''
    segments=[]
    ls=begin
    nnseps=separators
    for s in nnseps:
        segments.append((s,geta(data,ls,s)))
        ls=s
    return segments

if __name__ == '__main__':
    dfile=open(sys.argv[1],"r")
    sample=[]
    for line in dfile:
        if line[0]=="#":
            continue
        sample.append(float(line.split()[-1]))
    dfile.close()
    #sample = generate_normal_time_series(5)
    #print(sample)
    print(gettime()+"segmenting...",file=sys.stderr)
    verbose=True
    L = segment(sample)
    first=True
    for s in L:
        if first:
            print("\n",end="")
        print("%s %.6s"%(s[0],s[1]),end="")