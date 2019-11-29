import sys
mul=int(sys.argv[1])
pref=sys.argv[2]
for fn in sys.argv[3:]:
    SampleName=fn.split("\\")[-1].split("/")[-1][:-4]
    print("transforming %s..."%fn,file=sys.stderr)
    nfn="%s/%s%s.rdf"%(pref,mul,SampleName)
    nf=open(nfn,"w")
    first=True
    with open(fn,"r") as ifile:
        ConI=0
        data=0
        for line in ifile:
            if line[0]=="#":
                if ConI%mul!=0:
                    print("\n%s"%(data),end="",file=nf)
                chrname=line[1:].split()[0]
                length=int(line[1:].split()[1])
                newl=int(length/mul)+(1 if length%mul!=0 else 0)
                if not first:
                    print("\n",end="",file=nf)
                first=False
                print("#%s %s"%(chrname,newl),end="",file=nf)
                data=0
                ConI=0
                continue
            data+=int(line)
            ConI+=1
            if ConI%mul==0:
                print("\n%s"%(data),end="",file=nf)
                data=0
        if ConI%mul!=0:
            print("\n%s"%(data),end="",file=nf)
