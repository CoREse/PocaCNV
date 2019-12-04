import rpy2.robjects as robjects
datastring="1,2,3"
robjects.r("rddata=data.frame(mrd=c(%s))"%datastring)
sf=open("dnacopy_cbs.r","r")
script=str(sf.read())
#print(script)
sf.close()
x=robjects.r(script)
x:robjects.ListVector
b=x[1]
b: robjects.DataFrame
print(b.names)
v=b[3]
v:robjects.vectors.FloatVector
print(len(b[3]))
for d in b[3]:
    print(d)