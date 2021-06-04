import sys
ISAV=None
ISSD=None
for line in sys.stdin:
    sl=line.split("\t")
    if len(sl)<3:
        continue
    if ISAV!=None and ISSD!=None:
        print(ISAV,ISSD,end="")
        exit(0)
    if sl[1]=="insert size average:":
        ISAV=sl[2].strip()
        continue
    if sl[1]=="insert size standard deviation:":
        ISSD=sl[2].strip()
        continue