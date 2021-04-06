import globals as g
from utils import *

#def indexDPRs(TheContig):
#    TheContig.SortedDRPsBySampleI={}
#    for i in range(len(TheContig.RDWindows)):
#        TheContig.SortedDRPsBySampleI[i]=[]
#    for d in TheContig.DRPs:
#        TheContig.SortedDRPsBySampleI[g.SampleNameIndexes[d.SampleName]].append(d)
#    for i in range(len(TheContig.RDWindows)):
#        TheContig.SortedDRPsBySampleI[i].sort(key=lambda d: d.Start)

def processEvidenceWithDRPs(e, TheContig):
    for d in TheContig.DRPs[e.Sample]:
        if d.Start>e.Begin:
            break
        if (d.SupportedVariantType ==0 and e.SupportedSVType> 0) or (d.SupportedVariantType==1 and e.SupportedSVType==0):
            continue 
        if inclusion((e.Begin,e.End),(d.Start,d.End))>2:#d include e or identical
            e.SupportedDRPs.append(d)

def processEvidencesWithDRPs(Es, TheContig):
    for e in Es:
        processEvidenceWithDRPs(e,TheContig)
    return Es