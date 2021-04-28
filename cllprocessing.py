from utils import *

def processEvidencesWithCLL(Evidences,TheContig):
    for i in range(len(Evidences)):
        for j in range(len(Evidences[i])):
            ACL=0
            WBegin=int(Evidences[i][j].Begin/TheContig.RDWindowSize)
            WEnd=int(Evidences[i][j].End/TheContig.RDWindowSize)
            RDSum=TheContig.RDWindowsAcc[i][WEnd]-TheContig.RDWindowsAcc[i][WBegin]
            if RDSum==0:
                Evidences[i][j].ACL=ACL
                continue
            for k in range(WBegin,WEnd):
                ACL+=TheContig.AverageClipLengths[i][k]*TheContig.RDWindows[i][k]
            ACL/=RDSum
            Evidences[i][j].ACL=ACL
    return Evidences