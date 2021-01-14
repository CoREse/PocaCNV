import torch
import torch.nn as nn
import sys

class TwoLayerNet(torch.nn.Module):
    #  定义了init和forward就想到于给定了一个模型
    #  init相当于给了这个框架
    def __init__(self, D_in, H, D_out):
        super(TwoLayerNet, self).__init__()
        self.linear1 = nn.Linear(D_in, H)
        self.linear2 = nn.Linear(H, D_out)

    #  forward是
    def forward(self, x):
        y_pred = self.linear2(torch.sigmoid(self.linear1(x)))
        return y_pred

SavePath="data/ScoringTrainModelData"
#SavePath="ScoringTrain/Model0192/Model0192.pickle"
#SavePath="ScoringTrain/Model0179/Model0179.pickle"

DataFile=open(sys.argv[2],"r")
Data=[]
Labels=[]
D0=0
D1=0
D=0

CNPriors=[1.4818669642343562e-06, 0.0024241621501567266, 0.9914186177260093, 0.00550095102642564, 0.0006000125666241987, 4.0389660750366106e-05, 1.1302181811873978e-05, 2.3256457943417974e-06, 5.854614618138028e-07, 7.124415377290356e-08, 1.0010930675953374e-07, 3.210044872566197e-10, 3.565148390433508e-11, 2.993635914964604e-12, 7.20618688429682e-13, 1.3465236956844422e-13, 3.0109007332307176e-14, 3.4497272100312798e-15, 2.5051103812841497e-15]

#ChromNo=int(sys.argv[1].split("_")[1])
SampleCount=int(sys.argv[1])
for line in DataFile:
    sl=line.split()
    SegNum=int(sl[0])
    ContigLength=int(sl[1])
    SiblingCount=int(sl[2])
    Start=int(sl[3])
    End=int(sl[4])
    Mu=int(sl[5])
    MuS=int(sl[6])
    PassConfidence=float(sl[7])
    CN=float(sl[8])
    Confidence=float(sl[9])
    CScore=float(sl[10])
    Label=int(sl[11])
    #DataTensor=torch.Tensor([SegNum,ContigLength,SiblingCount,Start,End,Mu,MuS])
    #DataItem=[SegNum,ContigLength,SiblingCount,Start,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,ChromNo,CNPriors[int(CN)] if int(CN)<len(CNPriors) else 0,Label]
    #DataItem=[ContigLength,SiblingCount,Start,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,CNPriors[int(CN)] if int(CN)<len(CNPriors) else 0,Label]
    Length=End-Start
    StartPortion=Start/ContigLength
    EndPortion=End/ContigLength
    HasSibling=1 if SiblingCount>0 else 0
    MultiSibling=1 if SiblingCount>1 else 0
    SiblingRatio=SiblingCount/SampleCount
    DataItem=[SiblingRatio,HasSibling,MultiSibling,Length,StartPortion,EndPortion,Mu,MuS,PassConfidence,CN,Confidence,CScore,CNPriors[int(CN)] if int(CN)<len(CNPriors) else 0,Label]
    #DataItem=[SiblingCount,Length,StartPortion,EndPortion,Mu,MuS,PassConfidence,CN,Confidence,CScore,CNPriors[int(CN)] if int(CN)<len(CNPriors) else 0,Label]
    DataTensor=torch.Tensor(DataItem[:-1])
    DataTensor.resuires_grad=True
    LabelTensor=torch.LongTensor([Label])
    if Label==0:
        D0+=1
    else:
        D1+=1
    Data.append(DataTensor)
    Labels.append(LabelTensor)
    D+=1

model=torch.load(SavePath)

Correct=0
Correct0=0
Correct1=0
Pred0=0
Pred1=0
Real0=0
Real1=0
for i in range(D):
    #y_pred=1 if model(Data[i])>0.5 else 0
    y_pred=model(Data[i])
    y_pred=0 if y_pred[0]>y_pred[1] else 1
    if y_pred==1:
        Pred1+=1
    else:
        Pred0+=1
    if Labels[i]==0:
        Real0+=1
    else:
        Real1+=1
    if y_pred==Labels[i]:
        Correct+=1
        if y_pred==1:
            Correct1+=1
        else:
            Correct0+=1
#print(D,D0,D1,Correct0,Correct1,"%s %s %s"%(Correct0/D0,Correct1/D1,(Correct/D)))
PR=Correct1/Pred1 if Pred1!=0 else 0
Recall=Correct1/Real1 if Real1!=0 else 0
F1=2*PR*Recall/(PR+Recall) if (PR+Recall)!=0 else 0
print(Real0,Real1,Correct0,Correct1,"PR:%s, Recall:%s, F1:%s"%(PR,Recall,F1))