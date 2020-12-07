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

DataFile=open(sys.argv[1],"r")
Data=[]
Labels=[]
D0=0
D1=0
D=0

for line in DataFile:
    sl=line.split()
    SegNum=int(sl[0])
    ContigLength=int(sl[1])
    SiblingCount=int(sl[2])
    Start=int(sl[3])
    End=int(sl[4])
    Mu=int(sl[5])
    MuS=int(sl[6])
    Label=int(sl[7])
    DataTensor=torch.Tensor([SegNum,ContigLength,SiblingCount,Start,End,Mu,MuS])
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
for i in range(D):
    y_pred=1 if model(Data[i])>0.5 else 0
    if y_pred==Labels[i]:
        Correct+=1
        if y_pred==1:
            Correct1+=1
        else:
            Correct0+=1
print(D,D0,D1,Correct0,Correct1,"%s %s %s"%(Correct0/D0,Correct1/D1,(Correct/D)))