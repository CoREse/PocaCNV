import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt

import sys

SavePath="data/ScoringTrainModelData"

d_in=7
d_out=1
H=20

ValidatePortion=0.05

N=10000
D=0
DataFile=open(sys.argv[1],"r")
Data=[]
D0=0
D1=0

Data0=[]
Data1=[]

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
    DataTensor=torch.Tensor([SegNum,ContigLength,SiblingCount,Start,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,Label])
    DataTensor.resuires_grad=True
    LabelTensor=torch.LongTensor([Label])
    if Label==0:
        D0+=1
        Data0.append(DataTensor)
    else:
        D1+=1
        Data1.append(DataTensor)
    D+=1

TrainClassSize=min(D0,D1)
ValidateClassSize=int(TrainClassSize*ValidatePortion)
TrainClassSize-=ValidateClassSize

np.random.shuffle(Data0)
np.random.shuffle(Data1)

ValidateData=Data0[:ValidateClassSize]+Data1[:ValidateClassSize]
np.random.shuffle(ValidateData)

Data=Data0[ValidateClassSize:ValidateClassSize+TrainClassSize]+Data1[ValidateClassSize:ValidateClassSize+TrainClassSize]
np.random.shuffle(Data)

print("Label1:%s Label0:%s, Data:%s"%(D1,D0,len(Data)))
#print(Data,Labels)
#  随机创建一些训练数据
#x = torch.randn(N, d_in, requires_grad=True)
#y = torch.randn(N, d_out, requires_grad=True)

class TwoLayerNet(torch.nn.Module):
    #  定义了init和forward就想到于给定了一个模型
    #  init相当于给了这个框架
    def __init__(self, D_in, H, D_out):
        super(TwoLayerNet, self).__init__()
        self.linear1 = nn.Linear(D_in, H)
        self.linear2 = nn.Linear(H, D_out)

    #  forward是
    def forward(self, x):
        y_pred = torch.sigmoid(self.linear2(torch.sigmoid(self.linear1(x))))
        return y_pred


Epsilon=10e-10
Target=0.9
print("Target: ",Target)

while 1:
    # 初始化model
    model = TwoLayerNet(d_in, H, d_out)

    # 规定loss function
    loss_fn = nn.BCELoss(reduction='sum')

    # 规定学习率
    learning_rate = 1e-3

    #  定义optimizer做优化
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    loss_list = []
    AveLossList=[]
    LastAveLoss=0

    for it in range(N):
        AveLoss=0
        for i in range(TrainClassSize):
            #  forward pass
            y_pred = model(Data[i][:-1])
            #print(y_pred,Data[i][-1:])
            #exit(0)

            #  compute loss
            loss = loss_fn(y_pred, Data[i][-1:])
            #print(loss)
            #exit(0)
            AveLoss+= float(loss)
            #loss_list.append(loss)

            #  backward pass
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        AveLoss/=D
        print(it,AveLoss)
        if abs(AveLoss-LastAveLoss)<Epsilon:
            break
        LastAveLoss=AveLoss
        AveLossList.append(AveLoss)

    Correct=0
    Correct0=0
    Correct1=0
    for i in range(2*ValidateClassSize):
        y_pred=1 if model(ValidateData[i][:-1])>0.5 else 0
        if y_pred==ValidateData[i][-1]:
            Correct+=1
            if y_pred==1:
                Correct1+=1
            else:
                Correct0+=1
    print(ValidateClassSize,Correct,Correct0,Correct1)
    CR=Correct/(2*ValidateClassSize)
    if CR>Target:
        break
    print("CR:",CR)
print(model.state_dict())

torch.save(model,SavePath)
#plt.plot(range(200), loss_list)
#plt.show()