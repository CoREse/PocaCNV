import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt

import sys

SavePath="data/ScoringTrainModelData"

d_in=7
d_out=1
H=20

N=10000
D=0
DataFile=open(sys.argv[1],"r")
Data=[]
Labels=[]
D0=0
D1=0

Bads=np.arange(1000,50000)
np.random.shuffle(Bads)

BadNum=260
Bads=Bads[:BadNum]
Bads=np.sort(Bads)

BadCount=0

ValidateN=10
Selected0=0
Selected1=0
ValidateData=[]
ValidateLabel=[]
i=0
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
        if Selected0<ValidateN:
            if i%1000==347:
                ValidateData.append(DataTensor)
                ValidateLabel.append(LabelTensor)
                Selected0+=1
                continue
        if BadCount>=260:
            i+=1
            continue
        if D0!=Bads[BadCount]:
            i+=1
            continue
        BadCount+=1
    else:
        D1+=1
        if Selected1<ValidateN:
            if D1%10==7:
                ValidateData.append(DataTensor)
                ValidateLabel.append(LabelTensor)
                Selected1+=1
                continue
    Data.append(DataTensor)
    Labels.append(LabelTensor)
    D+=1
    i+=1

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
        y_pred = self.linear2(torch.sigmoid(self.linear1(x)))
        return y_pred


Epsilon=10e-10
Target=0.9
print("Target: ",Target)

while 1:
    # 初始化model
    model = TwoLayerNet(d_in, H, d_out)

    # 规定loss function
    loss_fn = nn.MSELoss(reduction='sum')

    # 规定学习率
    learning_rate = 1e-4

    #  定义optimizer做优化
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    loss_list = []
    AveLossList=[]
    LastAveLoss=0

    for it in range(N):
        AveLoss=0
        for i in range(D):
            #  forward pass
            y_pred = model(Data[i])

            #  compute loss
            loss = loss_fn(y_pred.float(), Labels[i].float())
            AveLoss+= loss.item()
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
    for i in range(Selected0+Selected1):
        y_pred=1 if model(ValidateData[i])>0.5 else 0
        if y_pred==ValidateLabel[i]:
            Correct+=1
            if y_pred==1:
                Correct1+=1
            else:
                Correct0+=1
    print(Selected0,Selected1,Correct,Correct0,Correct1)
    CR=Correct/(Selected0+Selected1)
    if CR>Target:
        break
    print("CR:",CR)
print(model.state_dict())

torch.save(model,SavePath)
#plt.plot(range(200), loss_list)
#plt.show()