import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt

import sys

SavePath="data/ScoringTrainWholeModelData"

#device=torch.device("cuda:1")
device=torch.device("cpu")

d_in=None
d_out=2
H=100
H2=100
H3=100
H4=100
H5=100
HiddenLayersN=20

# 规定学习率
learning_rate = 1e-3

BATCH_SIZE=1000
MaxBlock=20000

Even=False

ValidateTurn=50
SaveAlong=True

ValidatePortion=0.05

N=10000000
D=0
Data=[]
D0=0
D1=0

Data0=[]
Data1=[]
CNPriors=[1.4818669642343562e-06, 0.0024241621501567266, 0.9914186177260093, 0.00550095102642564, 0.0006000125666241987, 4.0389660750366106e-05, 1.1302181811873978e-05, 2.3256457943417974e-06, 5.854614618138028e-07, 7.124415377290356e-08, 1.0010930675953374e-07, 3.210044872566197e-10, 3.565148390433508e-11, 2.993635914964604e-12, 7.20618688429682e-13, 1.3465236956844422e-13, 3.0109007332307176e-14, 3.4497272100312798e-15, 2.5051103812841497e-15]

#SampleCount=int(sys.argv[1])
for DataFileName in sys.argv[1:]:
    DataFile=open(DataFileName,"r")
    ChromNo=int(DataFileName.split("_")[1])
    for line in DataFile:
        sl=line.split()
        SegNum=int(sl[0])
        ContigLength=int(sl[1])
        SiblingCount=int(sl[2])
        SiblingRatio=float(sl[3])
        Start=int(sl[4])
        End=int(sl[5])
        Mu=int(sl[6])
        MuS=int(sl[7])
        PassConfidence=float(sl[8])
        CN=float(sl[9])
        Confidence=float(sl[10])
        CScore=float(sl[11])
        HasSDRP=float(sl[12])
        HasMultiSDRP=float(sl[13])
        SDRPRatio=float(sl[14])
        Label=int(sl[15])
        #DataItem=[SegNum,ContigLength,SiblingCount,Start,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,ChromNo]+CNPriors+[Label]
        #DataItem=[SegNum,ContigLength,SiblingCount,Start,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,ChromNo,CNPriors[int(CN) if int(CN)<len(CNPriors) else -1],Label]
        #DataItem=[SegNum,ContigLength,SiblingCount,Start,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,ChromNo,CNPriors[int(CN)] if int(CN)<len(CNPriors) else 0,Label]
        #DataItem=[ContigLength,SiblingCount,Start,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,CNPriors[int(CN)] if int(CN)<len(CNPriors) else 0,Label]
        Length=End-Start
        StartPortion=Start/ContigLength
        EndPortion=End/ContigLength
        HasSibling=1 if SiblingCount>0 else 0
        MultiSibling=1 if SiblingCount>1 else 0
        DataItem=[SiblingRatio,HasSibling,MultiSibling,Length,StartPortion,EndPortion,Mu,MuS,PassConfidence,CN,Confidence,CScore,CNPriors[int(CN)] if int(CN)<len(CNPriors) else 0,HasSDRP,HasMultiSDRP,SDRPRatio,Label]
        #print(DataItem)
        #DataItem=[ContigLength,SiblingCount,Mu,MuS,PassConfidence,CN,Confidence,CScore,ChromNo,Label]
        #DataTensor=torch.Tensor([SegNum,ContigLength,SiblingCount,Start,End,Mu,MuS,PassConfidence,CN,Confidence,CScore,Label])
        #DataTensor.resuires_grad=True
        #LabelTensor=torch.LongTensor([Label])
        if Label==0:
            D0+=1
            Data0.append(DataItem)
        else:
            D1+=1
            Data1.append(DataItem)
        D+=1
    DataFile.close()

d_in=len(Data0[0])-1

ValidateClass0Size=int(D0*ValidatePortion)
ValidateClass1Size=int(D1*ValidatePortion)
if Even:
    TrainClass0Size=min(D0-ValidateClass0Size,D1-ValidateClass1Size)
    TrainClass1Size=TrainClass0Size
else:
    TrainClass0Size=D0-ValidateClass0Size
    TrainClass1Size=D1-ValidateClass1Size
ValidateSize=ValidateClass0Size+ValidateClass1Size

np.random.shuffle(Data0)
np.random.shuffle(Data1)

ValidateData=Data0[:ValidateClass0Size]+Data1[:ValidateClass1Size]
for i in range(len(ValidateData)):
    ValidateData[i]=torch.Tensor(ValidateData[i])
np.random.shuffle(ValidateData)
NPValidate=np.ndarray((ValidateSize,d_in+1))
for i in range(len(ValidateData)):
    NPValidate[i]=ValidateData[i]
TensorValidate=torch.from_numpy(NPValidate).to(device=device)

Data0=Data0[ValidateClass0Size:ValidateClass0Size+TrainClass0Size]
Data1=Data1[ValidateClass1Size:ValidateClass1Size+TrainClass1Size]
Data=Data0+Data1
#np.random.shuffle(Data)

NPData0=np.ndarray((len(Data0),d_in))
NPLabels0=np.ndarray((len(Data0)))
NPData1=np.ndarray((len(Data1),d_in))
NPLabels1=np.ndarray((len(Data1)))
NPData=np.ndarray((len(Data),d_in))
NPLabels=np.ndarray((len(Data)))
for i in range(len(Data0)):
    for j in range(d_in):
        NPData0[i][j]=Data0[i][j]
    NPLabels0[i]=Data0[i][-1]
for i in range(len(Data1)):
    for j in range(d_in):
        NPData1[i][j]=Data1[i][j]
    NPLabels1[i]=Data1[i][-1]

if BATCH_SIZE==0:
    BATCH_SIZE=min(len(Data0),len(Data1))

TensorData0=torch.from_numpy(NPData0).type(torch.FloatTensor).to(device=device)
TensorLabels0=torch.from_numpy(NPLabels0).type(torch.LongTensor).to(device=device)
TensorData1=torch.from_numpy(NPData1).type(torch.FloatTensor).to(device=device)
TensorLabels1=torch.from_numpy(NPLabels1).type(torch.LongTensor).to(device=device)
TensorData=torch.from_numpy(NPData).type(torch.FloatTensor).to(device=device)
TensorLabels=torch.from_numpy(NPLabels).type(torch.LongTensor).to(device=device)

train_set0 = torch.utils.data.TensorDataset(TensorData0,TensorLabels0)
train_set1 = torch.utils.data.TensorDataset(TensorData1,TensorLabels1)

if BATCH_SIZE!=0:
    train_loader0 = torch.utils.data.DataLoader(
            dataset=train_set0,
            batch_size=BATCH_SIZE,
            shuffle=True
            )
    train_loader1 = torch.utils.data.DataLoader(
            dataset=train_set1,
            batch_size=BATCH_SIZE,
            shuffle=True
            )

class TwoClassTrainLoader:
    def __init__(self, tl0,tl1):
        self.tl0=tl0
        self.tl1=tl1
    def __iter__(self):
        return self
    def __next__(self):
        step,(x1,y1)=next(self.tl0)
        step,(x2,y2)=next(self.tl1)
        return step,(x1+x2,y1+y2)
    def getEnu(self):
        Enu=[]
        for step,(x,y) in enumerate(self.tl0):
            Enu.append((step,(x,y)))
        i=0
        for step,(x,y) in enumerate(self.tl0):
            if i<len(Enu):
                Enu[i]=(step,(Enu[i][1][0]+x,Enu[i][1][1]+y))
            i+=1
        return Enu

train_loader=TwoClassTrainLoader(train_loader0,train_loader1)


print("Label1:%s Label0:%s, Data:%s, BATCH_SIZE:%s"%(D1,D0,len(Data),BATCH_SIZE))
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
        y_pred = self.linear2(torch.relu(self.linear1(x)))
        return y_pred


Epsilon=0#10e-10
Target=0.2
LossTarget=0#10e-10
print("Target: ",Target)

def CustomLoss(y_pred,y):
    loss=nn.CrossEntropyLoss(reduction='sum')(y_pred,y)
    return loss
    loss=loss+1/(torch.relu(0.693-loss)+0.1)
    return loss

def validate(model):
    #print("TrainingSet:")
    Correct=0
    Correct0=0
    Correct1=0
    Pred0=0
    Pred1=0
    if False:
        if len(TensorData)<=MaxBlock:
            y_preds=model(TensorData)
            for i in range(TrainClass0Size+TrainClass1Size):
                #y_pred=model(ValidateData[i][:-1])
                y_pred=y_preds[i]
                y_pred=0 if y_pred[0]>y_pred[1] else 1
                #y_pred=1 if model(ValidateData[i][:-1])>0.5 else 0
                if y_pred==0:
                    Pred0+=1
                else:
                    Pred1+=1
                if y_pred==TensorLabels[i]:
                    Correct+=1
                    if y_pred==1:
                        Correct1+=1
                    else:
                        Correct0+=1
        else:
            Start=0
            while Start<len(TensorData):
                #y_preds=torch.cat((y_preds,model(TensorData[Start:Start+MaxBlock])))
                y_preds=model(TensorData[Start:Start+MaxBlock])
                for i in range(len(y_preds)):
                    #y_pred=model(ValidateData[i][:-1])
                    y_pred=y_preds[i]
                    y_pred=0 if y_pred[0]>y_pred[1] else 1
                    #y_pred=1 if model(ValidateData[i][:-1])>0.5 else 0
                    if y_pred==0:
                        Pred0+=1
                    else:
                        Pred1+=1
                    if y_pred==TensorLabels[Start+i]:
                        Correct+=1
                        if y_pred==1:
                            Correct1+=1
                        else:
                            Correct0+=1
                Start+=MaxBlock
        print(TrainClass0Size,TrainClass1Size,Correct,Correct0,Correct1)
        CR=Correct/(TrainClass0Size+TrainClass1Size)
        #PR0=Correct0/(Pred0) if Pred0!=0 else 0
        #PR1=Correct1/(Pred1) if Pred1!=0 else 0
        #Recall0=Correct0/TrainClass0Size if TrainClass0Size!=0 else 0
        #Recall1=Correct1/TrainClass1Size if TrainClass1Size!=0 else 0
        #F10=2*Recall0*PR0/(Recall0+PR0) if (Recall0+PR0)!=0 else 0
        #F11=2*Recall1*PR1/(Recall1+PR1) if (Recall1+PR1)!=0 else 0
        #FF=2*F10*F11/(F10+F11) if (F10+F11)!=0 else 0
        PR=Correct1/Pred1 if Pred1!=0 else 0
        Recall=Correct1/TrainClass1Size if TrainClass1Size!=0 else 0
        F1=2*PR*Recall/(PR+Recall) if PR+Recall!=0 else 0
        print("CR:%s, PR:%s, Recall:%s, F1:%s"%(CR,PR,Recall,F1))
        if PR>Target and Recall>Target:
            return True
    #print("CR:%s, PR0:%s, PR1:%s, Recall0:%s, Recall1:%s, F1_0:%s, F1_1:%s, FF:%s"%(CR,PR0,PR1,Recall0,Recall1,F10,F11,FF))
    Correct=0
    Correct0=0
    Correct1=0
    Pred0=0
    Pred1=0
    if len(TensorValidate)<=MaxBlock:
        y_preds=model(TensorValidate[:,:-1].float())
        for i in range(ValidateSize):
            #y_pred=model(ValidateData[i][:-1])
            y_pred=y_preds[i]
            y_pred=0 if y_pred[0]>y_pred[1] else 1
            #y_pred=1 if model(ValidateData[i][:-1])>0.5 else 0
            if y_pred==0:
                Pred0+=1
            else:
                Pred1+=1
            if y_pred==TensorValidate[i,-1].long():
                Correct+=1
                if y_pred==1:
                    Correct1+=1
                else:
                    Correct0+=1
    else:
        Start=0
        while Start<len(TensorValidate):
            #y_preds=torch.cat((y_preds,model(TensorValidate[Start:Start+MaxBlock].float())))
            y_preds=model(TensorValidate[Start:Start+MaxBlock,:-1].float())
            for i in range(len(y_preds)):
                #y_pred=model(ValidateData[i][:-1])
                y_pred=y_preds[i]
                y_pred=0 if y_pred[0]>y_pred[1] else 1
                #y_pred=1 if model(ValidateData[i][:-1])>0.5 else 0
                if y_pred==0:
                    Pred0+=1
                else:
                    Pred1+=1
                if y_pred==TensorValidate[Start+i,-1].long():
                    Correct+=1
                    if y_pred==1:
                        Correct1+=1
                    else:
                        Correct0+=1
            Start+=MaxBlock
    print(ValidateClass0Size,ValidateClass1Size,Correct,Correct0,Correct1)
    CR=Correct/(ValidateSize)
    #PR0=Correct0/(Pred0) if Pred0!=0 else 0
    #PR1=Correct1/(Pred1) if Pred1!=0 else 0
    #Recall0=Correct0/ValidateClass0Size if ValidateClass0Size!=0 else 0
    #Recall1=Correct1/ValidateClass1Size if ValidateClass1Size!=0 else 0
    #F10=2*Recall0*PR0/(Recall0+PR0) if (Recall0+PR0)!=0 else 0
    #F11=2*Recall1*PR1/(Recall1+PR1) if (Recall1+PR1)!=0 else 0
    #FF=2*F10*F11/(F10+F11) if (F10+F11)!=0 else 0
    PR=Correct1/Pred1 if Pred1!=0 else 0
    Recall=Correct1/ValidateClass1Size if ValidateClass1Size!=0 else 0
    F1=2*PR*Recall/(PR+Recall) if (PR+Recall!=0) else 0
    print("CR:%s, PR:%s, Recall:%s, F1:%s"%(CR,PR,Recall,F1))
    if SaveAlong:
        torch.save(model,SavePath)
    if PR>Target and Recall>Target:
        return True
    #print("CR:%s, PR0:%s, PR1:%s, Recall0:%s, Recall1:%s, F1_0:%s, F1_1:%s, FF:%s"%(CR,PR0,PR1,Recall0,Recall1,F10,F11,FF))
    return False

#import swats
while 1:
    # 初始化model
    #model = TwoLayerNet(d_in, H, d_out)
    #model=nn.Sequential(nn.Linear(d_in,H),nn.ReLU(),nn.Linear(H,H2),nn.ReLU(),nn.Linear(H2,H3),nn.ReLU(),nn.Linear(H3,H4),nn.ReLU(),nn.Linear(H4,H5),nn.ReLU(),nn.Linear(H5,d_out))
    #exec("model=nn.Sequential(nn.Linear(d_in,H),nn.ReLU(),nn.Linear(H,H2),nn.ReLU(),nn.Linear(H2,H3),nn.ReLU(),nn.Linear(H3,H4),nn.ReLU(),nn.Linear(H4,H5),nn.ReLU(),nn.Linear(H5,d_out))")
    ModelString="model=nn.Sequential(nn.Linear(d_in,H),nn.ReLU()"
    for i in range(1,HiddenLayersN):
        ModelString+=",nn.Linear(H,H),nn.ReLU()"
    ModelString+=",nn.Linear(H,d_out))"
    exec(ModelString)
    """model=nn.Sequential(nn.Linear(d_in,H),nn.ReLU(),nn.Linear(H,H)\
        ,nn.ReLU(),nn.Linear(H,H)\
            ,nn.ReLU(),nn.Linear(H,H)\
                ,nn.ReLU(),nn.Linear(H,H)\
                    ,nn.ReLU(),nn.Linear(H,H)\
                        ,nn.ReLU(),nn.Linear(H,H)\
                            ,nn.ReLU(),nn.Linear(H,H)\
                                ,nn.ReLU(),nn.Linear(H,H)\
                                    ,nn.ReLU(),nn.Linear(H,H)\
                                        ,nn.ReLU(),nn.Linear(H,H)\
                                            ,nn.ReLU(),nn.Linear(H,H)\
                                                ,nn.ReLU(),nn.Linear(H,d_out))"""
    #model=torch.load("ScoringTrain/Model032/Model032.pickle")
    #model=torch.load("ScoringTrain/Model0179/Model0179.pickle")
    #model=torch.load("data/ScoringTrainModelData")
    model=torch.load(SavePath)
    #model=model.to(device=device)
    #validate(model)
    #exit(0)
    # 规定loss function
    #loss_fn = nn.BCELoss(reduction='sum')
    if Even:
        loss_fn = nn.CrossEntropyLoss(reduction='sum')
    else:
        loss_fn = nn.CrossEntropyLoss(reduction='sum',weight=torch.Tensor([D1/(D0+D1),D0/(D0+D1)]))
    loss_fn=loss_fn.to(device=device)
    #loss_fn=CustomLoss
    #loss_fn=nn.MSELoss(reduction='sum')


    #  定义optimizer做优化
    if Even:
        Decay=0.1
    else:
        Decay=0.1
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate,weight_decay=Decay)
    #optimizer=torch.optim.SGD(model.parameters(), lr=10e-7, momentum=0.001)
    #optimizer = swats.SWATS(model.parameters())
    #optimizer = torch.optim.Adagrad(model.parameters(), lr=learning_rate)
    #optimizer = torch.optim.ASGD(model.parameters(), lr=learning_rate)
    #optimizer = torch.optim.RMSprop(model.parameters(), lr=learning_rate)
    loss_list = []
    AveLossList=[]
    LastAveLoss=0
    LastLoss=0
    
    LastEpsilon=100

    Changed=True
    Changed2=True
    Changed3=True

    for it in range(N):
        AveLoss=0
        if BATCH_SIZE!=0:
            #for step,(x,y) in train_loader.getEnu():
            enu0=iter(train_loader0)
            enu1=iter(train_loader1)
            while 1:
                try:
                    (x0,y0)=next(enu0)
                    (x1,y1)=next(enu1)
                except:
                    break
                x=torch.cat((x0,x1))
                y=torch.cat((y0,y1))
                #print(step,x.shape)
                #y_pred=model(TensorData.float())
                y_pred=model(x)
                #y_preds+=y_pred.tolist()
                #print(y_pred.shape,TensorLabels.shape)
                #loss=loss_fn(y_pred,TensorLabels.long())
                loss=loss_fn(y_pred,y)
                AveLoss+=float(loss)
                optimizer.zero_grad()
                loss.backward()
                #print(loss)
                optimizer.step()
                #for l in loss:
                #    AveLoss+=float(l)
                '''
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
                    optimizer.step()'''
                #AveLoss/=D
                #LastLoss=loss
                #LastAveLoss=AveLoss
                #AveLossList.append(AveLoss)
        else:
            y_pred=model(TensorData)
            loss=loss_fn(y_pred,TensorLabels)
            AveLoss+=float(loss)
            optimizer.zero_grad()
            loss.backward()
            #print(loss)
            optimizer.step()
        AveLoss/=min(len(Data1),len(Data0))*2
        print(it,AveLoss)
        if not Changed and AveLoss<0.55:
            #optimizer=torch.optim.SGD(model.parameters(), lr=learning_rate,momentum=0.1)
            #print("change optimizer to SGDM!",file=sys.stderr)
            optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate/10)
            Changed=True
        if not Changed2 and AveLoss<0.35:
            optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate/100)
            Changed2=True
        if not Changed3 and AveLoss<0.2:
            optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate/1000)
            Changed3=True
        if abs(AveLoss-LastAveLoss)<Epsilon:
            if LastEpsilon<Epsilon:
                break
        LastEpsilon=abs(AveLoss-LastAveLoss)
        if AveLoss<LossTarget:
            break
        LastAveLoss=AveLoss
        if it%ValidateTurn==0:
            if validate(model):
                break
    if validate(model):
        break
    break

print(model.state_dict())

torch.save(model,SavePath)
#plt.plot(range(200), loss_list)
#plt.show()