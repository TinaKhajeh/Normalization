#Author : Tina Khajeh
import numpy as np
import numpy as pd
from math import exp, sqrt, pi 

#adaptive bandwidth recognize the length of bandwidth for kernel regression  
def adaptiveBandwidth(train_x, center, numOfPoints):
    if(numOfPoints > len(train_x) ):
        return 0;
    else:
        tmp = train_x
        tmp[:] = [abs(x-center) for x in tmp]
        tmp.sort()
        return tmp[numOfpoints] - tmp[0]
#Kernel Regression function
def kernelRegression (test_x, train_x, train_y,  kernelType):
    kernels = {
        'box': lambda x: 1/2 if (x<=1 and x>=-1) else 0,
        'gs': lambda x: 1/sqrt(2*pi)*exp(-x**2/2),
        'ep': lambda x: 3/4*(1-x**2) if (x<=1 and x>=-1) else 0
    }
    predict_y = []
    for entry in test_x:
        bandwidth = adaptiveBandwidth(train_x, entry, 10)
        nks = [np.sum((j-entry)**2)/bandwidth for j in train_x]
        ks = [kernels[kernelType](i) for i in nks]
        dividend = sum([ks[i]*train_y[i] for i in range(len(ks))])
        devisor = sum(ks)
        predict = dividend / divisor
        predict_y.extend(predict)
    return np.array(predict_y)[:,np.newaxis]

#inputs are from files shuffledMicData and shuffledRNAData which are the outputs of ExtractingTrainTestCaseAndControlSamples.R
micMatrix = np.matrix(pd.read_csv("/Users/Tina/Desktop/Thesis Code/shuffledMicData1_22.csv"))
RNAMatrix = np.matrix(pd.read_csv("/Users/Tina/Desktop/Thesis Code/shuffledRNAData1_22.csv"))
geneNumber= micMatrix.shape[1]

trainMic = micMatrix[0:95,]
trainRNA = RNAMatrix[0:95,]

testMic = micMatrix[95:104,]
testRNA = RNAMatrix[95:104,]

Nx = testRNA.shape[0]
Ny = testRNA.shape[1]
mylist= np.zeros((Nx,Ny)).tolist()
predictedMatrix = np.matrix(mylist)
predictedMatrixrid = np.matrix(mylist)

for i in range(0,(geneNumber)):
    predictedMatrix [:,i] = kernelRegression(testMic[:,i], trainMic[:,i], trainRNA[:,i],  'gs')
     
np.savetxt('test.txt', predictedMatrix , delimiter='\t')
