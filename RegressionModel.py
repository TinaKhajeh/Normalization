#Author: Tina Khajeh
#Implementation of regression model with weightd cost function 

import numpy as np
import numpy as pd

def WeightedlinearRegression(x_train, x_test, y_train, epochs, learningRate):
    teta1=0
    teta0 =0 
    N = len(y_train)
    mean_x = float(sum(x_train))/len(x_train)
    meanList = [mean_x for _ in xrange(N)]
    weights = np.array([data for data in (x_train-meanList)])
    for i in range(epochs):
        y_current = (teta1*x_train)+teta0 
        cost = sum(np.array([data**2 for data in (y_train-y_current)])*np.array(weights))/N
        teta1_gradient = -(2.0/N)*sum(x_train*(y_train - y_current)*weights)
        teta0_gradient = -(2.0/N)*sum((y_train - y_current)*weights)
        teta1 = teta1 -(learningRate *teta1_gradient)
        teta0 = teta0 -(learningRate *teta0_gradient)
        
    predict_y = []
    for entry in x_test:
            predict = teta1*entry+teta0
            predict_y.extend(predict)    
     
    return np.array(predict_y)[:,np.newaxis]

    
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
    predictedMatrix [:,i] = WeightedlinearRegression(trainMic[:,i], testMic[:,i], trainRNA[:,i], 1000, 0.0001)
np.savetxt('test.txt', predictedMatrix , delimiter='\t')
