#Author : Tina Khajeh
#This code will load data as input file and predict normalized values of test samples based on learned regression model on 
#training data output will be stored in test.txt file
from sklearn import linear_model
from scipy.stats.stats import pearsonr   
from sklearn.metrics import mean_squared_error
import pandas as pd
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import PolynomialFeatures


micMatrix = np.matrix(pd.read_csv("/Users/Tina/Desktop/Thesis Code/shuffledMicData1_22.csv"))
RNAMatrix = np.matrix(pd.read_csv("/Users/Tina/Desktop/Thesis Code/shuffledRNAData1_22.csv"))
geneNumber= micMatrix.shape[1]

poly = PolynomialFeatures(degree=2)
ridReg = Ridge(alpha=0.5)

trainMic = micMatrix[0:95,]
trainRNA = RNAMatrix[0:95,]

testMic = micMatrix[95:104,]
testRNA = RNAMatrix[95:104,]

reg = linear_model.LinearRegression()

Nx = testRNA.shape[0]
Ny = testRNA.shape[1]
mylist= np.zeros((Nx,Ny)).tolist()
predictedMatrix = np.matrix(mylist)

predictedMatrixrid = np.matrix(mylist)

coeficents = []

for i in range(0,(geneNumber)):
    reg.fit(trainMic[:,i],trainRNA[:,i])
    coeficents.append(reg.coef_)
    predictedMatrix [:,i] = reg.predict(testMic[:,i])
    
    ridReg.fit(poly.fit_transform(trainMic[:,i]),trainRNA[:,i])
    predictedMatrixrid [:,i] = ridReg.predict(poly.fit_transform(testMic[:,i]))    
    
MicError =  mean_squared_error(trainRNA, trainMic)
predictError = mean_squared_error(testRNA, predictedMatrixrid)

correlationMicRNA = pearsonr((testMic[2,:].tolist())[0],(testRNA[2,:].tolist())[0])[0]
correlationPredictRNA = pearsonr((predictedMatrix[2,:].tolist())[0],(testRNA[2,:].tolist())[0])[0]

#save to file csv file : 
np.savetxt('test.txt', predictedMatrix , delimiter='\t')
np.savetxt('test2.txt', predictedMatrixrid , delimiter='\t')
