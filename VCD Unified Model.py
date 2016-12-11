# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 19:15:27 2016

@author: Cameron
"""

import numpy as np
import scipy.integrate
import scipy.optimize
from matplotlib import pyplot as plt

import random
import SALib as sa
import SALib.sample

from collections import namedtuple

from numpy import array

Data = namedtuple('Data', ['Exp_Number', 'VCD_start','Initial_Cl', 'Initial_K', 'Initial_Lac', 'times', 'VCD'])

def plot_data(Values):

    plt.plot(0, Values.VCD_start, 'ko')
    plt.plot(Values.times, Values.VCD,':o', label=("Exp " + Values.Exp_Number))
    plt.ylim(0,)
    plt.ylabel('Viable Cell Density (10^6 Cells/mL)')
    plt.xlabel('time (Hours)')
    #plt.legend()
    
AllData = [

### D201 Experiment 065, 067, 068 Data ###

   Data(Exp_Number="D201 E067 FL1", 
        VCD_start=0.68, 
        Initial_Cl = 114,
        Initial_K = 7.3,
        Initial_Lac = 0.65,
        times=array([ 0, 14.83,  20.61,  39.06,  64.57,  87.20,  111.38]), 
        VCD=array([ 0.68, 0.82,  0.84,  0.94,  0.66,  0.45, 0.19])), 
   Data(Exp_Number="D201 E067 FL2", 
        VCD_start=0.7,
        Initial_Cl = 114,
        Initial_K = 7.3,
        Initial_Lac = 0.65,
        times=array([ 0, 14.91,  20.51,  39.00,  64.50,  87.21,  111.35]), 
        VCD=array([ 0.7, 0.89,  0.88,  0.91,  0.70,  0.44, 0.22])),
   Data(Exp_Number="D201 E068 FL1", 
        VCD_start=0.74,
        Initial_Cl = 114,
        Initial_K = 7.3,
        Initial_Lac = 0.65,
        times=array([ 0, 18.19,  43.39,  69.12,  94.12,  112.82]), 
        VCD=array([ 0.74, 0.82,  0.79,  0.70,  0.38,  0.15])),
   Data(Exp_Number="D201 E068 FL2", 
        VCD_start=0.82,
        Initial_Cl = 114,
        Initial_K = 7.3,
        Initial_Lac = 0.65,
        times=array([ 0, 17.99,  43.39,  69.11,  94.12,  112.96]), 
        VCD=array([ 0.82, 0.97,  0.89,  0.75,  0.40,  0.16])),
   Data(Exp_Number="D201 E065 FL1", 
        VCD_start=0.66,
        Initial_Cl = 114,
        Initial_K = 7.3,
        Initial_Lac = 0.65,
        times=array([ 0, 20.99,  41.58,  67.26,  90.35,  111.93]), 
        VCD=array([ 0.66, 0.88,  0.90,  0.60,  0.21,  0.10])),
   Data(Exp_Number="D201 E065 FL2", 
        VCD_start=0.68,
        Initial_Cl = 114,
        Initial_K = 7.3,
        Initial_Lac = 0.65,
        times=array([ 0, 20.99,  41.59,  67.26,  90.35,  112.02]), 
        VCD=array([ 0.68, 0.78,  0.79,  0.56,  0.22,  0.13])),
                  
### D401 Experiment 045 Data ###                  
                  
   Data(Exp_Number="D401 E045 FL1", 
         VCD_start=2.20,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.15,
         times=array([ 0, 2.73,  19.33,  26.52,  43.35,  68.03,  94.40, 117.88]), 
         VCD=array([ 2.20, 1.98,  1.35,  3.14,  3.34,  3.57, 2.96, 2.54])), 
    Data(Exp_Number="D401 E045 FL1", 
         VCD_start=1.95,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.15,
         times=array([ 0, 2.72,  19.33,  43.36,  68.03,  94.39,  117.87]), 
         VCD=array([ 1.95, 2.07,  2.82,  3.54,  3.76,  3.45, 2.55])),
    Data(Exp_Number="D401 E045 FL1", 
         VCD_start=2.05,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.15,
         times=array([ 0, 2.73,  19.33,  43.36,  68.03,  94.39,  117.87]), 
         VCD=array([ 2.05, 2.02,  2.65,  3.52,  3.61,  2.96, 2.52])),
    Data(Exp_Number="D401 E045 FL1", 
         VCD_start=2.03,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.15,
         times=array([ 0, 2.73,  19.33,  43.36,  68.03,  94.39,  117.87]), 
         VCD=array([ 2.03, 1.97,  2.86,  3.79,  4.00,  3.74, 2.79])),
    Data(Exp_Number="D401 E045 FL1", 
         VCD_start=1.99,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.15,
         times=array([ 0, 2.73,  19.33,  43.36,  68.03,  94.39,  117.87]), 
         VCD=array([ 1.99, 2.05,  2.69,  3.75,  4.17,  3.59, 3.23])),
    Data(Exp_Number="D401 E045 FL1", 
         VCD_start=1.96, 
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.15,
         times=array([ 0, 2.73,  19.33,  43.36,  68.03,  94.39,  117.87]),
         VCD=array([ 1.96, 2.08,  2.95,  3.89,  3.96,  3.43, 2.99])),
    Data(Exp_Number="D401 E045 FL1",
         VCD_start=2.17,
         Initial_Cl = 83,
         Initial_K = 7.3,
         Initial_Lac = 0.65,
         times=array([0, 2.74, 19.33, 43.36, 68.03, 94.39, 117.87]),
         VCD=array([2.17, 2.02, 2.87, 3.66, 4.28, 3.75, 3.20])),
    Data(Exp_Number="D401 E045 FL1",
         VCD_start=2.00,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.15,
         times=array([0, 2.74, 19.33, 43.36, 68.03, 94.39, 117.87]),
         VCD=array([2.00, 2.07, 2.87, 3.64, 3.88, 3.80, 3.25])),
    Data(Exp_Number="D401 E045 FL1",
         VCD_start=2.02,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.15,
         times=array([0, 2.74, 19.33, 43.36, 68.03, 94.39, 117.87]),
         VCD=array([2.02, 2.10, 2.80, 3.69, 3.89, 3.84, 3.28])),
                   
### D401 Experiment 046 Data ###       
            
    Data(Exp_Number="D401 E046 FL4", 
         VCD_start=1.73,
         Initial_Cl = 92,
         Initial_K = 8.1,
         Initial_Lac = 1.34,
         times=array([ 0, 5.93,  24.52,  46.54,  70.48,  96.21,  122.31]), 
         VCD=array([ 1.73, 1.77,  2.48,  2.87,  3.11,  2.88, 2.57])), 
    Data(Exp_Number="D401 E046 FL5", 
         VCD_start=1.69,
         Initial_Cl = 92,
         Initial_K = 8.1,
         Initial_Lac = 1.34,
         times=array([ 0, 5.93,  24.53,  46.54,  70.49,  96.22,  122.31]), 
         VCD=array([ 1.69, 1.81,  2.41,  2.95,  2.94,  2.82, 2.59])),
    Data(Exp_Number="D401 E046 FL6", 
         VCD_start=1.66,
         Initial_Cl = 92,
         Initial_K = 8.1,
         Initial_Lac = 1.34,
         times=array([ 0, 5.93,  24.53,  46.54,  70.49,  96.22,  122.31]), 
         VCD=array([ 1.66, 1.84,  2.43,  2.80,  3.04,  2.80, 2.49])),                   

### D301 Experiment 045 Data ###

    Data(Exp_Number="D301 E049 FL1", 
         VCD_start=1.43,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.58,  49.26,  70.18,  96.23,  119.95]), 
         VCD=array([ 1.43, 2.26,  3.03,  3.24,  2.77,  2.10])), 
    Data(Exp_Number="D301 E049 FL2", 
         VCD_start=1.37,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.59,  49.28,  70.18,  96.22,  119.95]), 
         VCD=array([ 1.37, 2.21,  3.05,  3.01,  2.75,  2.12])), 
    Data(Exp_Number="D301 E049 FL3", 
         VCD_start=1.46,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.58,  49.26,  70.18,  96.23,  119.95]), 
         VCD=array([ 1.46, 2.16,  3.25,  3.42,  2.79,  2.15])), 
    Data(Exp_Number="D301 E049 FL4", 
         VCD_start=1.52,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.63,  49.32,  70.18,  96.22,  119.95]), 
         VCD=array([ 1.52, 2.15,  3.30,  3.34,  3.04,  2.13])), 
    Data(Exp_Number="D301 E049 FL5", 
         VCD_start=1.50,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.58,  49.26,  70.18,  96.23,  119.95]), 
         VCD=array([ 1.50, 2.31,  3.29,  3.12,  3.21,  2.65])), 
    Data(Exp_Number="D301 E049 FL6", 
         VCD_start=1.52,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.67,  49.36,  70.50,  96.22,  119.95]), 
         VCD=array([ 1.52, 2.38,  3.30,  3.51,  3.10,  2.49])), 
    Data(Exp_Number="D301 E049 FL7",
         VCD_start=1.43,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.70,  49.35,  70.50,  96.22,  119.95]), 
         VCD=array([ 1.43, 2.30,  3.32,  3.40,  3.25,  2.38])), 
    Data(Exp_Number="D301 E049 FL8",
         VCD_start=1.47,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.72,  49.35,  70.50,  96.22,  119.95]), 
         VCD=array([ 1.47, 2.34,  3.07,  3.40,  3.12,  2.43])), 
    Data(Exp_Number="D301 E049 FL9",
         VCD_start=1.47,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.72,  49.35,  70.51,  96.21,  119.95]), 
         VCD=array([ 1.47, 2.26,  3.09,  3.18,  2.98,  2.37])), 
    Data(Exp_Number="D301 E049 FL10",
         VCD_start=1.45,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.71,  49.34,  70.49,  96.22,  119.99]), 
         VCD=array([ 1.45, 2.39,  3.10,  3.20,  2.95,  2.35])), 
    Data(Exp_Number="D301 E049 FL11",
         VCD_start=1.43,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.66,  49.29,  70.58,  96.27,  119.95]), 
         VCD=array([ 1.43, 2.42,  3.28,  3.43,  3.22,  2.66])), 
    Data(Exp_Number="D301 E049 FL12",
         VCD_start=1.40,
         Initial_Cl = 83,
         Initial_K = 4.4,
         Initial_Lac = 1.34,
         times=array([ 0, 23.62,  49.34,  70.53,  96.23,  120.10]), 
         VCD=array([ 1.40, 2.30,  3.15,  3.29,  3.25,  2.43])),

      ]
      
for i,datum in enumerate(AllData):
    print("Experiment {}, called {}, ran for {} hours".format(i, datum.Exp_Number , datum.times[-1]))
    plot_data(datum)
    
Variables = namedtuple('Variables', ['d', 'e', 'f', 'g', 'h' , 'i' , 'j' , 'k'])

InitialGuess = Variables(
    d = 13650000. , 
    e = -2.716 , 
    f = 1.981 , 
    g = 3.9872 ,
    h = 11.193 ,
    i = 33.294 ,
    j = 0.7036 ,
    k = 0.3646 ,
    )
    
FinalValues = Variables(0,0,0,0,0,0,0,0)

StandardErrors = Variables(0,0,0,0,0,0,0,0)

M = sum((len(Data.times) for Data in AllData))

print("In total, there will be M={} x_data entries".format(M))
print("each with k=6 values; Exp_Number, Initial_Cl, Initial_K, Initial_ Lac, VCD_Start, and t")
print("and M={} y_data entries, each being a VCD.".format(M))

x_data = np.zeros((6,M))
y_data = np.zeros(M)
i=0
for Data in AllData:
    for times, VCD in zip(Data.times, Data.VCD):
        #x_data[0,i] = Values.Exp_Number
        x_data[1,i] = Values.Initial_Cl
        x_data[2,i] = Values.Initial_K
        x_data[3,i] = Values.Initial_Lac
        x_data[4,i] = Values.VCD_start
        x_data[5,i] = times
        y_data[i] = VCD
        i += 1
        
print('x_data = ',repr(x_data))
print('y_data = ',repr(y_data))

def Model(x_data,
          d,  
          e,
          f,
          g,
          h,
          i,
          j,
          k,            
         ):
             
    Exp_Numbers, Initial_Cls, Initial_Ks, Initial_Lacs, VCD_starts, ts = x_data
    M = len(VCD_starts) 
    y_data = np.zeros(M)
    for i in range(M):
        ICl = Initial_Cls[i]
        IK = Initial_Ks[i]
        ILac = Initial_Lacs[i]
        t = ts[i]
        VCD_start = VCD_starts[i]
        
        alpha = f * (np.log(IK - ILac)) - g
        beta = h * np.log(-alpha) + i
        zeta = d * ( (ICl + IK)**(e))
        gamma = 3/2 * (VCD_start) - alpha/beta*(zeta)
        upsilon = j * ((np.e)**(k * gamma))
        

        y_data[i] = upsilon + alpha * np.sqrt(1 + (t - zeta)**2/beta**2)

    return y_data
    
ReturnValues, Covariance = scipy.optimize.curve_fit(Model,
                                                   x_data,
                                                   y_data,
                                                   p0=InitialGuess,
                                                   bounds=([13649999, -3, 1.95, 3.5,10,30,0.65,0.3], [13650001, -2.5, 2, 4.5,12,35,0.75,0.4]),
                                                   method='trf')
                                                   
print('fitted',FinalValues)
stdev = np.sqrt(np.diag(Covariance))
print('+/-',stdev,'(one sigma)')

print(Covariance)
FinalValues = Variables(*ReturnValues)
print(FinalValues)

StandardErrors = Variables(*stdev)
print(StandardErrors)                                              