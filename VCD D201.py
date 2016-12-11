# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 11:09:10 2016

@author: VCDmeron
"""

import numpy as np
import scipy.integrate
import scipy.optimize
from matplotlib import pyplot as plt

import random
import SALib as sa
import SALib.sample




#try:
#    import seaborn as sns
#except ImportError:
#    # This block will be run if there's an ImportError, i.e you don't have seaborn installed.
#    sns = False
#    print ("If you want to try different figure formatting, "
#           "type 'conda install seaborn' at an anaconda command prompt or terminal. "
#           "See https://stanford.edu/~mwaskom/software/seaborn/ for details")
#    # If not using seaborn, we can still control the size of the figures this way
#    from pylab import rcParams
#    rcParams['figure.figsize'] = 3, 3
#else:
#    # This block will be run if there is no ImportError
#    sns.set_style("ticks")
#    sns.set_context("paper",rc={"figure.figsize": (2, 2)})





from collections import namedtuple
ExperimentData = namedtuple('ExperimentData', ['Run', 'VCD_start', 'times', 'VCD'])

def plot_experiment(e):
    """
    Plots the experimental data provided in 'e' 
    which should be of the type ExperimentData.
    """
    plt.plot(0, e.VCD_start, 'ko')
    plt.plot(e.times, e.VCD,':o', label="Run {:.0f}".format(e.Run))
    plt.ylim(0,)
    plt.ylabel('$C_A$ (mol/L)')
    plt.xlabel('time (s)')
    plt.legend()
    
from numpy import array
experiments = [ExperimentData(Run=1, 
                              VCD_start=0.68, 
                              times=array([ 0, 14.83,  20.61,  39.06,  64.57,  87.20,  111.38]), 
                              VCD=array([ 0.68, 0.82,  0.84,  0.94,  0.66,  0.45, 0.19])), 
               ExperimentData(Run=2, 
                              VCD_start=0.7, 
                              times=array([ 0, 14.91,  20.51,  39.00,  64.50,  87.21,  111.35]), 
                              VCD=array([ 0.7, 0.89,  0.88,  0.91,  0.70,  0.44, 0.22])),
               ExperimentData(Run=3, 
                              VCD_start=0.74, 
                              times=array([ 0, 18.19,  43.39,  69.12,  94.12,  112.82]), 
                              VCD=array([ 0.74, 0.82,  0.79,  0.70,  0.38,  0.15])),
               ExperimentData(Run=4, 
                              VCD_start=0.82, 
                              times=array([ 0, 17.99,  43.39,  69.11,  94.12,  112.96]), 
                              VCD=array([ 0.82, 0.97,  0.89,  0.75,  0.40,  0.16])),
               ExperimentData(Run=5, 
                              VCD_start=0.66, 
                              times=array([ 0, 20.99,  41.58,  67.26,  90.35,  111.93]), 
                              VCD=array([ 0.66, 0.88,  0.90,  0.60,  0.21,  0.10])),
               ExperimentData(Run=6, 
                              VCD_start=0.68, 
                              times=array([ 0, 20.99,  41.59,  67.26,  90.35,  112.02]), 
                              VCD=array([ 0.68, 0.78,  0.79,  0.56,  0.22,  0.13])),

              ]
              
for i,e in enumerate(experiments):
    print("Experiment {}, called Run {}, ran for {} hours".format(i, e.Run, e.times[-1]))
    plot_experiment(e)
    
ParameterSet = namedtuple('ParameterSet', ['a', 'b', 'Xo', 'Yo'])

starting_guess = ParameterSet(
    a = -1. , 
    b = 60. , 
    Xo = 30. , 
    Yo = 1.  
    )
    
optimized_parameters = ParameterSet(0,0,0,0)

standard_errors = ParameterSet(0,0,0,0)

M = sum((len(e.times) for e in experiments))
print("In total will have M={} x_data entries".format(M))
print("each with k=6 values; Exp_Number, Initial_Cl, Initial_K, Initial_ Lac, VCD_Start, and t")
print("and M={} y_data entries, each being a VCD.".format(M))
x_data = np.zeros((2,M))
y_data = np.zeros(M)
i=0
for e in experiments:
    for time, VCD in zip(e.times, e.VCD):
        x_data[0,i] = e.VCD_start
        x_data[1,i] = time
        y_data[i] = VCD
        i += 1
print('x_data = ',repr(x_data))
print('y_data = ',repr(y_data))

#cA_start = 10.
def my_model(x_data,
          a,  # /s 
          b,  # kJ/mol
          Xo, # kJ/mol
          Yo,  # J/mol/K
         ):
    VCD_starts, ts = x_data
    M = len(VCD_starts) # number of data points
    y_data = np.zeros(M)
    for i in range(M):
        t = ts[i]
        VCD_start = VCD_starts[i]
        
        s = 30.568
        q = -0.2243
        r = 16.577
        
        #R = 8.314 # J/mol/K
        #kf = 10**logA * np.exp(-Ea*1e3 / (R * T))
        #dG = dH*1e3 - T * dS # J/mol
        #Ka = np.exp(-dG / (R * T))
        #Kc = Ka # ideal solution, unit activity is 1 M
        #kr = kf / Kc
#        def dcAdt(cA, t):
#            cB = cA_start - cA
#            return kf * cB - kr * cA
        #result = scipy.integrate.odeint(dcAdt, cA_start, [0,t])
        #cA = result[-1,0]

        y_data[i] = Yo + q * np.sqrt(1 + (t - s)**2/r**2)

    return y_data

#my_model(np.array([[298],[ 0, 14.83,  20.61,  39.06,  64.57,  87.20,  111.38]]), 7,50,-10,-40,)

optimal_parameters, covariance = scipy.optimize.curve_fit(my_model,
                                                          x_data,
                                                          y_data,
                                                         p0=starting_guess,
                                                         method='trf')

print('fitted',optimal_parameters)
stdev = np.sqrt(np.diag(covariance))
print('+/-',stdev,'(one sigma)')

print(covariance)
optimized_parameters = ParameterSet(*optimal_parameters)
print(optimized_parameters)

standard_errors = ParameterSet(*stdev)
print(standard_errors)





