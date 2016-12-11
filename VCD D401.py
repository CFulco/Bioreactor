# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 12:04:44 2016

@author: Cameron
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
    plt.ylabel('VCD (1E6 Cells/mL)')
    plt.xlabel('Time (hrs)')
    plt.legend()
    
from numpy import array
experiments1 = [ExperimentData(Run=1, 
                              VCD_start=2.20, 
                              times=array([ 0, 2.73, 26.52,  43.35,  68.03,  94.40, 117.88]), 
                              VCD=array([ 2.20, 1.98, 3.14,  3.34,  3.57, 2.96, 2.54])), 
               ExperimentData(Run=2, 
                              VCD_start=1.95, 
                              times=array([ 0, 2.72,  19.33,  43.36,  68.03,  94.39,  117.87]), 
                              VCD=array([ 1.95, 2.07,  2.82,  3.54,  3.76,  3.45, 2.55])),
               ExperimentData(Run=3, 
                              VCD_start=2.05, 
                              times=array([ 0, 2.73,  19.33,  43.36,  68.03,  94.39,  117.87]), 
                              VCD=array([ 2.05, 2.02,  2.65,  3.52,  3.61,  2.96, 2.52])),
               ExperimentData(Run=4, 
                              VCD_start=2.03, 
                              times=array([ 0, 2.73,  19.33,  43.36,  68.03,  94.39,  117.87]), 
                              VCD=array([ 2.03, 1.97,  2.86,  3.79,  4.00,  3.74, 2.79])),
               ExperimentData(Run=5, 
                              VCD_start=1.99, 
                              times=array([ 0, 2.73,  19.33,  43.36,  68.03,  94.39,  117.87]), 
                              VCD=array([ 1.99, 2.05,  2.69,  3.75,  4.17,  3.59, 3.23])),
               ExperimentData(Run=6, 
                              VCD_start=1.96, 
                              times=array([ 0, 2.73,  19.33,  43.36,  68.03,  94.39,  117.87]),
                              VCD=array([ 1.96, 2.08,  2.95,  3.89,  3.96,  3.43, 2.99])),
               ExperimentData(Run=7,
                              VCD_start=2.17,
                              times=array([0, 2.74, 19.33, 43.36, 68.03, 94.39, 117.87]),
                              VCD=array([2.17, 2.02, 2.87, 3.66, 4.28, 3.75, 3.20])),
               ExperimentData(Run=8,
                              VCD_start=2.00,
                              times=array([0, 2.74, 19.33, 43.36, 68.03, 94.39, 117.87]),
                              VCD=array([2.00, 2.07, 2.87, 3.64, 3.88, 3.80, 3.25])),
               ExperimentData(Run=9,
                              VCD_start=2.02,
                              times=array([0, 2.74, 19.33, 43.36, 68.03, 94.39, 117.87]),
                              VCD=array([2.02, 2.10, 2.80, 3.69, 3.89, 3.84, 3.28])),

              ]
              
              
for i,e in enumerate(experiments1):
    print("Experiment {}, called Run {}, ran for {} hours".format(i, e.Run, e.times[-1]))
    plot_experiment(e)
    
ParameterSet = namedtuple('ParameterSet', ['a', 'b', 'Xo', 'Yo'])

starting_guess = ParameterSet(
    a = -1 , 
    b = 40. , 
    Xo = 73. , 
    Yo = 6.  
    )
    
optimized_parameters = ParameterSet(0,0,0,0)

standard_errors = ParameterSet(0,0,0,0)

M = sum((len(e.times) for e in experiments1))
print("In total will have M={} x_data entries".format(M))
print("each with k=2 values, Run# and t")
print("and M={} y_data entries, each being a VCD.".format(M))
x_data = np.zeros((2,M))
y_data = np.zeros(M)
i=0
for e in experiments1:
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
        
        s = 73.8
        q = -1.646
        r = 39.121
        
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
                                                         #bounds=([-3, 30, 40, 0], [0, 100, 80, 10]),
                                                         method='trf')

print('fitted',optimal_parameters)
stdev = np.sqrt(np.diag(covariance))
print('+/-',stdev,'(one sigma)')

print(covariance)
optimized_parameters = ParameterSet(*optimal_parameters)
print(optimized_parameters)

standard_errors = ParameterSet(*stdev)
print(standard_errors)