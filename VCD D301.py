# -*- coding: utf-8 -*-
"""
Created on Sun Dec  4 14:17:10 2016

@author: Cameron
"""

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
    #plt.legend()
    
from numpy import array
experiments1 = [ExperimentData(Run=1, 
                              VCD_start=1.43, 
                              times=array([ 0, 23.58,  49.26,  70.18,  96.23,  119.95]), 
                              VCD=array([ 1.43, 2.26,  3.03,  3.24,  2.77,  2.10])), 
               ExperimentData(Run=2, 
                              VCD_start=1.37, 
                              times=array([ 0, 23.59,  49.28,  70.18,  96.22,  119.95]), 
                              VCD=array([ 1.37, 2.21,  3.05,  3.01,  2.75,  2.12])), 
               ExperimentData(Run=3, 
                              VCD_start=1.46, 
                              times=array([ 0, 23.58,  49.26,  70.18,  96.23,  119.95]), 
                              VCD=array([ 1.46, 2.16,  3.25,  3.42,  2.79,  2.15])), 
               ExperimentData(Run=4, 
                              VCD_start=1.52, 
                              times=array([ 0, 23.63,  49.32,  70.18,  96.22,  119.95]), 
                              VCD=array([ 1.52, 2.15,  3.30,  3.34,  3.04,  2.13])), 
               ExperimentData(Run=5, 
                              VCD_start=1.50, 
                              times=array([ 0, 23.58,  49.26,  70.18,  96.23,  119.95]), 
                              VCD=array([ 1.50, 2.31,  3.29,  3.12,  3.21,  2.65])), 
               ExperimentData(Run=6, 
                              VCD_start=1.52, 
                              times=array([ 0, 23.67,  49.36,  70.50,  96.22,  119.95]), 
                              VCD=array([ 1.52, 2.38,  3.30,  3.51,  3.10,  2.49])), 
               ExperimentData(Run=7,
                              VCD_start=1.43,
                              times=array([ 0, 23.70,  49.35,  70.50,  96.22,  119.95]), 
                              VCD=array([ 1.43, 2.30,  3.32,  3.40,  3.25,  2.38])), 
               ExperimentData(Run=8,
                              VCD_start=1.47,
                              times=array([ 0, 23.72,  49.35,  70.50,  96.22,  119.95]), 
                              VCD=array([ 1.47, 2.34,  3.07,  3.40,  3.12,  2.43])), 
               ExperimentData(Run=9,
                              VCD_start=1.47,
                              times=array([ 0, 23.72,  49.35,  70.51,  96.21,  119.95]), 
                              VCD=array([ 1.47, 2.26,  3.09,  3.18,  2.98,  2.37])), 
               ExperimentData(Run=10,
                              VCD_start=1.45,
                              times=array([ 0, 23.71,  49.34,  70.49,  96.22,  119.99]), 
                              VCD=array([ 1.45, 2.39,  3.10,  3.20,  2.95,  2.35])), 
               ExperimentData(Run=11,
                              VCD_start=1.43,
                              times=array([ 0, 23.66,  49.29,  70.58,  96.27,  119.95]), 
                              VCD=array([ 1.43, 2.42,  3.28,  3.43,  3.22,  2.66])), 
               ExperimentData(Run=12,
                              VCD_start=1.40,
                              times=array([ 0, 23.62,  49.34,  70.53,  96.23,  120.10]), 
                              VCD=array([ 1.40, 2.30,  3.15,  3.29,  3.25,  2.43])),                                         

              ]
              
              
for i,e in enumerate(experiments1):
    print("Experiment {}, called Run {}, ran for {} hours".format(i, e.Run, e.times[-1]))
    plot_experiment(e)
    
ParameterSet = namedtuple('ParameterSet', ['a', 'b', 'Xo', 'Yo'])

starting_guess = ParameterSet(
    a = -3 , 
    b = 50. , 
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
        
        s = 73
        q = -1.766
        
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

        y_data[i] = Yo + q * np.sqrt(1 + (t - s)**2/b**2)

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