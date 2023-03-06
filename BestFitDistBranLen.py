# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 12:44:55 2019

@author: Suha Naser-Khdour
"""

from scipy import stats as st
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

class Distribution(object):

    def __init__(self,dist_names_list = []):
        self.dist_names = ['alpha', 'beta', 'bradford', 'chi', 'chi2', 'dgamma',
                           'dweibull', 'erlang', 'exponnorm', 'exponweib', 
                           'exponpow', 'gamma', 'genlogistic', 'genpareto',
                           'gennorm', 'genexpon', 'gengamma', 'halflogistic', 
                           'halfnorm', 'halfgennorm', 'invgamma', 'invgauss',
                           'invweibull', 'laplace', 'loggamma', 'logistic',
                           'loglaplace', 'lognorm', 'maxwell', 'norm', 'pareto',
                           'powerlaw', 'powerlognorm', 'powernorm', 'uniform',
                           'weibull_max', 'weibull_min']
        self.dist_results = []
        self.params = {}

        self.DistributionName = ""
        self.PValue = 0
        self.Param = None

        self.isFitted = False


    def Fit(self, y):
            self.dist_results = []
            self.params = {}
            for dist_name in self.dist_names:
                dist = getattr(st, dist_name)
                param = dist.fit(y)
                
                self.params[dist_name] = param
                #Applying the Kolmogorov-Smirnov test
                D, p = st.kstest(y, dist_name, args=param);
                self.dist_results.append((dist_name,p))
    
            #select the best fitted distribution
            sel_dist,p = (max(self.dist_results,key=lambda item:item[1]))
            #store the name of the best fit and its p value
            self.DistributionName = sel_dist
            self.Param = self.params[sel_dist]
            self.PValue = p
            
            self.isFitted = True
            return self.DistributionName,self.PValue,self.Param

    def Random(self, n = 1):
        if self.isFitted:
            dist_name = self.DistributionName
            param = self.params[dist_name]
            #initiate the scipy distribution
            dist = getattr(st, dist_name)
            return dist.rvs(*param[:-2], loc=param[-2], scale=param[-1], size=n)
        else:
            raise ValueError('Must first run the Fit method.')

    def Plot(self,y,path):
        x = self.Random(n=len(y))
        hist, bins = np.histogram(y, bins=8)
        logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
        plt.hist(x, alpha=0.5, label='fitted', range=[0,max(y)], bins=logbins)
        plt.hist(y, alpha=0.5, label='empirical', range=[0,max(y)], bins=logbins)
        plt.xscale('log')
        plt.legend(loc='upper right')
        plt.savefig(path, dpi=600)

def BranchLenDist(BranchLen):
    '''
    Parameters
    ----------
    BranchLen : path
        the csv file that contains the empirical distribution
    Returns
    ----------
    fitted_dist : the best fitted distribution
    '''
    df = pd.read_csv(BranchLen)
    emp_dist = [float(i) for i in df['branchLen'].values if i > 0.00001]
    dst = Distribution()
    dst.Fit(emp_dist)
    dst.Plot(emp_dist,os.path.splitext(BranchLen)[0]+'.png')
    return 

if __name__ == '__main__':
    BranchLen = '\BranchLen_leaf.csv'
    BranchLenDist(BranchLen)