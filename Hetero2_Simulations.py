# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 15:23:10 2019

@author: Suha Naser-Khdour
"""
import os
from ete3 import Tree
import pandas as pd
import numpy as np
import scipy.stats as st

class Distribution(object):
    
    def __init__(self,dist_names_list = []):
        self.dist_names = ['alpha', 'beta', 'bradford', 'chi', 'chi2', 'dgamma', 'dweibull', 'erlang', 'exponnorm', 'exponweib', 
                           'exponpow', 'gamma', 'genlogistic', 'genpareto', 'gennorm', 'genexpon', 'gengamma', 'halflogistic', 
                           'halfnorm', 'halfgennorm', 'invgamma', 'invgauss', 'invweibull', 'laplace', 'loggamma', 'logistic', 'loglaplace', 'lognorm', 
                           'maxwell', 'norm', 'pareto', 'powerlaw', 'powerlognorm', 'powernorm', 'uniform', 'weibull_max', 'weibull_min']
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
        
            
def empirical_dist(df):
    '''
    Given a dataframe (df) with empirical distributions of parameters (GTR, nucleotide frequencies, and ivariant sites' proportion)
    return a list of all parameters
    '''
    AC = [float(i) for i in df['A-C'].values]
    AG = [float(i) for i in df['A-G'].values]
    AT = [float(i) for i in df['A-T'].values]
    CG = [float(i) for i in df['C-G'].values]
    CT = [float(i) for i in df['C-T'].values]
    A = [float(i) for i in df['A'].values]
    C = [float(i) for i in df['C'].values]
    G = [float(i) for i in df['G'].values]
    T = [float(i) for i in df['T'].values]
    
    return [AC, AG, AT, CG, CT, A, C, G, T]

def bestFit_paramDist(paramDist):
    '''
    Given an empirical distributions of parameters (GTR and nucleotide frequencies)
    return a dictionary of best-fit distribution for each parameter
    '''
    parameters = {'A-C':paramDist[0], 'A-G':paramDist[1], 'A-T':paramDist[2], 'C-G':paramDist[3], 'C-T':paramDist[4], 'A':paramDist[5],
                  'C':paramDist[6], 'G':paramDist[7], 'T':paramDist[8]}
    fitted_dist = {}
    for parameter, dist in parameters.items():
        dst = Distribution()
        fitted_dist[parameter] = dst.Fit(dist)
    return fitted_dist

def rand_GTR_dist(paramProbabilityDist):
    '''
    Randomly choose a rate matrix (S0) and base frequencies (pi0) for the SRH alignment   
    '''
    for parameter, dist in paramProbabilityDist.items():
        dis = getattr(st, dist[0])
        if parameter == 'A-C':
            AC = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
            while AC > 100:
                AC = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'A-G':
            AG = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
            while AG > 100:
                AG = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'A-T':
            AT = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
            while AT > 100:
                AT = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'C-G':
            CG = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
            while CG > 100:
                CG = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'C-T':
            CT = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
            while CT > 100:
                CT = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'A':
            A = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'C':
            C = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'G':
            G = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'T':
            T = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
    freqSum = A + C + G + T
    Q = np.array([AC,AG,AT,CG,CT,1,A/freqSum,C/freqSum,G/freqSum,T/freqSum])
    return Q

def branch_dist(f):
    '''
    Given an empirical distributions of trees (topology + branch length)
    return a list of the non-zero branch lengths
    '''
    df = pd.read_csv(f)
    dist = [float(i) for i in df['branchLen'].values]
    return dist

def bestFit_branchDist(brachLenDist):
    '''
    Given an empirical distributions of branch lengths from different datasets
    return a dictionary of best-fit distribution for each dataset
    '''
    td = [i for i in brachLenDist if i>0.00001]
    dst = Distribution()
    fitted_dist = dst.Fit(td)
    return fitted_dist

def rand_branch_dist(branchProbabilityDist, n):
    '''
    Randomly choose a distribution for branch lengths
    return n sampled branch lengths from this distribution
    '''
    dis = getattr(st, branchProbabilityDist[0])
    branches = []
    for i in range(0,n):
        b = dis.rvs(*branchProbabilityDist[2][:-2], loc=branchProbabilityDist[2][-2], scale=branchProbabilityDist[2][-1], size=1)[0]
        while b <= 0:
            b = dis.rvs(*branchProbabilityDist[2][:-2], loc=branchProbabilityDist[2][-2], scale=branchProbabilityDist[2][-1], size=1)[0]
        branches.append(b)
    return branches

def topology_dist(simTrees, branchProbabilityDist):
    '''
    Given a file with simulated trees
    Fro each tree, assign branch lengths from a random distribution of branch lengths and return 
    '''
    trees1 = []
    trees2 = []
    with open(simTrees, 'r') as inf:
        for i, line in enumerate(inf):
            t= Tree(line)
            n = 2*len(t)-2
            branches = rand_branch_dist(branchProbabilityDist, n)
            j = 0
            k = 1
            for n in t.traverse():
                if not n.is_root():
                    if not n.is_leaf():
                        n.add_features(name="n"+str(k))
                        k += 1
                    n.dist = branches[j]
                    j += 1
            if i < 3630: #Assuming that we simulated 3960 trees for each number of taxa
                trees1.append(t.write(format=1))
            else:
                trees2.append(t.write(format=1))
    return [trees1, trees2]
        
  
def Create_Hetero_Files(site_info_file, parameter_file, param_file_list, tree, paramProbabilityDist, m, n, v, w, main=True):
    '''
    A function that creates the input files for Hetero2 for SRH, non-SRH, and non-homogeneous conditions
    Input: the paths for all Hetero2 files, the best-fit distribution of parameters, and the dials
    '''
#determine the initial parameters
    Q0 = rand_GTR_dist(paramProbabilityDist)
    S0 = Q0[:-4]
    pi0 = Q0[-4:]

#determine the root
    Q = rand_GTR_dist(paramProbabilityDist)
    if main:
        piRoot = np.round(pi0*(v/10) + Q[-4:]*(1-(v/10)),5)
        pii = np.round(pi0,5)
    else:
        piRoot = np.round(pi0,5)

#write the site information file
    invariant = ['0.0'] + list(map(str,piRoot))
    invariant = ','.join(invariant)
    invariant = invariant.replace(',','\t')
    variant = ['1.0'] + list(map(str,piRoot))
    variant = ','.join(variant)
    variant = variant.replace(',','\t')

    with open(site_info_file, 'w') as output_file:
        output_file.writelines('# Format: [name of site category]	[variant/invariant]	[proportion]	[freq(A)]	[freq(C)]	[freq(G)]	[freq(T)]\n')
        output_file.writelines('Constant_site	invariant\t'+invariant+'\n')
        output_file.writelines('Category_1	variant\t'+variant)
    
#write the parameter file
    t = Tree(tree, format=1)
    tips = []
    nodes = []
    for node in t.traverse():
        if node.is_leaf():
            tips.append(node.name)
        elif not node.is_root():
            nodes.append(node.name)

    nodes = list(filter(None, nodes))
    if main:
        with open(parameter_file, 'w') as output_file:
            output_file.write('#Node	S1	S2	S3	S4	S5	S6	Pi_1	Pi_2	Pi_3	Pi_4\n')
            for tip in tips:
                Q = rand_GTR_dist(paramProbabilityDist)
                S = S0*(w/10) + Q[:-4]*(1-(w/10))
                Q = np.append(S,pi0)
                Q = np.round(Q, 5)
                Q = list(map(str, Q))
                output_file.writelines(tip+'\t')
                output_file.writelines('%s\t' % p  for p in Q)
                output_file.writelines('\n')
            for node in nodes:
                Q = rand_GTR_dist(paramProbabilityDist)
                S = S0*(w/10) + Q[:-4]*(1-(w/10))
                Q = np.append(S,pi0)
                Q = np.round(Q, 5)
                Q = list(map(str, Q))
                output_file.writelines(node+'\t')
                output_file.writelines('%s\t' % p  for p in Q)
                output_file.writelines('\n')
    else:
        with open(parameter_file, 'w') as output_file:
            output_file.write('#Node	S1	S2	S3	S4	S5	S6	Pi_1	Pi_2	Pi_3	Pi_4\n')
            for tip in tips:
                Q = rand_GTR_dist(paramProbabilityDist)
                S = S0*(w/10) + Q[:-4]*(1-(w/10))
                pii = pi0*(v/10) + Q[-4:]*(1-(v/10))
                Q = np.append(S,pii)
                Q = np.round(Q, 5)
                Q = list(map(str, Q))
                output_file.writelines(tip+'\t')
                output_file.writelines('%s\t' % p  for p in Q)
                output_file.writelines('\n')
            for node in nodes:
                Q = rand_GTR_dist(paramProbabilityDist)
                S = S0*(w/10) + Q[:-4]*(1-(w/10))
                pii = pi0*(v/10) + Q[-4:]*(1-(v/10))
                Q = np.append(S,pii)
                Q = np.round(Q, 5)
                Q = list(map(str, Q))
                output_file.writelines(node+'\t')
                output_file.writelines('%s\t' % p  for p in Q)
                output_file.writelines('\n')

#write the parameter list file   
    with open(param_file_list, 'w') as output_file:
        output_file.writelines('# Format: [name of variant site category]	[parameter file name]\n')
        output_file.writelines('Category_1\t')
        output_file.writelines(parameter_file)
    
    return

if __name__ == '__main__': 
    rootDir = '/data/Suha/GTR_parameters_dist'
    paramDist = empirical_dist(pd.read_csv(os.path.join(rootDir,'GTRparam.csv'))) #param_dist.csv contains all the empirical parameters for Q matrix, base frequencies and %invariant sites
    treesDist = branch_dist(os.path.join(rootDir,'BranchLen.csv')) #BranchLen.csv contains the brach lengths from all empirical trees
    paramProbabilityDist = bestFit_paramDist(paramDist) #the best-fit probability distribution for each parameter
    branchProbabilityDist = bestFit_branchDist(treesDist) #the best-fit probability distribution of branch lengths for each dataset
    simulatedTrees20taxa = topology_dist(os.path.join(rootDir,'simTrees20taxa.txt'), branchProbabilityDist) #simulated trees with branch lengths
    simulatedTrees40taxa = topology_dist(os.path.join(rootDir,'simTrees40taxa.txt'), branchProbabilityDist) #simulated trees with branch lengths
    simulatedTrees60taxa = topology_dist(os.path.join(rootDir,'simTrees60taxa.txt'), branchProbabilityDist) #simulated trees with branch lengths
    simulatedTrees80taxa = topology_dist(os.path.join(rootDir,'simTrees80taxa.txt'), branchProbabilityDist) #simulated trees with branch lengths
    simulatedTrees100taxa = topology_dist(os.path.join(rootDir,'simTrees100taxa.txt'), branchProbabilityDist) #simulated trees with branch lengths
    
    '''Hetero2 files'''
    treeFile = os.path.join(rootDir, 'tree.txt')
    site_info_file = os.path.join(rootDir, 'site_info_file.txt')
    param_file_list = os.path.join(rootDir, 'param_file_list.txt')
    parameter_file = os.path.join(rootDir, 'parameter_1.txt')

    '''Start simulation'''    
    simulatedTrees1 = {'20':simulatedTrees20taxa[0], '40':simulatedTrees40taxa[0], '60':simulatedTrees60taxa[0],
                      '80':simulatedTrees80taxa[0], '100':simulatedTrees100taxa[0]}
    simulatedTrees2 = {'20':simulatedTrees20taxa[1], '40':simulatedTrees40taxa[1], '60':simulatedTrees60taxa[1],
                      '80':simulatedTrees80taxa[1], '100':simulatedTrees100taxa[1]}

    for  m, trees in simulatedTrees1.items():
        i = 0
        for n in [100,1000,10000]:      
            for w in range(0,11):
                for v in range(0,11):
                    for c in range(0,10):
                        with open (treeFile, 'w') as inf:
                            inf.writelines('# Format: [name of variant site category]	[newick tree format]\n')
                            inf.writelines('Category_1\t')                       
                            inf.writelines(trees[i])
                        Create_Hetero_Files(site_info_file, parameter_file, param_file_list, trees[i], paramProbabilityDist, m, n, v, w, main=True)
                        bashCommand = " ".join(["Hetero2", treeFile, site_info_file, param_file_list, "-l", str(n), "-o", os.path.join(rootDir,'m'+str(m)+'n'+str(n)+'w'+str(w)+'v'+str(v)+'_'+str(c))])
                        os.system(bashCommand) #run Hetero2 software that is available on https://github.com/thomaskf/Hetero/releases
                        i += 1
    for  m, trees in simulatedTrees2.items():
        i = 0
        for n in [100,1000,10000]:
            for p in range(0,11):
                for c in range(0,10):
                    with open (treeFile, 'w') as inf:
                        inf.writelines('# Format: [name of variant site category]	[newick tree format]\n')
                        inf.writelines('Category_1\t')
                        inf.writelines(trees[i])
                    Create_Hetero_Files(site_info_file, parameter_file, param_file_list, trees[i], paramProbabilityDist, m, n, p, p, main=False)
                    i += 1
                    bashCommand = " ".join(["Hetero2", treeFile, site_info_file, param_file_list, "-l", str(n), "-o", os.path.join(rootDir,'m'+str(m)+'n'+str(n)+'p'+str(p)+'_'+str(c))])
                    os.system(bashCommand) #run Hetero2software that is available on https://github.com/thomaskf/Hetero/releases
