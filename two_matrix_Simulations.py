# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:47:10 2019

@author: Suha Naser-Khdour
"""

import os
from ete3 import Tree
import pandas as pd
import numpy as np
import scipy.stats as st
from random import choice, sample
import time

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
    Given an empirical distributions of parameters (GTR, nucleotide frequencies, and ivariant sites' proportion)
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
        if parameter == 'A-G':
            AG = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'A-T':
            AT = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'C-G':
            CG = abs(round(dis.rvs(*dist[2][:-2], loc=dist[2][-2], scale=dist[2][-1], size=1)[0],5))
        if parameter == 'C-T':
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
    dst = Distribution()
    fitted_dist = dst.Fit(brachLenDist)
    return fitted_dist

def rand_branch_dist(branchProbabilityDist, n):
    '''
    Randomly choose a distribution for branch lengths
    return n sampled branch lengths from this distribution
    '''
    dis = getattr(st, branchProbabilityDist[0])
    branches = []
    for i in range(0,n):
        branches.append(dis.rvs(*branchProbabilityDist[2][:-2], loc=branchProbabilityDist[2][-2], scale=branchProbabilityDist[2][-1], size=1)[0])
    return branches

def topology_dist(tree, nodes, r_nodes, r_tips, branchProbabilityDist):
    '''
    Given a simulated tree and branches for inheritance
    assign branch lengths from a random distribution of branch lengths
    '''
    n = 2*len(tree)-2
    branches = rand_branch_dist(branchProbabilityDist, n)
    j = 0
    for n in tree.traverse():
        if not n.is_root():
            n.dist = abs(branches[j])
            j += 1            
    t = tree.write(format=1)
    return t

def random_tree(trees):
    '''
    Randomly choose a tree and find two nodes for inheritance
    '''
#Randomly choose a tree   
    while True:
        tree = choice(open(trees).readlines())
        t = Tree(tree, format=1)
        tips = []
        nodes = []
        k = 1
        for node in t.traverse():
            if node.is_leaf():
                tips.append(node.name)
            elif not node.is_root():
                node.add_features(name='n'+str(k))
                nodes.append(node.name)  
                k += 1
        nodes = list(filter(None, nodes))

#Randomly choose two nodes for inheritance    
        timeout1 = time.time() + 60
        timeout2 = time.time() + 90
        while True:
            rn2 = Tree(tree, format=1)
            rn = sample(nodes, 2)
            rn1 = t.search_nodes(name=rn[0])[0]
            rn2 = t.search_nodes(name=rn[1])[0]
            if time.time() <= timeout1:
                if (len(rn1.get_leaves()) <= 2) or (len(rn2.get_leaves()) <= 2):
                    continue
                elif rn2 in rn1.get_descendants():
                    continue
                elif rn1 in rn2.get_descendants():
                    continue
                elif rn2 in rn1.get_sisters():
                    continue
                else:
                    r_tips = []
                    r_nodes = []
                    for node in rn1.traverse():
                        if node.is_leaf():
                            r_tips.append(node.name)
                        else:
                            r_nodes.append(node.name)
                    root1 =  t.get_common_ancestor(r_tips)
                    root2 = []
                    for node in rn2.traverse():
                        if node.is_leaf():
                            r_tips.append(node.name)
                            root2.append(node.name)
                        else:
                            r_nodes.append(node.name)
                    root2 =  t.get_common_ancestor(root2)
                    dist = t.get_distance(root1,root2, topology_only=True)
                    tree = topology_dist(t, nodes, r_nodes, r_tips, branchProbabilityDist)
                    return [tree, nodes, tips, r_nodes, r_tips, dist]
            elif time.time() <= timeout2:
                if (len(rn1.get_leaves()) < 2) or (len(rn2.get_leaves()) < 2):
                    continue
                elif rn2 in rn1.get_descendants():
                    continue
                elif rn1 in rn2.get_descendants():
                    continue
                elif rn2 in rn1.get_sisters():
                    continue
                else:
                    r_tips = []
                    r_nodes = []
                    root1 = []
                    for node in rn1.traverse():
                        if node.is_leaf():
                            r_tips.append(node.name)
                            root1.append(node.name)
                        else:
                            r_nodes.append(node.name)
                    root1 =  t.get_common_ancestor(root1)
                    root2 = []
                    for node in rn2.traverse():
                        if node.is_leaf():
                            r_tips.append(node.name)
                            root2.append(node.name)
                        else:
                            r_nodes.append(node.name)
                    root2 =  t.get_common_ancestor(root2)
                    dist = t.get_distance(root1,root2, topology_only=True)
                    tree = topology_dist(t, nodes, r_nodes, r_tips, branchProbabilityDist)
                    return [tree, nodes, tips, r_nodes, r_tips, dist]
            else:
               break
  
def Create_Hetero_Files(site_info_file, parameter_file, param_file_list, treeFile, trees, paramProbabilityDist):
    '''
    A function that creates the input files for Hetero2 for SRH, non-SRH, and non-homogeneous conditions
    Input: the paths for all Hetero2 files, the best-fit distribution of parameters, and the dials
    '''
#determine the initial parameters
    Q0 = rand_GTR_dist(paramProbabilityDist)
    Q0 = np.round(Q0, 5)
    piRoot = np.round(Q0[-4:],5)       

#determine the other substitution model
    Q1 = rand_GTR_dist(paramProbabilityDist)
    Q1 = np.round(Q1, 5)

#write the site information file
    invariant = ['0.0'] + list(map(str,piRoot))
    invariant = ','.join(invariant)
    invariant = invariant.replace(',','\t')
    variant = ['1.0'] + list(map(str,piRoot))
    variant = ','.join(variant)
    variant = variant.replace(',','\t')

#write the site information file
    with open(site_info_file, 'w') as output_file:
        output_file.writelines('# Format: [name of site category]	[variant/invariant]	[proportion]	[freq(A)]	[freq(C)]	[freq(G)]	[freq(T)]\n')
        output_file.writelines('Constant_site	invariant\t'+invariant+'\n')
        output_file.writelines('Category_1	variant\t'+variant)
    
# Randomly choose a tree
    tree, nodes, tips, r_nodes, r_tips, dist = random_tree(trees)
    
#assign substitution model for each node and tip       
    with open(parameter_file, 'w') as output_file:
        output_file.write('#Node	S1	S2	S3	S4	S5	S6	Pi_1	Pi_2	Pi_3	Pi_4\n')
        Q = list(map(str, Q0))
        for tip in np.setdiff1d(tips,r_tips):
            output_file.writelines(tip+'\t')
            output_file.writelines('%s\t' % p  for p in Q)
            output_file.writelines('\n')
        for node in np.setdiff1d(nodes,r_nodes):
            output_file.writelines(node+'\t')
            output_file.writelines('%s\t' % p  for p in Q)
            output_file.writelines('\n')
#assign the other substitution model for each node and tip
        Q = list(map(str, Q1))
        for tip in r_tips:
            output_file.writelines(tip+'\t')
            output_file.writelines('%s\t' % p  for p in Q)
            output_file.writelines('\n')
        for node in r_nodes:
            output_file.writelines(node+'\t')
            output_file.writelines('%s\t' % p  for p in Q)
            output_file.writelines('\n')

#write the parameter list file   
    with open(param_file_list, 'w') as output_file:
        output_file.writelines('# Format: [name of variant site category]	[parameter file name]\n')
        output_file.writelines('Category_1\t')
        output_file.writelines(parameter_file)

#write the tree file
    with open (treeFile, 'w') as inf:
        inf.writelines('# Format: [name of variant site category]	[newick tree format]\n')
        inf.writelines('Category_1\t')
        inf.writelines(tree)
    return dist

if __name__ == '__main__': 
    rootDir = '/data/Suha/GTR_parameters_dist' #root directory
    paramDist = empirical_dist(pd.read_csv(os.path.join(rootDir,'GTRparam.csv'))) #GTRparam.csv contains all the empirical parameters for Q matrix and base frequencies
    treesDist = branch_dist(os.path.join(rootDir,'BranchLen.csv')) #BranchLen.csv contains the empirical branch lengths from all trees
    paramProbabilityDist = bestFit_paramDist(paramDist) #the best-fit probability distribution for each parameter
    branchProbabilityDist = bestFit_branchDist(treesDist) #the best-fit probability distribution of branch lengths for each dataset
    simulatedTrees20taxa = topology_dist(os.path.join(rootDir,'simTrees20taxa.txt'), branchProbabilityDist) #simulated trees with 20 taxa
    simulatedTrees40taxa = topology_dist(os.path.join(rootDir,'simTrees40taxa.txt'), branchProbabilityDist) #simulated trees with 40 taxa
    simulatedTrees60taxa = topology_dist(os.path.join(rootDir,'simTrees60taxa.txt'), branchProbabilityDist) #simulated trees with 60 taxa
    simulatedTrees80taxa = topology_dist(os.path.join(rootDir,'simTrees80taxa.txt'), branchProbabilityDist) #simulated trees with 80 taxa
    simulatedTrees100taxa = topology_dist(os.path.join(rootDir,'simTrees100taxa.txt'), branchProbabilityDist) #simulated trees with 100 taxa
    
    '''Hetero2 files'''
    treeFile = os.path.join(rootDir, 'tree.txt')
    site_info_file = os.path.join(rootDir, 'site_info_file.txt')
    param_file_list = os.path.join(rootDir, 'param_file_list.txt')
    parameter_file = os.path.join(rootDir, 'parameter_1.txt')

    '''Start simulation'''    
    simulatedTrees = {'20':simulatedTrees20taxa, '40':simulatedTrees40taxa, '60':simulatedTrees60taxa,
                      '80':simulatedTrees80taxa, '100':simulatedTrees100taxa}

    ls = []
    for  m, trees in simulatedTrees.items():
        for n in [100,1000,10000]:
            for c in range(0,1000):
                dist = Create_Hetero_Files(site_info_file, parameter_file, param_file_list, treeFile, trees, paramProbabilityDist)
                ls.append([m,n,c,dist])
                bashCommand = " ".join(["Hetero2", treeFile, site_info_file, param_file_list, "-l", str(n), "-o", os.path.join(rootDir,'m'+str(m)+'n'+str(n)+'_'+str(c))])
                os.system(bashCommand) #run Hetero2 software that is available on https://github.com/thomaskf/Hetero/releases
    pd.DataFrame(ls, columns=('m', 'n', 'k', 'distance')).to_csv(os.path.join(rootDir,'distances.csv'))
