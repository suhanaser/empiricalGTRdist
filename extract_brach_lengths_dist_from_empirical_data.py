# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 11:38:28 2019

@author: Suha Naser-Khdour
"""

import pandas as pd
import os
from ete3 import Tree
import sys

def branchLen(treefile):
    '''
    Creats a dataframe with all branch lengths for each dataset
    This function works with tree files in newick format 
    '''
    t = Tree(treefile)
    branchLenDist = []
    for n in t.traverse():
        if n.dist > 0:
            branchLenDist.append(n.dist)
    df = pd.DataFrame({'branchLen':branchLenDist})
    df['dataset'] = os.path.basename(os.path.dirname(treefile))
    return df


if __name__ == '__main__':
    print("In order to run this script all files must have the same name and extension and they should be saved in directories that have the datasets name. Please see an example below")
    diagram = Tree("((----->treeFileName.treefile)----->dataset1Dir, (----->treeFileName.treefile)----->dataset2Dir, (----->treeFileName.treefile)----->dataset3Dir)rootDir;", format=1)
    print(diagram.get_ascii(show_internal=True))
    rootDir = '/data/Suha/GTR_parameters_dist' #the rootDir name to the directories that contain the tree files
    treeFileName = 'branches.treefile' #the name of the tree file with .treefile extension (any newick format file can be used)
    branchLenFile = 'BranchLen.csv' #the name of the branch lengths output file with .csv extension
    proceed = input("do you want to proceed? Y/N\n")

    if proceed == 'Y':
    
        df = pd.DataFrame()
    
        for DirName, subdirList, fileList in os.walk(rootDir):
            if treeFileName in fileList:
                treeFile = os.path.join(DirName,treeFileName)
                df = df.append(branchLen(treeFile))
    
        df.to_csv(os.path.join(rootDir, branchLenFile))
    else:
        sys.exit()
