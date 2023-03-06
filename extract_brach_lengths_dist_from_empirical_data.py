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
    branchLenDist_leaf = []
    branchLenDist_int = []
    for n in t.traverse():
        if n.dist > 0:
            if n.is_leaf():
                branchLenDist_leaf.append(n.dist)
            else:
                branchLenDist_int.append(n.dist)
    df_leaf = pd.DataFrame({'branchLen':branchLenDist_leaf})
    df_leaf['dataset'] = os.path.basename(os.path.dirname(treefile))
    df_int = pd.DataFrame({'branchLen':branchLenDist_int})
    df_int['dataset'] = os.path.basename(os.path.dirname(treefile))
    return df_leaf, df_int


if __name__ == '__main__':
    print("In order to run this script all files must have the same name and extension"\
          " and they should be saved in directories that have the datasets name"\
              " Please see an example below")
    diagram = Tree("((----->treeFileName.treefile)----->dataset1Dir, "\
                   "(----->treeFileName.treefile)----->dataset2Dir, "\
                       "(----->treeFileName.treefile)----->dataset3Dir)rootDir;"
                       , format=1)
    print(diagram.get_ascii(show_internal=True))
    rootDir = '/MY_PATH/' #the root directory that contain the tree files
    treeFileName = 't.treefile' #the name of the tree file with .treefile extension 
    branchLenFile_leaf = 'BranchLen_leaf.csv' #the name of the leaves branch lengths output file
    branchLenFile_int = 'BranchLen_int.csv' #the name of the internal branch lengths output file 
    proceed = input("do you want to proceed? Y/N\n")

    if proceed == 'Y':
    
        df_leaf = pd.DataFrame()
        df_int = pd.DataFrame()
    
        for DirName, subdirList, fileList in os.walk(rootDir):
            if treeFileName in fileList:
                treeFile = os.path.join(DirName,treeFileName)
                df_leaf = df_leaf.append(branchLen(treeFile)[0])
                df_int = df_int.append(branchLen(treeFile)[1])
    
        df_leaf.to_csv(os.path.join(rootDir, branchLenFile_leaf))
        df_int.to_csv(os.path.join(rootDir, branchLenFile_int))
    else:
        sys.exit()
