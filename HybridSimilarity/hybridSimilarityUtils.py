import sys
sys.path.append('/home/d4gl0s/diploma/Experiment')
import json
import random
# data analysis tools
import pandas as pd
import numpy as np
# graph analysis tools
import networkx as nx
# ontology analysis tools
import ontobio as ob
# parallel programming
from dask.distributed import Client,LocalCluster,progress
import dask
# custom modules
from EdgeBasedSimilarity import edgeBasedMethods as ebm
from InformationContentSimilarity import informationContentUtility as icu
#
#
#
def getTransitionProb(ic, ont):
    # 1. find transition probabilities
    transitionProb = {}
    tcounter = 0
    for i in ic['terms'].values.tolist():
        tcounter+=1
        icu.progressBar(tcounter, len(ic['terms'].values.tolist()))
        children = ont.children(i)
        validChildren = set(children)&set(ic['terms'].values.tolist())# SOME CHILDREN ARE XCLUDED BECAUSE THEY DONT EXIST IN ANNOTATION DATASET
        if len(validChildren)==0:# IF LEAF NODE , THEN NO CHILDREN
            transitionProb[i] = []
            continue
        #endif
        #  filter the DataFrame for children
        children_df = ic[ic['terms'].isin(validChildren)]
        # calculate the total frequency of children
        Nsum = children_df['frequency'].sum()
        # calculate the parent frequency without children
        pfrequency = ic[ic['terms'] == i]['frequency'].values[0]
        pfrequencySub = pfrequency - Nsum
        # calculate transition probabilities for each valid child
        for c in validChildren :
            if c in transitionProb.keys():
                transitionProb[c].append( (i , ( (1 - pfrequencySub/pfrequency) * (ic[ic['terms'] == c]['frequency']/Nsum) ).values[0]) )
            else:
                transitionProb[c]=[(i , ( (1 - (pfrequencySub/pfrequency) ) * (ic[ic['terms'] == c]['frequency']/Nsum) ).values[0]) , ]
        #endfor
    #endfor
    '''
    for i in transitionProb.keys():
        for item in transitionProb[i]:
            print(f'Transition Probabilities from term {i} to term {item[0]} : {item[1]}')
        #endfor
    #endfor
    '''
    return transitionProb
