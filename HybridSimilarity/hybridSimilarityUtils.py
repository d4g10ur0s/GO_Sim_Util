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
    for i in ic.transpose().columns:
        children = ont.children(i)
        validChildren = set(children)&set(ic.transpose().columns)# SOME CHILDREN ARE XCLUDED BECAUSE THEY DONT EXIST IN ANNOTATION DATASET
        if len(validChildren)==0:# IF LEAF NODE , THEN NO CHILDREN
            transitionProb[i] = []
            continue
        #endif
        # find parent frequency without its childrens frequency
        Nsum = sum([ic['frequency'][c] for c in validChildren])
        pfrequency = ic['frequency'][i] - Nsum
        # calculate transition probabilities for each valid child
        for c in validChildren :
            # debug
            if (1 - (pfrequency/ic['frequency'][i]) ) * (ic['frequency'][c]/Nsum) > 1.1:
                print(sum([ic['frequency'][c] for c in validChildren]))
                print(pfrequency)
                print((1 - (pfrequency/ic['frequency'][i]) ))
                input( (ic['frequency'][c]/Nsum) )
            if i in transitionProb.keys():
                transitionProb[i].append( (c , (1 - pfrequency/ic['frequency'][i]) * (ic['frequency'][c]/Nsum) ) )
            else:
                transitionProb[i]=[(c , (1 - (pfrequency/ic['frequency'][i]) ) * (ic['frequency'][c]/Nsum) ) , ]
        #endfor
    #endfor
    for i in transitionProb.keys():
        for item in transitionProb[i]:
            print(f'Transition Probabilities from term {i} to term {item[0]} : {item[1]}')
        #endfor
    #endfor
    return transitionProb
