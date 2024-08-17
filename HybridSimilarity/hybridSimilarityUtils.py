import os
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
#
#
#
def integrated_similarity_measure(ic , ont):
    P = None
    if os.path.exists('/home/d4gl0s/diploma/Experiment/Datasets/transitionProb.csv'):
        P = pd.read_csv('/home/d4gl0s/diploma/Experiment/Datasets/transitionProb.csv')
        P.index=ic['terms'].values.tolist()
    else:
        # 1. get transition prob for each parent-child
        tp = getTransitionProb(ic, ont)
        # 2. initialize matrix P
        P = pd.DataFrame(0, index=ic['terms'].values.tolist(), columns=ic['terms'].values.tolist(), dtype=np.float64)
        # 3. construct matrix P (filling rows)
        tcounter = 0
        for i in ic['terms'].values.tolist():# 3.1 for each column term
            tcounter+=1
            icu.progressBar(tcounter , len(ic['terms'].values.tolist())**2)
            for j in ic['terms'].values.tolist():# 3.2 for each row term
                tcounter+=1
                # 3.2.1 get roots to see if there is a reachable path
                #rooti , namespacei = findRoot(i,ont)
                #rootj , namespacej = findRoot(j,ont)
                if i==j:# 3.2.2 if same term , then value is 1
                    P.loc[i , j]=1
                #elif not (rooti==rootj):# 3.2.3 there is no reachable path
                #    P.loc[i, j]=0
                else:# 3.2.4 find transition probability from i to j
                    p = 0
                    if i in tp.keys():# roots have no transition prob
                        for anc in tp[i]:
                            if anc[0]==j:
                                p=anc[1]
                                break
                            #endif
                        #endfor
                    #endif
                    P.loc[i , j]=p
                #endif
            #endfor
        #endfor
        P.to_csv('/home/d4gl0s/diploma/Experiment/Datasets/transitionProb.csv',index=False)
    #endif
    # 4. construct matrix W rows are leaf nodes and
    identity_matrix = np.eye( len(ic['terms'].values.tolist()) )
    W = pd.DataFrame(identity_matrix, index=ic['terms'].values.tolist(), columns=ic['terms'].values.tolist(),dtype=np.float64)
    epsilon = .001
    W_star = W
    while 1 :
        W=W_star
        W_star=P*W
        if (np.linalg.norm(W_star-W, axis=1)<epsilon).any() :
            break
    #endwhile
    W=W_star
    # 5. get leaves - parents matrix
    W_ll = W.loc[list(terms)][ic['terms'].tolist()]
    # 6. Get leaves hosted similarity matrix
    HSM = None
    if os.path.exists('/home/d4gl0s/diploma/Experiment/Datasets/HSM.csv'):
        HSM = pd.read_csv('/home/d4gl0s/diploma/Experiment/Datasets/HSM.csv')
        HSM.index=ic['terms'].values.tolist()
    else:
        HSM = pd.DataFrame(0, index=ic['terms'].values.tolist(), columns=ic['terms'].values.tolist(), dtype=np.float64)
        tcounter=0
        G=ont.get_graph()
        for i in ic['terms'].values.tolist() :
            for j in ic['terms'].values.tolist() :
                tcounter+=1
                icu.progressBar(tcounter , len(ic['terms'].values.tolist())**2)
                HSM.loc[i,j] = icu.simResnikMICA(i,j,ont,ic,G=G)
            #endfor
        #endfor
        HSM.to_csv('/home/d4gl0s/diploma/Experiment/Datasets/HSM.csv')
    #endif
    RWC = (W.loc[ic['terms'].tolist()][list(terms)].dot(HSM.loc[terms][terms])).dot(W.loc[list(terms)][ic['terms'].tolist()])
    ICM = .5 * (RWC + HSM)
    return
