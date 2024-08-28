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
from Parsing import graphUtilities as gu
#
#
#
def getTransitionProb(nonL , frequency , ont):
    # 1. find transition probabilities
    terms = frequency['terms'].values.tolist()
    transitionProb = {}
    tcounter = 0
    for t in terms:
        tcounter+=1
        icu.progressBar(tcounter, len(frequency['terms'].values.tolist()))
        # 1.1 add artificial node
        if t in nonL:
            # 1.2 get all children in dataset
            validChildren = set(terms)&set(list(ont.children(t)))
            childFrequency = sum(frequency[frequency['terms'].isin(list(validChildren))]['frequency'].values.tolist())
            # 1.3 create artificial node
            pFrequency = frequency[frequency['terms']==t]['frequency'].values[0]
            tProbArt = (pFrequency-childFrequency)/pFrequency# transition probability for artificial node
            transitionProb[t]={}
            transitionProb[t][str(t)+'_art']=tProbArt# store node and transition probability for each child
            # 1.4 transition probability for children
            for c in list(validChildren):
                cFrequency=frequency[frequency['terms']==c]['frequency'].values[0]
                transitionProb[t][c]=(1 - tProbArt) * (cFrequency/childFrequency)# add in dictionary
            #endfor
        else:
            pass# do nothing for children
        #endif
        print(transitionProb)
    #endfor
    return transitionProb
#
#
#
def getNonLeaves(terms, ont):
    nonL = []
    # 1. for each term
    for t in terms:
        # 2. get immidiate children
        imChildren = list(ont.children(t))
        if len(imChildren)==0:
            continue
        #endif
        if len(list(set(terms)&set(imChildren)))==0:
            continue
        #endif
        nonL.append(t)
    #endfor
    return nonL
#
#
#
def integratedSM(frequency , prob , ont):
    # 0. extract terms from gene data
    terms=prob['terms'].values.tolist()
    # 1. find non-leaf nodes
    nonL = getNonLeaves(terms, ont)
    P = None
    # 2. get transition prob for each parent-child & add artificial nodes
    if os.path.exists('/home/d4gl0s/diploma/Experiment/Datasets/transitionProb.csv'):
        # 2.1 transition probability is already calculated
        P = pd.read_csv('/home/d4gl0s/diploma/Experiment/Datasets/transitionProb.csv')
        P.index=P.columns
    else:
        # 2.1 transition probability has not yet calculated
        # get transition probabilities and artificial nodes based on non-leaves
        tProb = getTransitionProb(nonL, frequency , ont)
        # 2.2 initialize matrix P
        artNodeNames = [str(i)+'_art' for i in nonL]# get artificial nodes names
        P = pd.DataFrame(0, index=terms+artNodeNames, columns=terms+artNodeNames, dtype=np.float64)# P is a n x n matrix
        print(P)
        # 3. construct matrix P (filling rows)
        tcounter = 0
        for i in terms+artNodeNames:# 3.1 for each row term
            tcounter+=1
            icu.progressBar(tcounter , len(terms+artNodeNames)**2)
            for j in terms+artNodeNames:# 3.2 for each column term
                tcounter+=1
                if i==j:# 3.3 if same term , then value is 1
                    P.loc[i , j]=1
                else:# 3.4 find transition probability , parent is in rows
                    tParent=tProb[j]
                    if j in list(tParent.keys()):# 3.4.1 if j is a child of tParent
                        P.loc[i , j]=tParent[j]# add the transition prob
                    else:# 3.4.2 if j is not a child of tParent
                        P.loc[i , j]=0# transition prob is 0
                    #endif
                #endif
            #endfor
            print(P)
        #endfor
        P.to_csv('/home/d4gl0s/diploma/Experiment/Datasets/transitionProb.csv',index=False)
    #endif
    # 4. construct matrix W rows are leaf nodes and
    identity_matrix = np.eye( len(P.index) )
    W = pd.DataFrame(identity_matrix, index=P.index, columns=P.index,dtype=np.float64)
    epsilon = .001
    W_star = W
    while 1 :
        W=W_star
        W_star=P*W
        if (np.linalg.norm(W_star-W, axis=1)<epsilon).any() :
            break
    #endwhile
    W=W_star
    print(W)
    # 5. get leaves - parents matrix
    W_ll = W.loc[list(lterms)][ic['terms'].tolist()]
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
    RWC = (W.loc[ic['terms'].tolist()][list(lterms)].dot(HSM.loc[lterms][lterms])).dot(W.loc[list(lterms)][ic['terms'].tolist()])
    ISM = .5 * (RWC + HSM)
    return ISM
#
#
#
def hybridRSS(t1, t2 , ont, ic, G=None):
    if G==None:
        G=ont.get_graph()
    #endif
    # 1. find MICA
    print('Calculating alpha')
    dists_1 , anc1 , namespace1 = gu.allAncestors(t1 , ont , G=None)
    dists_2 , anc2 , namespace2 = gu.allAncestors(t1 , ont , G=None)
    # 1.2 get all ancestors
    commonAnc = set(anc1)&set(anc2)
    # 1.3 find the MOST INFORMATIVE
    mica = None
    for a in commonAnc :
        if mica==None :
            mica=a
        elif ic[ic['terms']==mica]['IC'].values[0]<ic[ic['terms']==a]['IC'].values[0] :
            mica=a
        #endif
    #endfor
    alpha = ic[ic['terms']==mica]['IC'].values[0]
    # 2. find MIL
    print('Calculating beta')
    dist_1 , n1 = gu.findAllChildrenInGraph(t1 , ont)
    dist_2 , n2 = gu.findAllChildrenInGraph(t2 , ont)
    beta = None
    if n1==0 and n2==0 :
        beta=0
    elif n1==0 :
        # get t2's mil
        index_2 = list(dist_2.keys())[-1]
        leaves_2 = dist_2[index_2]
        mil_2 = ic[ic['terms'].isin(leaves_2)]['IC'].max()
        if np.isnan(mil_2) :# some leaves do not exist , because of handling subontology
            beta=0
        else:
            beta = ((ic[ic['terms']==t2]['IC'] - mil_2)/2).values[0]
    elif n2==0 :
        # get t1's mil
        index_1 = list(dist_1.keys())[-1]
        leaves_1 = dist_1[index_1]
        mil_1 = ic[ic['terms'].isin(leaves_1)]['IC'].max()
        if np.isnan(mil_1) :# some leaves do not exist , because of handling subontology
            beta=0
        else:
            beta = ((ic[ic['terms']==t1]['IC'] - mil_1)/2).values[0]
    else:
        # get t1's mil
        index_1 = list(dist_1.keys())[-1]
        leaves_1 = dist_1[index_1]
        mil_1 = ic[ic['terms'].isin(leaves_1)]['IC'].max()
        if np.isnan(mil_1) :# some leaves do not exist , because of handling subontology
            mil_1=0
        # get t2's mil
        index_2 = list(dist_2.keys())[-1]
        leaves_2 = dist_2[index_2]
        mil_2 = ic[ic['terms'].isin(leaves_2)]['IC'].max()
        if np.isnan(mil_2):# some leaves do not exist , because of handling subontology
            mil_2=0
        # calculate beta
        beta = (((ic[ic['terms']==t1]['IC'] - mil_1)+(ic[ic['terms']==t2]['IC'] - mil_2))/2).values[0]
    # 3. calculate relevant distance from MICA
    print('Calculating gamma')
    G = ont.get_graph()
    path_1 , rdist_1 = ebm.findMinimumPath(mica , t1 , G)
    path_2 , rdist_2 = ebm.findMinimumPath(mica , t2 , G)
    gamma = rdist_1 + rdist_2
    sim_HRSS = ( 1/(1+gamma) ) * ( alpha/(alpha + beta) )
    print(f'HRSS of terms {t1} , {t2} : {sim_HRSS}')
    return sim_HRSS
