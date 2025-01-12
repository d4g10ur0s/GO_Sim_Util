import os
import sys
sys.path.append('/home/d4gl0s/diploma/Experiment')
import json
import random
import time
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
def integratedSM(frequency, prob , ont):
    # 0. extract terms from gene data
    terms=frequency['terms'].values.tolist()
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
        # 3. construct matrix P (filling rows)
        tcounter = 0
        for i in terms+artNodeNames:# 3.1 for each row term
            tcounter+=1
            icu.progressBar(tcounter , len(terms+artNodeNames)**2)
            for j in terms+artNodeNames:# 3.2 for each column term
                tcounter+=1
                if i==j:# 3.3 if same term , then value is 1
                    P.loc[i , j]=1
                elif i in tProb.keys():# 3.4 find transition probability , parent is in rows
                    tParent=tProb[i]
                    if j in list(tParent.keys()):# 3.4.1 if j is a child of tParent
                        P.loc[i , j]=tParent[j]# add the transition prob
                    else:# 3.4.2 if j is not a child of tParent
                        P.loc[i , j]=0# transition prob is 0
                    #endif
                else:# 3.5 term i is not a parent term
                    P.loc[i , j]=0# transition prob is 0
                #endif
            #endfor
        #endfor
        P.to_csv('/home/d4gl0s/diploma/Experiment/Datasets/transitionProb.csv',index=False)
    #endif
    # 4. construct matrix W rows are leaf nodes and
    identity_matrix = np.eye( len(P.index) )
    W = pd.DataFrame(identity_matrix, index=P.index, columns=P.index,dtype=np.float64)
    epsilon = 1e-5
    W_star = W
    while 1 :
        W=W_star
        W_star=P*W
        if (np.linalg.norm(W_star-W, axis=1)<epsilon).any() :
            break
    #endwhile
    W=W_star
    # 5. get leaves - parents matrix
    lterms = list(set(P.columns)-set(nonL))
    # 6. Get leaves hosted similarity matrix
    HSM = None
    if os.path.exists('/home/d4gl0s/diploma/Experiment/Datasets/HSM.csv'):
        HSM = pd.read_csv('/home/d4gl0s/diploma/Experiment/Datasets/HSM.csv')
        HSM.drop(columns=['Unnamed: 0'],inplace=True)
        HSM.index=HSM.columns
    else:
        HSM = pd.DataFrame(0, index=P.index, columns=P.columns, dtype=np.float64)
        tcounter=0
        G=ont.get_graph()
        for i in HSM.index :
            for j in HSM.columns :
                tcounter+=1
                icu.progressBar(tcounter , len(HSM.index)**2)
                HSM.loc[i,j] = icu.simResnikMICA(i,j,ont,prob,G=G)
            #endfor
        #endfor
        HSM.to_csv('/home/d4gl0s/diploma/Experiment/Datasets/HSM.csv')
    #endif
    RWC = (((W.loc[lterms , list(nonL)]).T).dot(HSM.loc[lterms,lterms])).dot(W.loc[lterms , list(nonL)])
    print(RWC)
    ISM = .5 * (RWC + HSM)
    return ISM
#
#
#
def hybridRSS(t1, t2 , ont, prob, G=None):
    if G==None:
        G=ont.get_graph()
    #endif
    # 1. find MICA
    #print('Calculating alpha')
    anc1 = gu.getAllParents(t1 , ont)
    anc2 = gu.getAllParents(t2 , ont)
    # 1.2 get all ancestors
    commonAnc = set(anc1)&set(anc2)
    # 1.3 find the MOST INFORMATIVE
    mica = None
    for a in commonAnc :
        if mica==None :
            mica=a
        elif -np.log(prob[prob['terms']==mica]['probability'].values[0])<-np.log(prob[prob['terms']==a]['probability'].values[0]) :
            mica=a
        #endif
    #endfor
    alpha = None
    if mica==None:
        return 0# there are no common ancestors , 0 similarity
    else:
        alpha = -np.log(prob[prob['terms']==mica]['probability'].values[0])
    # 2. find MIL
    #print('Calculating beta')
    leaves_1 = gu.getAllLeaves(t1 , ont)
    leaves_2 = gu.getAllLeaves(t2 , ont)
    beta = None
    if len(leaves_1)==0 and len(leaves_2)==0 :
        beta=0
    elif len(leaves_1)==0 :
        # get t2's mil
        mil_2 = prob[prob['terms'].isin(leaves_2)]['probability'].max()
        mil_2=-np.log(mil_2)
        if np.isnan(mil_2) :# probability is 2 small , make beta = 1
            beta=1
        else:
            beta = ((-np.log(prob[prob['terms']==t2]['probability']) - mil_2)/2).values[0]
    elif len(leaves_2)==0 :
        # get t1's mil
        mil_1 = prob[prob['terms'].isin(leaves_1)]['probability'].max()
        mil_1=-np.log(mil_1)
        if np.isnan(mil_1) :# probability is 2 small , make beta = 1
            beta=1
        else:
            beta = ((-np.log(prob[prob['terms']==t1]['probability']) - mil_1)/2).values[0]
    else:
        # get t1's mil
        mil_1 = prob[prob['terms'].isin(leaves_1)]['probability'].max()
        if np.isnan(mil_1) :# probability is to small , make beta = 1
            mil_1=1
        else:
            mil_1=-np.log(mil_1)
        #endif
        # get t2's mil
        mil_2 = prob[prob['terms'].isin(leaves_2)]['probability'].max()
        if np.isnan(mil_2):# probability is to small , make beta = 1
            mil_2=1
        else:
            mil_2=-np.log(mil_2)
        #endif
        # calculate beta
        beta = ( ((-np.log(prob[prob['terms']==t1]['probability']) - mil_1)+(-np.log(prob[prob['terms']==t2]['probability']) - mil_2) )/2 ).values[0]
    # 3. calculate relevant distance from MICA
    #print('Calculating gamma')
    if np.isnan(beta):# probability is to small , make beta = 1
        beta=1
    #endif
    path_1 , rdist_1 = ebm.findMinimumPath(list(nx.all_simple_paths(G ,mica ,t1)))
    path_2 , rdist_2 = ebm.findMinimumPath(list(nx.all_simple_paths(G ,mica ,t2)))
    gamma = rdist_1 + rdist_2
    if abs(alpha)==0:
        return 0
    #endif
    sim_HRSS = ( 1/(1+gamma) ) * ( alpha/(alpha + beta) )
    #endif
    #print(f'HRSS of terms {t1} , {t2} : {sim_HRSS}')
    return sim_HRSS
#
#
#
def calculateHRSS(prob ,ont):
    start_time=time.time()
    G=ont.get_graph()# get ontology graph
    # 1. create the matrix if it does not exist
    if os.path.exists('/home/d4gl0s/diploma/Experiment/Datasets/HRSS.csv'):
        HRSS = pd.read_csv('/home/d4gl0s/diploma/Experiment/Datasets/HRSS.csv')
        HRSS.drop(columns=['Unnamed: 0'],inplace=True)
        HRSS.index=HRSS.columns
    else:
        HRSS = pd.DataFrame(0, index=prob['terms'].values.tolist(), columns=prob['terms'].values.tolist(), dtype=np.float64)
        # calculate hrss
        tcounter=0
        for i in HRSS.index:
            for j in HRSS.columns:
                tcounter+=1
                icu.progressBar(tcounter , len(HRSS.index)**2,start_time)
                if i==j:
                    HRSS.loc[i,j]=1
                else:
                    HRSS.loc[i,j]=hybridRSS(i, j , ont, prob , G)
            #endfor
        #endfor
        # save in csv format
        HRSS.to_csv('/home/d4gl0s/diploma/Experiment/Datasets/HRSS.csv')
    #endif
    print(HRSS)
    return HRSS
#
#
#
def calculateSemanticWeights(prob , ont):
    # 1. calculate knowledge for all terms
    start_time=time.time()
    K=prob.copy()
    K['probability']=-1/np.log(K['probability'])
    # 2. calculate semantic weights for all terms
    SW=K.copy()
    SW['probability']=1/(1+np.exp(-SW['probability']))
    # 3. calculate semantic value
    SV=K['probability'].copy()
    SV.columns=['semanticValue']
    SV.index=SW['terms'].values.tolist()
    # 3.1 for each term calculate SV[term]
    tcounter=0
    for t in SW['terms'].values.tolist():
        tcounter+=1
        icu.progressBar(tcounter , len(SW['terms'].values.tolist()), start_time)
        ps = gu.getAllParents(t , ont)
        SV.loc[t]=(SW[SW['terms'].isin(list(set(ps)))]['probability'].sum())
    #endfor
    return SV , SW
#
#
#
def calculateSemanticWeightsHybridSimilarity(prob , ont):
    start_time=time.time()
    if os.path.exists('/home/d4gl0s/diploma/Experiment/Datasets/semanticWeightsHybridSimilarity.csv'):
        semanticWeightsHybridSimilarity = pd.read_csv('/home/d4gl0s/diploma/Experiment/Datasets/semanticWeightsHybridSimilarity.csv')
        semanticWeightsHybridSimilarity.drop(columns=['Unnamed: 0'],inplace=True)
        semanticWeightsHybridSimilarity.index=semanticWeightsHybridSimilarity.columns
        return semanticWeightsHybridSimilarity
    else:
        semanticWeightsHybridSimilarity = pd.DataFrame(0, index=prob['terms'].values.tolist(), columns=prob['terms'].values.tolist(), dtype=np.float64)
        # 1. calculate semantic weights and semantic values for each term
        print(f'Calculating Weights')
        SV , SW = calculateSemanticWeights(prob , ont)
        tcounter=0
        print(f'Calculating Similarity')
        for i in semanticWeightsHybridSimilarity.index:
            for j in semanticWeightsHybridSimilarity.columns:
                tcounter+=1
                icu.progressBar(tcounter , len(semanticWeightsHybridSimilarity.index)**2,start_time)
                if i==j:
                    semanticWeightsHybridSimilarity.loc[i,j]=1
                    continue
                #endif
                # 2. get all parents for each term
                anc1=gu.getAllParents(i,ont)
                anc2=gu.getAllParents(j,ont)
                commonParents=set(anc1)&set(anc2)
                if len(commonParents)==0:# 2.1 if no common parents , then 0 similarity
                    semanticWeightsHybridSimilarity.loc[i,j]=0
                else:# 2.2 calculate semanticWeightsHybridSimilarity
                    ancSum = SW[SW['terms'].isin(list(commonParents))]['probability'].sum()
                    semanticWeightsHybridSimilarity.loc[i,j]=ancSum/(SV.loc[i]+SV.loc[j])
                #endif
            #endfor
        #endfor
        print(semanticWeightsHybridSimilarity)
        return semanticWeightsHybridSimilarity
    #endif
