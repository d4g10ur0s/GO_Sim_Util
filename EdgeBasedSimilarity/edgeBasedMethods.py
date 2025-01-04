import os
import sys
import json
import time
# data analysis modules
import numpy as np
import pandas as pd
# ontology modules
import ontobio as ob
import networkx as nx
# custom modules
from Parsing import graphUtilities as gu
from Parsing import parsing as pu
from InformationContentSimilarity import informationContentUtility as icu
#
#
#
def getTValues(t , root , ont , tValues , G=None):
    # 0. preprocessing steps
    tValues[root]=1# 0.1 root has t-value 1
    if G==None :
        G=ont.get_graph()
    #endif
    # 1. find all ancestors to root term
    pdist = {}
    dist = 1
    parents = ont.parents(t)
    while len(parents)>0:
        pdist[dist]=parents
        dist+=1
        tparents=[]
        for p in parents:
            tparents+=list(ont.parents(p))
        #endfor
        parents=list(set(tparents))
    #endwhile
    # 2. go down the subgraph
    while dist-1 > 0:
        dist-=1
        for pterm in pdist[dist]:
            if pterm in tValues.keys():
                pass
            else:
                parents=ont.parents(pterm)# 2.1 get parents
                tchildren , tn = gu.findAllChildrenInGraph(pterm,ont)# 2.3 find number of descendants in ontology graph for term pterm
                val=0
                for p in parents:# 2.4 calculate t-value
                    pchildren , pn = gu.findAllChildrenInGraph(p,ont)# 2.5 find number of descendants in ontology graph for term p
                    omega=tn/pn
                    val+=omega*tValues[p]
                #endfor
                tValues[pterm]=val/len(parents)# 2.6 tvalue is the mean of summmation
            #endif
        #endfor
    #endwhile
    # 3. do it for term t
    parents=ont.parents(t)# 3.1 get parents
    val=0
    tchildren , tn = gu.findAllChildrenInGraph(t,ont)# 3.2 find number of descendants in ontology graph for term t
    for p in parents:# 3.3 calculate t-value
        pchildren , pn = gu.findAllChildrenInGraph(p,ont)# 3.4 find number of descendants in ontology graph for term p
        omega=tn/pn
        val+=omega*tValues[p]
    #endfor
    if len(parents)==0:
        tValues[t]=1
    else:
        tValues[t]=val/len(parents)# 3.5 tvalue is the mean of summmation
#
#
#
def geneSemanticValue(g1 , g2 , gene1 , gene2 , ont):
    # 1. for each possible term pair get the maximum similarity
    sim = 0
    for term1 in gene1 :
        tsim = 0
        for term2 in gene2 :
            psv=ebm.pairSemanticValue(t1, t2, ont)
            if tsim < psv:
                tsim=psv
            #endif
        #endfor
        sim+=tsim
    #endfor
    for term2 in gene2 :
        tsim = 0
        for term1 in gene1 :
            psv=ebm.pairSemanticValue(t1, t2, ont)
            if tsim < psv:
                tsim=psv
            #endif
        #endfor
        sim+=tsim
    #endfor
    # use the formula to calculate gene similarity
    print(f'Similarity between genes {g1} , {g2} : {sim/(len(gene1) + len(gene2))}')
    return sim/(len(gene1) + len(gene2))
#
#
#
def averageLengthToLCA(t1 , t2 , lca , ont, G=None):
    if G==None :
        G=ont.get_graph()
    path1 = list(nx.all_simple_paths(G, lca, t1))
    path2 = list(nx.all_simple_paths(G, lca, t2))
    return sum([len(i) for i in path1 + path2])/len(path1+path2)
#
#
#
def lowest_common_ancestor(t1, t2,nt1 , nt2):
    '''
    1. Given two terms t1 , t2 find their LCA .
    '''
    for i in t1:
        for j in t2:
            intersection_1=list(set(t1[i])&set(t2[j]))
            # check if t1 or t2 is an ancestor
            intersection_2=list(set(t1[i])&set([nt1,nt2]))
            intersection_3=list(set(t2[j])&set([nt1,nt2]))
            if len(intersection_1)>0 :
                return [i+j , intersection_1]
            elif len(intersection_2)>0 :
                return [i+j , intersection_2]
            elif len(intersection_3)>0:
                return [i+j , intersection_3]
        #endfor
    #endfor
    return [np.nan,None]
#
#
#
def findMinimumPath(paths):
    dist = 10e10
    path = None
    for p in paths:
        if dist > len(p):
            dist=len(p)
            path=p
        #endif
    #endfor
    return path,dist
def findMinimumPathLCA(t1, t2 , lca , G):
    # 1. List all simple paths
    paths_1 = list(nx.all_simple_paths(G ,lca ,t1))
    paths_2 = list(nx.all_simple_paths(G ,lca ,t2))
    # 2. Find the minimum path
    path_1 , dist_1 = findMinimumPath(paths_1)
    path_2 , dist_2 = findMinimumPath(paths_2)
    path_1.pop(0)
    return path_1[::-1]+path_2 , dist_1+dist_2
#
#
#
def minimumPathLength(t1 , t2):
    if not(t1[1] == t2[1]):
        return 0
    dist_1 = 0
    dist_2 = 0
    # 1. for each layer of t1
    for i in t1[0]:
        # 2. for each layer of t2
        for j in t2[0]:
            intersection=list(set(t1[0][i])&set(t2[0][j]))
            if len(intersection)>0:
                dist_1 = i
                dist_2 = j
                break
            else:
                continue
        #endfor
    #endfor
    return dist_1 + dist_2
#
#
#
def pairSemanticValue(t1 , t2):
    semanticValue1 = sum([t1[t] for t in t1.keys()])
    semanticValue2 = sum([t2[t] for t in t2.keys()])
    # Calculate semantic similarity of t1 , t2
    inTerms = set(t1.keys()) & set(t2.keys())
    svalPair=0
    if not len(inTerms)==0:
        svalPair = sum([t1[i]+t2[i] for i in inTerms])
    #endif
    pairSimilarity = svalPair/(semanticValue1+semanticValue2)
    return pairSimilarity
#
#
#
def getSvalue(term , ont):
    # 1. calculate all semantic values for the graph
    sval = {}
    weightFactor = .815
    # 2. calculate S-Values
    sval[term]=1# 2.1 for term t S-Value calculation is trivial
    # 2.2 get parents of term
    # set last term to be the last child
    parents = ont.parents(term)
    lchildren = [term]
    while not(len(parents)==0):
        # 2.3 for each parent get its children
        newParents=[]
        for p in parents :
            newParents+=list(set(newParents)|set(ont.parents(p)))
            # 2.3.1 children must be in the subgraph
            valid_children = list(set(ont.children(p))&set(lchildren))
            maxSV=0
            for vc in valid_children:# 2.3.2 calculate max sValue
                if maxSV<weightFactor*sval[vc]:
                    maxSV=weightFactor*sval[vc]
                #endif
            #endfor
            # 2.3.3 set svalue for parent
            sval[p]=maxSV
        #endfor
        # 2.4 set parents last children
        # get the parents of parents
        del lchildren
        lchildren=parents
        del parents
        parents=list(set(newParents))
        print(len(parents))
        print('mlkeia')
    #endwhile
    return sval
#
#
#
def findLCAsAnc(t1, t2):
    # 1. for each of distances
    for d1 in sorted(list(t1.keys())):
        for d2 in sorted(list(t2.keys())):
            # 2. check if there is an intersection between ancestor set
            intersection = list(set(t1[d1])&set(t2[d2]))
            if len(intersection)>0:
                return intersection
            #endif
        #endfor
    #endfor
    return None# there is no LCA
#
#
#
def getAllParents(t, ont):
    ret = list()
    parents=list(ont.parents(t))
    while len(parents)>0:
        print(str(parents))
        tparents=[]
        for p in parents:
            tparents+=list(ont.parents(p))
        #endfor
        ret.extend(parents)
        parents=list(set(tparents))
    #endwhile
    return set(ret)
#
#
#
def findCommonParents(t1, t2 , ont):
    p1 = getAllParents(t1,ont)
    p2 = getAllParents(t2,ont)
    commonParents=list(p1&p2)
    if len(commonParents)==0:
        return None
    else:
        return commonParents
#
#
#
def findLCAs(t1, t2 , ont):
    # 1. get all parents and distances
    anc_1 = gu.allAncestors(t1, ont)
    anc_2 = gu.allAncestors(t2, ont)
    # 2. find lowest common ancestor
    return findLCAsAnc(anc_1 , anc_2)
#
#
#
def shortestSemanticDifferentiationDistance(geneData , ont):
    if os.path.exists(os.getcwd()+'/Datasets/ssddSimilarity.csv'):
        sim = pd.read_csv(os.getcwd()+'/Datasets/ssddSimilarity.csv')
        sim.drop(columns=['Unnamed: 0'] , inplace=True)
        sim.index=sim.columns
        return sim
    #endif
    start_time=time.time()
    # 1. extract all terms from gene data
    terms = pu.extractTerms(geneData)
    # 0. preprocessing steps
    G = ont.get_graph()
    ancestors = None
    if os.path.exists(os.getcwd()+'/Datasets/allAncestors.json'):
        ancestors = gu.read_json_file(os.getcwd()+'/Datasets/allAncestors.json')
    else:# create annotation file
        ancestors=gu.allAncestorsAllTerms(terms , ont)
        with open(os.getcwd()+'/Datasets/allAncestors.json', "w") as f:
            json.dump(ancestors, f, indent=4)
    #endif
    # 2. for each term calculate term's t-value
    tValues = {}
    counter=0
    print('Getting all T-Values !')
    if os.path.exists(os.getcwd()+'/Datasets/tValues.json'):
        tValues = gu.read_json_file(os.getcwd()+'/Datasets/tValues.json')
    else:
        for t in terms :
            counter+=1
            icu.progressBar(counter, len(terms), start_time)
            root ,namespace = gu.findRoot(t, ont)# 2.1 find root term
            getTValues(t, root , ont , tValues , G)# 2.2 calculate tvalues
        #endfor
        with open(os.getcwd()+'/Datasets/tValues.json', "w") as f:
            json.dump(tValues, f, indent=4)
    # 3. calculate term distance for each pair of terms in geneData
    ssddSim = pd.DataFrame(0, index=terms, columns=terms, dtype=np.float64)
    counter=0
    print('Calculating Shortest Semantic Differentation Distance !')
    for t1 in terms :
        for t2 in terms :
            counter+=1
            icu.progressBar(counter, len(terms)**2,start_time)
            if t1==t2 :# 3.1 trivial similarity
                ssddSim.loc[t1 , t2]=0
                continue
            #endif
            # 3.3 find LCAs
            lcas = findLCAsAnc(ancestors[t1], ancestors[t2])
            if lcas==None:# 3.3 trivial distance , must be 1 to turn into 0
                ssddSim.loc[t1 , t2]=1
                continue
            #endif
            # 3.4 find minimum path from lca to t1, t2
            minPathTvalue = 1
            for lca in lcas :
                mpath , n = findMinimumPathLCA(t1, t2 , lca , G)
                # 3.5 sum all tValues in path
                tvalSum = sum([tValues[node] for node in mpath])
                if tvalSum<minPathTvalue:
                    minPathTvalue=tvalSum
                #endif
            #endfor
            ssddSim.loc[t1 , t2]=minPathTvalue
        #endfor
    #endfor
    ssddSim = 1 - ssddSim
    print(ssddSim)
    ssddSim.to_csv(os.getcwd()+'/Datasets/ssddSimilarity.csv')
    return ssddSim
#
#
#
def semanticValueSimilarity(geneData , ont):
    if os.path.exists(os.getcwd()+'/Datasets/sValueSimilarity.csv'):
        sim = pd.read_csv(os.getcwd()+'/Datasets/sValueSimilarity.csv')
        sim.drop(columns=['Unnamed: 0'] , inplace=True)
        sim.index=sim.columns
        return sim
    #endif
    start_time=time.time()
    G = ont.get_graph()
    # 1. extract all terms from gene data
    terms = pu.extractTerms(geneData)
    # 2. for each subgraph calculate term's s-value
    sval = {}
    counter=0
    print('Getting s-values for every term')
    for t in terms:
        counter+=1
        icu.progressBar(counter, len(terms),start_time)
        sval[t] = getSvalue(t , ont)
    #endfor
    with open(os.getcwd()+'/Datasets/svalues.json', "w") as f:
        json.dump(sval, f, indent=4)
    # 4. for each term pair calculate their semantic similarity
    sValueSimilarity = pd.DataFrame(0, index=terms, columns=terms, dtype=np.float64)
    print('Calculate term similarity between all term pairs')
    counter=0
    for t1 in terms :
        for t2 in terms :
            counter+=1
            icu.progressBar(counter, len(terms)**2,start_time)
            if t1==t2 :# terms are the same term , trivial calculation of similarity
                sValueSimilarity.loc[t1 ,t2]=1
            else:
                sValueSimilarity.loc[t1 ,t2]=pairSemanticValue(sval[t1] , sval[t2])
            #endif
        #endfor
    #endfor
    print(sValueSimilarity)
    sValueSimilarity.to_csv(os.getcwd()+'/Datasets/sValueSimilarity.csv')
    return sValueSimilarity
#
#
#
def simpleWeightedDistance(geneData , ont):
    '''
    Input : a dictionary of gene names and corresponding go terms
            an ontology object from ontobio module
    - - - - -
    Process : for each term pair get their lowest common ancestor
              calculate the weighted distance of lca from root
              if they are not under the same subontology then their normalized distance is 1
              (they are semantically dissimilar).
    - - - - -
    Output : Weighted Distance
    '''
    if os.path.exists(os.getcwd()+'/Datasets/simSimpleWeights.csv'):
        sim = pd.read_csv(os.getcwd()+'/Datasets/simSimpleWeights.csv')
        sim.drop(columns=['Unnamed: 0'] , inplace=True)
        sim.index=sim.columns
        return sim
    #endif
    G = ont.get_graph()
    start_time=time.time()
    # 1. extract all terms from gene data
    terms = pu.extractTerms(geneData)
    # 2. for each term get all of its ancestors
    termDict = {}
    counter=0
    if os.path.exists(os.getcwd()+'/Datasets/allAncestors.json'):
        termDict = gu.read_json_file(os.getcwd()+'/Datasets/allAncestors.json')
    else:# create annotation file
        termDict=gu.allAncestorsAllTerms(terms , ont)
        with open(os.getcwd()+'/Datasets/allAncestors.json', "w") as f:
            json.dump(termDict, f, indent=4)
    #endif
    simSimpleWeights = pd.DataFrame(0, index=terms, columns=terms, dtype=np.float64)
    counter=0
    for t1 in terms :
        for t2 in terms:
            counter+=1
            icu.progressBar(counter, len(terms)**2,start_time)
            if t1==t2 :# trivial similarity t1,t2 are the same term
                simSimpleWeights.loc[t1 ,t2]=1
                continue
            #endif
            dist , lca = lowest_common_ancestor(termDict[t1.strip()], termDict[t2.strip()],t1 , t2)
            if not lca==None :
                try :
                    dist = termDict[lca[0]]
                    avgL = averageLengthToLCA(t1,t2,lca[0],ont,G)
                    simSimpleWeights.loc[t1 , t2] = (avgL + sum([.815,] + [.815**(i+1) for i in range(1,len(dist.keys()))]))/avgL
                except KeyError :
                    simSimpleWeights.loc[t1 , t2] = 0
            else:
                simSimpleWeights.loc[t1 , t2] = 0
            #endif
        #endfor
    #endfor
    # 4. save similarity
    print(simSimpleWeights)
    simSimpleWeights.to_csv(os.getcwd()+'/Datasets/simSimpleWeights.csv')
    return simSimpleWeights
#
#
#
def simRada(geneData, ont):
    '''
    Input : a dictionary of gene names and corresponding go terms
            an ontology object from ontobio module
    - - - - -
    Process : for each term pair get the minimum distance in graph
              if they are not under the same subontology then their normalized distance is 1
              (they are semantically dissimilar).
    - - - - -
    Output : Rada Distance
    '''
    if os.path.exists(os.getcwd()+'/Datasets/simRada.csv'):
        sim = pd.read_csv(os.getcwd()+'/Datasets/simRada.csv')
        sim.drop(columns=['Unnamed: 0'] , inplace=True)
        sim.index=sim.columns
        return sim
    #endif
    G = ont.get_graph()
    start_time = time.time()
    # 1. extract all terms from gene data
    terms = pu.extractTerms(geneData)
    # 2. create a dataframe full of 0s
    similarityRada = pd.DataFrame(0, index=terms, columns=terms, dtype=np.float64)
    counter=0
    for t1 in terms :
        for t2 in terms:
            counter+=1
            icu.progressBar(counter, len(terms)**2 , start_time)
            if t1==t2 :# trivial similarity t1,t2 are the same term
                similarityRada.loc[t1 ,t2]=1
                continue
            #endif
            try :
                spl = nx.shortest_path_length(G , t1, t2)
                similarityRada.loc[t1 ,t2]=spl
            except nx.exception.NetworkXNoPath :
                similarityRada.loc[t1 ,t2]=np.nan
        #endfor
    #endfor
    # 4. save similarity
    nrm = similarityRada.max().max()
    similarityRada.fillna(nrm*2,inplace=True)
    similarityRada=similarityRada/(nrm*2)
    print(similarityRada)
    similarityRada.to_csv(os.getcwd()+'/Datasets/simRada.csv')
    return similarityRada
#
#
#
