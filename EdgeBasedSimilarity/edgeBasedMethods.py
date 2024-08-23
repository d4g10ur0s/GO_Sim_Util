import os
import sys
import json
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
def lowest_common_ancestor(t1, t2):
    '''
    1. Given two terms t1 , t2 find their LCA .
    2. Find maximum distance of LCA from root .
    '''
    if not(t1[1] == t2[1]):
        return [0,None]
    dist_1 = 0
    dist_2 = 0
    for i in t1[0]:
        for j in t2[0]:
            intersection=list(set(t1[0][i])&set(t2[0][j]))
            if len(intersection)>0:
                return [i+j , intersection]
        #endfor
    #endfor
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
    sval[term]=1# 2.1 for term t S-Value calculation is trivial
    children = [term]
    parents = ont.parents(term)
    # 2. calculate S-Values
    while not(len(parents)==0):
        # 2.1 get parents of parents
        ochildren = []# children in subgraph
        oparents = []
        for p in parents :
            oparents+=list(ont.parents(p))# 2.2 get parents of p in subgraph
            tchildren = set(children)&set(ont.children(p))# 2.2 get children of p in subgraph
            maxSval=0
            for c in tchildren :# children are already in sval
                if sval[c]*weightFactor>maxSval*weightFactor:# 2.3 find maximum svalue
                    maxSval=sval[c]*weightFactor
                #endif
            #endfor
            sval[p]=maxSval# 2.4 set svalue
            ochildren+=list(tchildren)# update parents' children
        #endfor
        parents = list(set(oparents))# 2.5 set new parents
        children = list(set(ochildren))# 2.5 set new children
    #endwhile
    return sval
#
#
#
def findLCAs(t1, t2):
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
def shortestSemanticDifferentiationDistance(geneData , ont):
    # 1. extract all terms from gene data
    terms = pu.extractTermsFromGenes(geneData)
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
            icu.progressBar(counter, len(terms))
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
            icu.progressBar(counter, len(terms)**2)
            if t1==t2 :# 3.1 trivial similarity
                ssddSim.loc[t1 , t2]=0
                continue
            #endif
            # 3.3 find LCAs
            lcas = findLCAs(ancestors[t1], ancestors[t2])
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
    G = ont.get_graph()
    # 1. extract all terms from gene data
    terms = pu.extractTermsFromGenes(geneData)
    # 2. for each subgraph calculate term's s-value
    sval = {}
    counter=0
    print('Getting s-values for every term')
    for t in terms:
        counter+=1
        icu.progressBar(counter, len(terms))
        sval[t] = getSvalue(t , ont)
    #endfor
    with open(os.getcwd()+'/Datasets/svalues.json', "w") as f:
        json.dump(sval, f, indent=4)
    # 4. for each term pair calculate their semantic similarity
    sValueSimilarity = pd.DataFrame(0, index=terms, columns=terms, dtype=np.float64)
    print('Calculate term similarity between all term pairs')
    for t1 in terms :
        for t2 in terms :
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
    G = ont.get_graph()
    # 1. extract all terms from gene data
    terms = pu.extractTermsFromGenes(geneData)
    # 2. for each term get all of its ancestors
    termDict = {}
    counter=0
    for t in terms :
        counter+=1
        icu.progressBar(counter, len(terms))
        anc = gu.allAncestors(t , ont , G)
        termDict[t] = anc
    #endfor
    simSimpleWeights = pd.DataFrame(0, index=terms, columns=terms, dtype=np.float64)
    counter=0
    for t1 in terms :
        for t2 in terms:
            counter+=1
            icu.progressBar(counter, len(terms)**2)
            if t1==t2 :# trivial similarity t1,t2 are the same term
                simSimpleWeights.loc[t1 ,t2]=1
                continue
            #endif
            dist , lca = lowest_common_ancestor( (termDict[t1][0],termDict[t1][2]) , (termDict[t2][0],termDict[t2][2]) )
            if not lca==None :
                try :
                    dist , anc , namespace = gu.allAncestors(lca[0], ont)
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
    G = ont.get_graph()
    '''
    1. for every possible combination of nodes find the minimum path between them .
    2. sum all the distances and divide by product of number of elements of each set .
    '''
    # 1. extract all terms from gene data
    terms = pu.extractTermsFromGenes(geneData)
    # 2. for each term get all of its ancestors
    termDict = {}
    counter=0
    for t in terms :
        counter+=1
        icu.progressBar(counter, len(terms))
        anc = gu.allAncestors(t , ont , G)
        termDict[t] = anc
    #endfor
    # 3. create a dataframe full of 0s
    similarityRada = pd.DataFrame(0, index=terms, columns=terms, dtype=np.float64)
    counter=0
    for t1 in terms :
        for t2 in terms:
            counter+=1
            icu.progressBar(counter, len(terms)**2)
            if t1==t2 :# trivial similarity t1,t2 are the same term
                similarityRada.loc[t1 ,t2]=1
                continue
            #endif
            dist , lca = lowest_common_ancestor( (termDict[t1][0],termDict[t1][2]) , (termDict[t2][0],termDict[t2][2]) )
            similarityRada.loc[t1 , t2] = dist
        #endfor
    #endfor
    # 4. save similarity
    print(similarityRada)
    similarityRada.to_csv(os.getcwd()+'/Datasets/simRada.csv')
    return similarityRada
#
#
#
