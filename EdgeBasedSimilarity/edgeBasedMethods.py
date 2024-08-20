import os
import sys
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
def getTValues(t , root , ont):
    # 1. get sub graph with all parents
    anc = gu.allAncestors(t,root,ont)
    # 2. for each parent calculate T-value
    # 2.1 for each parent get all children till leaf nodes
    nChildren = {}
    for i in anc :
        for p in anc[i]:
            children , n = gu.findAllChildrenInGraph(p,ont)
            nChildren[p]=n
        #endfor
    #endfor
    # 2.2 calculate T-values
    tValues = {}
    for i in list(anc.keys())[::-1]:# start from root
        for p in anc[i]:
            if p == root:
                tValues[root]=1
            else :# if not root
                pp = ont.parents(p)# get parents of parent
                omega = [tValues[o]*(nChildren[p]/nChildren[o]) for o in pp]# calculate omega using kids
                # t-value is the average of parents' t-values multiplied by omega
                tValues[p]=sum(omega)/len(omega)
                #endfor
            #endif
        #endfor
    #endfor
    # Do the same for t
    # Add children for t
    children , n = gu.findAllChildrenInGraph(t,ont)
    nChildren[t]=n
    # Calculate t's T-value
    pp = ont.parents(t)# get parents of t
    omega = [tValues[o]*(nChildren[t]/nChildren[o]) for o in pp]# calculate omega using kids
    tValues[t]=sum(omega)/len(omega)
    return tValues
#
#
#
def shortestSemanticDifferentiationDistance(t1 , t2 , ont):
    # 1. find roots
    root1 ,namespace1 = gu.findRoot(t1, ont)
    root2 ,namespace2 = gu.findRoot(t2, ont)
    # NO COMMON ROOT THE SIMILARITY IS 0
    if not(root1==root2):
        return 0
    #endif
    # 2. calculate T - Values based on root
    tval1 = getTValues(t1, root1 , ont)
    tval2 = getTValues(t2, root2 , ont)
    # 3. get distances to all ancestors
    anc1 = gu.allAncestors(t1 , ont)
    anc2 = gu.allAncestors(t2 , ont)
    # 4. find the minimum path
    mpath = findMinimumPath(t1, t2, (anc1 , root1) , (anc2 , root2) ,ont)
    sim = 0
    if not mpath==None :# 5. if common ancestor exists
        for node in mpath : # 6. add up ancestor's T-value
            if node in tval1.keys():
                sim+=tval1[node]
            else:
                sim+=tval2[node]
        #endfor
        return 1 - np.arctan(sim)/(np.pi/2) # 7. return semantic similarity
    else :# 7. There is no common ancestor
        return 0
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
def findMinimumPath(t1, t2 ,G):
    # 1. List all simple paths
    paths = list(nx.all_simple_paths(G ,t1 ,t2))
    # 2. Find the minimum path
    dist = 10e10
    path = None
    for p in paths :
        if dist > len(p):
            dist=len(p)
            path=p
        #endif
    #endfor
    return path,dist
#
#
#
def minimumPathLength(t1 , t2):
    '''
    1.
    '''
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
    semanticValue2 = sum([t1[t] for t in t2.keys()])
    # Calculate semantic similarity of t1 , t2
    inTerms = set(t1.keys()) & set(t2.keys())
    svalPair = sum([t1[i]+t2[i] for i in inTerms])
    pairSimilarity = svalPair/(semanticValue1+semanticValue2)
    return pairSimilarity
#
#
#
def getSvalue(term , tAnc , ont):
    subgraph = tAnc[0]
    # 1. calculate all semantic values for the graph
    sval = {}
    weightFactor = .815
    sval[t]=1# 2.1 for term t S-Value calculation is trivial
    # 2. calculate S-Values
    for i in subgraph:# 2.2 for parent terms , S-Value depends on its children in subgraph
        if i==1 :# 2.2.1 for layer 1 terms , S-Value calculation is trivial
            for p in subgraph[i]:
                sval[t]=1*weightFactor
            #endfor
            continue
        #endif
        # 2.2.2 for each child term , find the intersection with subgraph's nodes from previous layer
        for p in subgraph[i]:
            sval[p]=0
            childrenIntersection = set(subgraph[i-1])&set(ont.children(p))
            for c in list(childrenIntersection):
                if sval[p]<sval[c]*weightFactor:# find the maximum
                    sval[p]=sval[c]*weightFactor
                #endif
            #endfor
        #endfor
    #endfor
    return sval
#
#
#
def semanticValueSimilarity(geneData , ont):
    G = ont.get_graph()
    # 1. extract all terms from gene data
    terms = pu.extractTermsFromGenes(geneData)
    # 2. for each term get its subgraph
    termDict = {}
    counter=0
    print('Getting all ancestors for every term')
    for t in terms :
        counter+=1
        icu.progressBar(counter, len(terms))
        anc = gu.allAncestors(t , ont , G)
        termDict[t] = anc
    #endfor
    # 3. for each subgraph calculate term's s-value
    sval = {}
    counter=0
    print('Getting s-values for every term')
    for t in termDict.keys():
        counter+=1
        icu.progressBar(counter, len(terms))
        sval[t] = getSvalue(t , termDict[t] , ont)
    #endfor
    # 4. for each term pair calculate their semantic similarity
    sValueSimilarity = pd.DataFrame(0, index=terms, columns=terms, dtype=np.float64)
    print('Calculate term similarity between all term pairs')
    for t2 in terms :
        for t2 in terms :
            if t1==t2 :# terms are the same term , trivial calculation of similarity
                sValueSimilarity.loc[t1 ,t2]=1
            elif not termDict[t1][2]==termDict[t2][2]:# terms are in different subontologies , trivial calculation of similarity
                sValueSimilarity.loc[t1 ,t2]=0
            else:
                sValueSimilarity.loc[t1 ,t2]=pairSemanticValue(sval[t1] , sval[t2])
            #endif
        #endfor
    #endfor
    print(sValueSimilarity)
    sValueSimilarity.to_csv(os.getcwd()+'/Datasets/sValueSimilarity.csv')
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
#
#
#
