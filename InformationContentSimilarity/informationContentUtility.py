# data analysis modules
import numpy as np
import pandas as pd
# ontology modules
import ontobio as ob
import networkx as nx
# custom modules
from Parsing import graphUtilities as gu
def grasm(t1, t2 , ont , ic):
    # 0. get roots
    rootNodes={}
    namespace = []
    for r in ont.get_roots():
        rootNodes[str(ont.node(r)['label'])] = r
        namespace.append(str(ont.node(r)['label']))
    # 1. find roots
    root1 = gu.findRoot(t1, ont , namespace, rootNodes)
    root2 = gu.findRoot(t2, ont , namespace, rootNodes)
    anc1 = gu.allAncestors(t1 , root1[0] , ont)
    anc2 = gu.allAncestors(t2 , root2[0] , ont)
    # 2. get all common ancestors
    ancestor1 = []
    for i in anc1 :
        ancestor1+=anc1[i]
    ancestor2 = []
    for i in anc2 :
        ancestor2+=anc2[i]
    commonAnc = set(ancestor1)&set(ancestor2)
    # 3. for each ancestor find the number of paths
    G = ont.get_graph()
    ca = {}
    for i in commonAnc:
        npaths1 = len(list(nx.all_simple_paths(G, i, t1)))
        npaths2 = len(list(nx.all_simple_paths(G, i, t2)))
        if abs(npaths1 - npaths2) in ca.keys():
            ca[abs(npaths1 - npaths2)].append(i)
        else:
            ca[abs(npaths1 - npaths2)] = [i]
        #endif
    #endfor
    # 4. exclude multiple ancestors from each level
    dca = {}
    for i in ca :
        p = None
        tempic = -1
        for a in ca[i]:
            if ic['IC'][a] > tempic :
                tempic = ic['IC'][a]
                p = a
        #endfor
        dca[i] = (p , tempic)
    #endfor
    sm = 0
    for i in dca :
        sm+=dca[i][1]
    if len(dca.keys())==0:
        print(f'DCA similarity of {t1} , {t2} is {sm}')
    else:
        sm /= len(dca.keys())
        print(f'DCA similarity of {t1} , {t2} is {sm}')
#
#
#
def parentFrequency(terms, ont):
    # 1. for each term get its parents' frequency
    tKeys = list(terms.keys())
    counter = 0
    for t in tKeys :
        parents = ont.parents(t)# start with the immidiate parents
        while len(parents) > 0:
            tparents = []
            # 2. for each parent add one to parent frequency
            for p in parents:
                if p in list(terms.keys()):
                    terms[p]+=1
                else:
                    terms[p]=1
                #endif
                # 3. get parents of parents
                tparents += ont.parents(p)# u need to have multiple times one parent , so not use of set
            #endfor
            parents = tparents
        #endwhile
    #endfor
    tf = pd.DataFrame.from_dict(terms,orient='index')
    return tf
#
#
#
def frequencyANDprobability(geneData , ont):
    termFrequency = {}
    for g in geneData :
        for t in geneData[g] :
            if t[0] in termFrequency.keys():
                termFrequency[t[0]]+=1
            else:
                termFrequency[t[0]]=1
            #endif
        #endfor
    #endfor
    df = parentFrequency(termFrequency, ont)
    df = pd.concat([df, df/df[0].max()], axis=1)
    new_columns = ['frequency', 'probability']
    df.columns = new_columns
    return df
#
#
#
def simResnik(t1, t2 , ont , df):
    # 0. get roots
    rootNodes={}
    namespace = []
    for r in ont.get_roots():
        rootNodes[str(ont.node(r)['label'])] = r
        namespace.append(str(ont.node(r)['label']))
    root1 ,namespace1 = gu.findRoot(t1, ont , namespace, rootNodes)
    root2 ,namespace2 = gu.findRoot(t2, ont , namespace, rootNodes)
    anc1 = gu.allAncestors(t1 , root1 , ont)
    anc2 = gu.allAncestors(t2 , root2 , ont)
    # find the common ancestor with the minimum IC
    ancset1 = []
    for i in anc1.keys():
        ancset1 += anc1[i]
    ancset2 = []
    for i in anc2.keys():
        ancset2 += anc2[i]
    ancIntersection = set(ancset1)&set(ancset2)
    if len(ancIntersection)<1:
        return [None, 0]
    else:
        pic = -1
        anc = None
        for p in ancIntersection:
            if df['IC'][p] > pic:
                pic=df['IC'][p]
                anc = p
            #endif
        #endfor
        print('*'*25)
        print(f'Term {t1} , {t2} Resnik similarity given by common ancestor {anc} is : {pic}')
        return [anc, pic]
    #endif
#
# Similarities using Resnik Measure
#
def similarityJiang(simRes , t1 , t2 , ic , anc):
    simJiang = (2 * simRes) / (ic['IC'][t1]+ic['IC'][t2])
    print('-'*25)
    print(f'Term {t1} , {t2} Jiang similarity given by common ancestor {anc} is : {simJiang}')
    print('-'*25)
def similarityLin(simRes , t1 , t2 , ic , anc):
    simLin = (2 * simRes) - ic['IC'][t1] -ic['IC'][t2]
    print(f'Term {t1} , {t2} Lin similarity given by common ancestor {anc} is : {simLin}')
    print('*'*25)
def similarityRelevance(t1 , t2 , ic , anc):
    simRel = ( (2 * np.log(ic['probability'][anc]) ) / (np.log(ic['probability'][t1])+np.log(ic['probability'][t2])) ) * (1-ic['probability'][anc])
    print(f'Term {t1} , {t2} Relevance similarity given by common ancestor {anc} is : {simRel}')
    print('*'*25)
def similarityIC(t1 , t2 , ic , anc):
    simIC = ( (2 * ic['probability'][anc]) / (np.log(ic['probability'][t1])+np.log(ic['probability'][t2])) ) * (1- ( 1/(1+ic['probability'][anc]) ) )
    print(f'Term {t1} , {t2} IC similarity given by common ancestor {anc} is : {simIC}')
    print('*'*25)
