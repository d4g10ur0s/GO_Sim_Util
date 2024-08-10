# data analysis modules
import numpy as np
import pandas as pd
# ontology modules
import ontobio as ob
import networkx as nx
# custom modules
from Parsing import graphUtilities as gu

def parentFrequency(terms, ont):
    pf = {}
    # 1. for each term get its parents' frequency
    tKeys = terms.columns.tolist()
    counter = 0
    for t in tKeys :
        print(f'Number of terms processed : {counter}')
        counter+=1
        parents = ont.parents(t)
        while len(parents) > 0:
            tparents = []
            # 2. for each parent add one to parent frequency
            for p in parents:
                if p in pf.keys():
                    pf[p]+=1
                else:
                    pf[p]=1
                #endif
                # 3. get parents of parents
                tparents += ont.parents(p)# u need to have multiple times one parent , so not use of set
            #endfor
            parents = tparents
        #endwhile
    #endfor
    tf = pd.DataFrame.from_dict(pf,orient='index')
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
    tf = pd.DataFrame.from_dict(termFrequency,orient='index')
    df = pd.concat([tf.transpose(), parentFrequency(tf.transpose() , ont).transpose()], axis=1)
    df = df.transpose()
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
        print(f'Term {t1} , {t2} similarity is : {0}')
        return [None, 0]
    else:
        pic = 10e10
        anc = None
        for p in ancIntersection:
            if df['IC'][p] < pic:
                pic=df['IC'][p]
                anc = p
            #endif
        #endfor
        return [anc, pic]
    #endif
