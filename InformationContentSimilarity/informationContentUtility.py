# system modules
import sys
import time
# data analysis modules
import numpy as np
import pandas as pd
# ontology modules
import ontobio as ob
import networkx as nx
# custom modules
from Parsing import graphUtilities as gu
def progressBar(count_value, total, suffix=''):
    bar_length = 120
    filled_up_Length = int(round(bar_length* count_value / float(total)))
    percentage = round(100.0 * count_value/float(total),1)
    bar = '=' * filled_up_Length + '>' + '-' * (bar_length - filled_up_Length)
    sys.stdout.write('[%s] %s%s Nodes left : %s\r' %(bar, percentage, '%', total-count_value))
    sys.stdout.flush()
#
#
#
def calcICT(terms, ont):
    # 0. get roots
    rootNodes={}
    namespace = []
    for r in ont.get_roots():
        rootNodes[str(ont.node(r)['label'])] = r
        namespace.append(str(ont.node(r)['label']))
    # 1. calculate ICT values
    # 1.1 get number of terms
    nTerms = ont.get_graph().number_of_nodes()
    # 1.2 for each term add its parents to term dataset
    for a in anc :
        for i in anc[a][0].keys():
            terms+=anc[a][0][i]
        #endfor
    #endfor
    # 1.3 for each term in the genes dataset find number of children
    ict = {}
    for t in set(terms) :
        dists , n = gu.findAllChildrenInGraph(t, ont)
        if n==0:
            ict[t] = 0
        else:
            ict[t] = -np.log(n/nTerms)
    #endfor
    # 1. get roots
    # 2. reconstruct each ontology graph
    # 2.1 get roots
    rootNodes={}
    namespace = []
    for r in ont.get_roots():
        rootNodes[str(ont.node(r)['label'])] = r
        namespace.append(str(ont.node(r)['label']))
    # 2.2 get all nodes from each ontology
    ontologies = {}
    for i in ict.keys() :
        root , subont = gu.findRoot(i,ont,namespace,rootNodes)
        if subont in ontologies.keys():
            ontologies[subont].append(i)
        else:
            ontologies[subont]=[root,i]
        #endif
    #endfor
    # 2.3 find cutoff to be least 20% of values
    cutoffValues = {}
    for sub in ontologies:
        avg = 0
        for t in ontologies[sub]:
            avg+=ict[t]
        #endfor
        avg/=len(ontologies[sub])
        avg=avg - avg * 0.8
        cutoffValues[sub]=avg
    #endfor
    print(cutoffValues)
    # 2.4 construct meta roots
    # 2.4.1 find dist 1 nodes with similar ict values
#
#
#
def grasm(t1, t2 , ont , ic):
    # 1. find all ancestors roots
    anc1 = gu.allAncestors(t1 , ont)
    anc2 = gu.allAncestors(t2 , ont)
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
    G=ont.get_graph()
    # 0. get all terms
    tKeys = list(terms.keys())
    # 1. for each term get all of its parents
    tcounter=0
    for t in tKeys :
        tcounter+=1
        print(f'Processing term : {t}\n')
        progressBar(tcounter, len(tKeys))
        # 2. get root
        rterm , namespace = gu.findRoot(t,ont)
        # 3. get all paths from root to term
        paths = list( nx.all_simple_paths(G , rterm , t) )
        # 4. for each path
        counter=0
        totalTerms=sum([len(p) for p in paths])
        for p in paths :
            # 5. for each parent term in path
            for pt in p :
                if pt == t :# DO NOT COUNT THE TERM t
                    pass
                elif pt in terms.keys():# add one more to pt
                    terms[pt]+=1
                else:
                    terms[pt]=1# first time seeing parent
                #endif
            #endfor
        #endfor
    #endfor
#
#
#
def frequencyANDprobability(geneData , ont):
    termFrequency = {}
    # 1. Get term frequency from annotations
    for g in geneData :
        for t in geneData[g] :
            if t[0] in termFrequency.keys():
                termFrequency[t[0]]+=1
            else:
                termFrequency[t[0]]=1
            #endif
        #endfor
    #endfor
    # 2. Get parent frequency from children term frequency
    parentFrequency(termFrequency , ont)
    df = pf.DataFrame.from_dict(data, orient='index',)
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
