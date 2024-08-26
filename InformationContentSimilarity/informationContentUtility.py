# basic modules
import sys
import os
import time
import random
# data analysis modules
import numpy as np
import pandas as pd
# ontology modules
import ontobio as ob
import networkx as nx
# custom modules
from Parsing import graphUtilities as gu
from Parsing import parsing as pu
from EdgeBasedSimilarity import edgeBasedMethods as ebm
#
#
#
def progressBar(count_value, total, suffix=''):
    emoji = [
      '👿','🤡','👽','💩','😤','🤬','🤯','🥶','🥵',
      '🤒','🤕','🤢','🤮','🤧','🥴','🤥','🤫','🤭',
      '🥸','🥲','🥱','🤫','🤭','😬','🤒','🤮','🤧',
      '🥴','🙁','😤','🙄',
      ]
    #emoji = ['\u1F615','\u1FAE4','\u1F61F','\u1F641','\u2639','\u1F62E','\u1F62F','\u1F632','\u1F633','\u1F97A','\u1F979','\u1F626','\u1F627','\u1F628','\u1F630','\u1F625','\u1F622','\u1F62D','\u1F631','\u1F616']
    #emoji = [e.decode('utf-8') for e in encoded_emojis]
    bar_length = 80
    filled_up_Length = int(round(bar_length* count_value / float(total)))
    percentage = round(100.0 * count_value/float(total),1)
    bstring = ''.join(random.choices(emoji, k=filled_up_Length))
    bar = bstring + '-' * (bar_length - filled_up_Length)
    sys.stdout.write('[%s] %s%s Nodes left : %s\n' %(bar, percentage, '%', total-count_value))
    #sys.stdout.flush()
#
#
#
def getAllParents(t , ont):
    tparents = ont.parents(t)
    parents = tparents
    while len(tparents)>0:
        temp = []
        for p in tparents :
            temp+=ont.parents(p)
        #endfor
        tparents=temp
        parents+=tparents
    #endwhile
    return parents
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
def parentFrequency(terms, ont):
    G=ont.get_graph()
    # 0. get all terms
    tKeys = list(terms.keys())
    # 1. for each term get all of its parents
    tcounter=0
    for t in tKeys :
        tcounter+=1
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
def calcParentFrequency(t , tf , ont):
    parents = ont.parents(t)
    # 1. for each immidiate parent add parent's frequency
    for p in parents:
        if p in tf.keys():
            tf[p] += tf[t]
        else:
            tf[p] = tf[t]
    #endfor
    # 2. do the same as long as parents exist
    nextKids=parents
    while len(nextKids) > 0:
        tparents=set([])
        #tcounter=0
        for p in nextKids:
            #tcounter+=1
            #progressBar(tcounter, len(nextKids))
            # 3. parent becomes the kid
            temp = ont.parents(p)
            tparents = tparents | set(temp)
            # 4. for each parent of parent add parent's frequency
            for pp in temp:
                if pp in tf.keys():
                    tf[pp] += tf[p]
                else:
                    tf[pp] = tf[p]
            #endfor
        #endfor
        nextKids = list(tparents)
    #endwhile
#
#
#
def calcIC(annoTerms , ont):
    G = ont.get_graph()
    tf = {}
    # 0. turn terms to dictionary
    for i in annoTerms['terms'].values.tolist():
        tf[i]=annoTerms[annoTerms['terms']==i]['frequency'].values[0]
    #endfor
    # 1. for each term get his parents
    parents = []
    for t in annoTerms['terms'].values.tolist() :
        tp = ont.parents(t)
        calcParentFrequency(t , tf , ont)
    #endfor
    df = pd.DataFrame.from_dict(tf, orient='index',)
    new_columns = ['frequency']
    df.columns = new_columns
    df['probability'] = df['frequency']/df['frequency'].max()
    df['IC']= - np.log(df['probability'])
    return df
#
#
#
def termFrequency(geneData , ancestors , ont ):
    tFrequency = {}
    # 1. Get term frequency of leaf terms from annotations
    for g in geneData :
        for t in geneData[g] :
            if t[0] in tFrequency.keys():
                tFrequency[t[0]]+=1
            else:
                tFrequency[t[0]]=1
            #endif
        #endfor
    #endfor
    # 2. Get ancestors' frequencies
    counter=0
    for t in list(ancestors.keys()):# 2.1 for each term
        counter+=1
        progressBar(counter, len(ancestors.keys()))
        tanc = ancestors[t]
        for dist in tanc:# 2.3 for each ancestor list
            for a in tanc[dist]:# 2.4 for each ancestor in ancestor list
                achildren=set(ont.children(a))# 2.5 get his children
                validChildren=set(tFrequency.keys())&achildren# 2.6 find intersection with dataset's terms valid children are already stored in termFrequency
                tFrequency[a]=sum([tFrequency[kid] for kid in validChildren])# 2.7 calculate frequency using children's frequency
            #endfor
        #endfor
    #endfor
    df = pd.DataFrame.from_dict(tFrequency, orient='index',)
    new_columns = ['frequency']
    df.columns = new_columns
    return df
#
#
#
def l_MICA(t1 , t2 , ont , ic):
    # 1. find all parents
    parents_1 = getAllParents(t1 , ont)
    parents_2 = getAllParents(t2 , ont)
    # 2. find all common parents
    commonParents = set(parents_1)&set(parents_2)
    if len(commonParents)==0:
        return 0
    commonParentsIC = ic[ic['terms'].isin(commonParents)]['IC']
    return commonParentsIC.max()
#
#
#
def simResnikMICA(t1, t2 , ont , df, G=None):
    if t1==t2:
        return 1
    # 0. get roots
    '''
    root1 ,namespace1 = gu.findRoot(t1, ont)
    root2 ,namespace2 = gu.findRoot(t2, ont)
    if not root1==root2 :
        print('No common root .')
        return 0
    '''
    mica = l_MICA(t1 , t2 , ont , df)
    return mica
#
#
#
def diShin(t1 , t2 , prob , ont, G=None):
    if G==None:
        G=ont.get_graph()
    #endif
    # 1. find all common ancestors
    anc_1 = getAllParents(t1, ont)
    anc_2 = getAllParents(t2, ont)
    commonParents = (set(anc_1)&set(anc_2))
    if len(commonParents)==0 :
        return 0 , {}
    #endif
    # 2. for each ancestor in common ancestors find the number of paths from ancestor to terms
    ancPath = {}
    for par in commonParents:
        n1 = len(list(nx.all_simple_paths(G , par , t1)))
        n2 = len(list(nx.all_simple_paths(G , par , t2)))
        diff = abs(n1-n2)
        if diff in ancPath.keys():
            ancPath[diff].append(par)
        else:
            ancPath[diff]=[par]
        #endif
    #endfor
    # 3. find the most informative from each difference
    simDiShin=0
    disjunctionPath = {}
    for diff in ancPath.keys():
        pic = None
        mica = 0
        for par in ancPath[diff]:
            tic = -np.log(prob[prob['terms']==par]['probability']).values[0]
            if mica<=tic:
                mica=tic
                pic=par
            #endif
        #endfor
        disjunctionPath[diff]=pic
        simDiShin+=mica
    #endfor
    simDiShin=simDiShin/len(ancPath.keys())# similarity is the average of disjunctive ICs
    return simDiShin ,disjunctionPath
#
#
#
def calculateDiShin(prob , ont):
    G=ont.get_graph()
    dishinSimilarity = pd.DataFrame(0, index=prob['terms'].values.tolist(), columns=prob['terms'].values.tolist(), dtype=np.float64)
    counter=0
    for t1 in dishinSimilarity.index :
        for t2 in dishinSimilarity.columns:
            counter+=1
            progressBar(counter , len(dishinSimilarity.index)**2)
            dishinSimilarity.loc[t1, t2], ancPath = diShin(t1 , t2 , prob , ont, G)
        #endfor
    #endfor
    dishinSimilarity.to_csv(os.getcwd()+'/Datasets/dishinSimilarity.csv')
    print(dishinSimilarity)
    return dishinSimilarity
#
#
#
def simGic(t1 , t2 , prob , ont):
    # 1. get all parents for t1
    anc_1 = getAllParents(t1, ont)
    anc_2 = getAllParents(t2, ont)
    commonParents = (set(anc_1)&set(anc_2))|set([t1 , t2])
    allParents = set(anc_1)|set(anc_2)|set([t1 , t2])
    # 2. calculate similarity based on all and common parents
    sum_1 = sum(prob[prob['terms'].isin(commonParents)]['probability'].values.tolist())
    sum_2 = sum(prob[prob['terms'].isin(allParents)]['probability'].values.tolist())
    graphSim = sum_1/sum_2
    return graphSim
#
#
#
def calculateSimGIC(prob , ont):
    graphICSimilarity = pd.DataFrame(0, index=prob['terms'].values.tolist(), columns=prob['terms'].values.tolist(), dtype=np.float64)
    counter=0
    for t1 in graphICSimilarity.index :
        for t2 in graphICSimilarity.columns:
            counter+=1
            progressBar(counter , len(graphICSimilarity.index)**2)
            graphICSimilarity.loc[t1, t2]=simGic(t1 , t2 , prob , ont)
        #endfor
    #endfor
    graphICSimilarity.to_csv(os.getcwd()+'/Datasets/graphICSimilarity.csv')
    print(graphICSimilarity)
    return graphICSimilarity
#
#
#
def simResnik(t1, t2 , prob , ont):
    if t1==t2:#trivial similarity
        return 1
    #endif
    # 1. find lcas
    lcas = ebm.findLCAs(t1, t2, ont)
    # 2. calculate similarity using lca
    # 2.1 if no lca , then 0 similarity
    if lcas==None:
        return 0
    else:
        mica = 0
        for p in lcas :
            ic = -np.log(prob[prob['terms']==p]['probability']).values[0]
            if ic > mica:
                mica = ic
            #endif
        #endfor
        return mica
    #endif
#
#
#
def calculateSimResnik(prob , ont):
    if os.path.exists(os.getcwd()+'/Datasets/resnikSimilarity.csv'):
        resnikSimilarity = pd.read_csv(os.getcwd()+'/Datasets/resnikSimilarity.csv')
        resnikSimilarity.index=resnikSimilarity.columns
        return resnikSimilarity
    #endif
    resnikSimilarity = pd.DataFrame(0, index=prob['terms'].values.tolist(), columns=prob['terms'].values.tolist(), dtype=np.float64)
    counter=0
    for t1 in resnikSimilarity.index :
        for t2 in resnikSimilarity.columns:
            counter+=1
            progressBar(counter , len(resnikSimilarity.index)**2)
            resnikSimilarity.loc[t1, t2]=simResnik(t1 , t2 , prob , ont)
        #endfor
    #endfor
    resnikSimilarity.to_csv(os.getcwd()+'/Datasets/resnikSimilarity.csv')
    print(resnikSimilarity)
    return resnikSimilarity
#
# Similarities using Resnik Measure
#
def similarityJiang(simRes , t1 , t2 , prob):
    ic1 = -np.log(prob[prob['terms']==t1]['probability']).values[0]
    ic2 = -np.log(prob[prob['terms']==t2]['probability']).values[0]
    simLin = -(2 * simRes)+ic1+ic2
    return simLin
#
#
#
def calculateSimJiang(prob , ont):
    # 1. calculate resnik similarity
    simRes = None
    if os.path.exists(os.getcwd()+'/Datasets/resnikSimilarity.csv'):
        simRes = pd.read_csv(os.getcwd()+'/Datasets/resnikSimilarity.csv')
        simRes.drop(columns=['Unnamed: 0'],inplace=True)
    else:
        # 1.1 calculate resnik similarity
        simRes = calculateSimResnik(prob , ont)
    #endif
    simRes.index=simRes.columns# 1.2 fix index
    print(simRes)
    jiangSimilarity = pd.DataFrame(0, index=prob['terms'].values.tolist(), columns=prob['terms'].values.tolist(), dtype=np.float64)
    counter=0
    for t1 in jiangSimilarity.index :
        for t2 in jiangSimilarity.columns:
            counter+=1
            progressBar(counter , len(jiangSimilarity.index)**2)
            jiangSimilarity.loc[t1, t2]=similarityJiang(simRes.loc[t1, t2] , t1 , t2 , prob)
        #endfor
        print(jiangSimilarity)
    #endfor
    jiangSimilarity.to_csv(os.getcwd()+'/Datasets/jiangSimilarity.csv')
    return jiangSimilarity
#
#
#
def similarityLin(simRes , t1 , t2 , prob):
    ic1 = -np.log(prob[prob['terms']==t1]['probability']).values[0]
    ic2 = -np.log(prob[prob['terms']==t2]['probability']).values[0]
    simJiang = (2 * simRes) / (ic1 + ic2)
    return simJiang
#
#
#
def calculateSimLin(prob , ont):
    # 1. calculate resnik similarity
    simRes = None
    if os.path.exists(os.getcwd()+'/Datasets/resnikSimilarity.csv'):
        simRes = pd.read_csv(os.getcwd()+'/Datasets/resnikSimilarity.csv')
        simRes.drop(columns=['Unnamed: 0'],inplace=True)
    else:
        # 1.1 calculate resnik similarity
        simRes = calculateSimResnik(prob , ont)
    #endif
    simRes.index=simRes.columns# 1.2 fix index
    linSimilarity = pd.DataFrame(0, index=prob['terms'].values.tolist(), columns=prob['terms'].values.tolist(), dtype=np.float64)
    counter=0
    for t1 in linSimilarity.index :
        for t2 in linSimilarity.columns:
            counter+=1
            progressBar(counter , len(linSimilarity.index)**2)
            linSimilarity.loc[t1, t2]=similarityLin(simRes.loc[t1, t2] , t1 , t2 , prob)
        #endfor
    #endfor
    linSimilarity.to_csv(os.getcwd()+'/Datasets/linSimilarity.csv')
    print(linSimilarity)
    return linSimilarity
#
#
#
def similarityRelevance(t1 , t2 , prob , ont):
    if t1==t2:# trivial similarity
        return 1
    #endif
    # 1. find lcas of t1 , t2
    lcas = ebm.findLCAs(t1, t2, ont)
    if lcas==None :# trivial similarity
        return 0
    #endif
    # 2. find maximum ic
    maxIC = 0
    maxLCA = None
    for lca in lcas:
        lcaIC = -np.log(prob[prob['terms']==lca]['probability']).values[0]
        if maxIC<=lcaIC:
            maxIC=lcaIC
            maxLCA=lca
        #endif
    #endfor
    ic1 = -np.log(prob[prob['terms']==t1]['probability']).values[0]
    ic2 = -np.log(prob[prob['terms']==t2]['probability']).values[0]
    simRel = ( (2 * maxIC ) / (ic1+ic2) ) * (1-prob[prob['terms']==lca]['probability']).values[0]
    return simRel
#
#
#
def calculateSimRel(prob , ont):
    relSimilarity = pd.DataFrame(0, index=prob['terms'].values.tolist(), columns=prob['terms'].values.tolist(), dtype=np.float64)
    counter=0
    for t1 in relSimilarity.index :
        for t2 in relSimilarity.columns:
            counter+=1
            progressBar(counter , len(relSimilarity.index)**2)
            relSimilarity.loc[t1, t2]=similarityRelevance(t1 , t2 , prob, ont)
        #endfor
        print(relSimilarity)
    #endfor
    relSimilarity.to_csv(os.getcwd()+'/Datasets/relSimilarity.csv')
    print(relSimilarity)
    return relSimilarity
#
#
#
def similarityIC(t1 , t2 , prob , ont):# its MICA not LCA
    if t1==t2:# trivial similarity
        return 1
    #endif
    # 1. find lcas of t1 , t2
    lcas = ebm.findCommonParents(t1, t2, ont)
    if lcas==None :# trivial similarity
        return 0
    #endif
    # 2. find maximum ic
    maxIC = 0
    maxLCA = None
    for lca in lcas:
        lcaIC = -np.log(prob[prob['terms']==lca]['probability']).values[0]
        if maxIC<=lcaIC:
            maxIC=lcaIC
            maxLCA=lca
        #endif
    #endfor
    ic1 = -np.log(prob[prob['terms']==t1]['probability']).values[0]
    ic2 = -np.log(prob[prob['terms']==t2]['probability']).values[0]
    pic = -np.log(prob[prob['terms']==maxLCA]['probability']).values[0]
    simIC =( (2 * pic) / (ic1+ic2) ) * ( 1 - (1/(1+pic)) )
    return simIC
#
#
#
def calculateSimInfoCoeff(prob , ont):
    simIC = pd.DataFrame(0, index=prob['terms'].values.tolist(), columns=prob['terms'].values.tolist(), dtype=np.float64)
    counter=0
    for t1 in simIC.index :
        for t2 in simIC.columns:
            counter+=1
            progressBar(counter , len(simIC.index)**2)
            simIC.loc[t1, t2]=similarityIC(t1 , t2 , prob , ont)
        #endfor
        print(simIC)
    #endfor
    print(simIC)
    simIC.to_csv(os.getcwd()+'/Datasets/simIC.csv')
    return simIC
#
#
#
def calculateInformationContent(geneData , ont):
    # 1. extract all terms from gene data
    terms = pu.extractTermsFromGenes(geneData)
    # 2. find all ancestors for each term
    ancestors = None
    if os.path.exists(os.getcwd()+'/Datasets/allAncestors.json'):
        ancestors = gu.read_json_file(os.getcwd()+'/Datasets/allAncestors.json')
    else:# create annotation file
        ancestors=gu.allAncestorsAllTerms(terms , ont)
        with open(os.getcwd()+'/Datasets/allAncestors.json', "w") as f:
            json.dump(ancestors, f, indent=4)
    #endif
    # 3. create a set with all terms
    allTerms = []
    for t in terms :
        tanc = ancestors[t]
        for k in list(tanc.keys()):# 3.1 get all ancestors of t
            allTerms+=list(tanc[k])
        #endfor
    #endfor
    allTerms=list(set(allTerms))# 3.2 turn it into a set for all terms to be unique
    # 4. find the frequency using child terms
    tFrequency = None
    if os.path.exists(os.getcwd()+'/Datasets/termFrequency.csv'):
        tFrequency = pd.read_csv(os.getcwd()+'/Datasets/termFrequency.csv')
    else:# create annotation file
        tFrequency = termFrequency(geneData , ancestors , ont)# returns a pandas dataframe
        # 4.1 save term frequency
        tFrequency.to_csv(os.getcwd()+'/Datasets/termFrequency.csv')
    # 5. calculate propabilities based on sub ontology
    df_prob = None
    if os.path.exists(os.getcwd()+'/Datasets/termProbability.csv'):
        df_prob = pd.read_csv(os.getcwd()+'/Datasets/termProbability.csv')
    else:# create annotation file
        tProb = {}
        counter=0
        for t in tFrequency.index:
            counter+=1
            progressBar(counter, len(tFrequency.index))
            root , namespace = gu.findRoot(t,ont)
            tProb[t]=tFrequency.loc[t]/tFrequency.loc[root]
        #endfor
        df_prob = pd.DataFrame.from_dict(tProb, orient='index',)
        new_columns = ['probability']
        df_prob.columns = new_columns
        df_prob.to_csv(os.getcwd()+'/Datasets/termProbability.csv')
    #endif
    return df_prob
