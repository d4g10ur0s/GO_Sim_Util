import sys
import os
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
from HybridSimilarity import hybridSimilarityUtils as hsu
#
#
#
def findRoot(t,ont):
    # 0. get roots
    rootNodes={}
    namespaces = []
    for r in ont.get_roots():
        rootNodes[str(ont.node(r)['label'])] = r
        namespaces.append(str(ont.node(r)['label']))
    # 1. get node information
    node = ont.node(t)
    root = None
    namespace = None
    for i in node['meta']['basicPropertyValues']:
        if i['val'] in namespaces :
            namespace=i['val']
            root = rootNodes[i['val']]
        #endif
    #endfor
    return [root,namespace]
#
#
#
def findAllChildrenInGraph(t , ont):
    dists = {}
    n=0#number of children in ontology graph
    dist=1
    children = ont.children(t)
    while len(children) > 0 :
        dists[dist]=children
        n+=len(children)
        tempChildren = []
        for c in dists[dist]:# for each last child
            tempChildren+=ont.children(c)# get its children and add them to temp
        #endfor
        dist+=1
        children=tempChildren
    #endwhile
    return dists , n
#
#
#
def extractTermsFromGenes(genes):
    terms = []
    for g in genes:
        gene = genes[g]
        for tset in gene :
            terms.append(tset[0])
        #endfor
    #endfor
    return list(set(terms))
#
#
#
def read_json_file(file_path):
  """Reads a JSON file and returns its contents as a Python object.

  Args:
    file_path: The path to the JSON file.

  Returns:
    The parsed JSON data as a Python object.
  """

  with open(file_path, 'r') as f:
    data = json.load(f)
  return data
#
#
#
def allAncestors(t , ont):
    '''
    Input : id of a go term
            ontology graph
    - - - - -
    Process : 1. get parents
              2. add parents
              3. if parent root , stop
              4. if parent !root , step 1
    - - - - -
    Output : a dictionary with all ancestors and distance from them
    '''
    dists = {}
    # 1. get root term
    root , namespace = findRoot(t,ont)
    # 2.  get all paths from root to term
    G = ont.get_graph()
    paths = list(nx.all_simple_paths(G , root , t))
    # 3. for each path get uncommon nodes only
    a = set()
    rootDist = 10e10
    for p in paths:# 4.1 add minimum path from root to t
        if len(p)<rootDist:
            rootDist = len(p)
        #endif
        a = a|set(p)
    #endfor
    #print(f'Uncommon nodes : {a}')
    dists[rootDist] = [root]# 4.2 add minimum dist to root
    # 4. exclude t and root
    a.remove(t)
    a.remove(root)
    # 5. for each parent find the minimum path to t
    for anc in a :
        path , ancDist = ebm.findMinimumPath(anc , t , G)
        if ancDist in dists.keys():# 5.3 add dist with the correct form
            dists[ancDist].append(anc)
        else :
            dists[ancDist] = [anc]
        #endif
    #endfor
    #print(dists)
    anc = (dists , list(a) , namespace)
    return anc

def allAncestorsAllTerms(terms , ont):
    '''
    Input : terms , ontology
    - - - - - - - -
    Algorithm :
        1. for each term get all paths to every ancestor
        2. find the minimum path to every ancestor
        3. add the distance to term's dictionary
    - - - - - - - -
    Output :
    {
        termID :
                (
                    {distanceFromAncestor : ancestorList},
                    ancestorList ,
                    namespace
                )

    }
    '''
    anc = {}
    counter=0
    # 1. for each term
    for t in terms :
        #dists = allAncestors(t)
        dists = {}
        # 2. get root term
        root , namespace = findRoot(t,ont)
        # 3.  get all paths from root to term
        G = ont.get_graph()
        paths = list(nx.all_simple_paths(G , root , t))
        # 4. for each path get uncommon nodes only
        a = set()
        rootDist = 10e10
        for p in paths:# 4.1 add minimum path from root to t
            if len(p)<rootDist:
                rootDist = len(p)
            #endif
            a = a|set(p)
        #endfor
        dists[rootDist] = [root]# 4.2 add minimum dist to root
        # 5. exclude t and root
        a.remove(t)
        a.remove(root)
        # 6. for each parent find the minimum path to t
        for anc in a :
            paths = list(nx.all_simple_paths(G , anc , t))
            ancDist = 10e10
            for p in paths : # 6.1 for each path
                if len(p) < ancDist :# 6.2 find minimum
                    ancDist=len(p)
                #endif
            #endfor
            if ancDist in dists.keys():# 6.3 add dist with the correct form
                dists[ancDist].append(anc)
            else :
                dists[ancDist] = [anc]
            #endif
        #endfor
        anc[t] = (dists , list(a) , namespace)
    #endfor
    return anc

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <obo_file_path> <ga_file_path>")
        sys.exit(1)
    # path for .obo
    obo_path = sys.argv[1]
    # read json file
    json_file = './Datasets/genes.json'
    geneData = read_json_file(json_file)
    genes = list(geneData.keys())
    # extract terms from each gene
    terms = extractTermsFromGenes(geneData)
    # create ontology
    print('Creating Ontology from file .')
    ont = ob.OntologyFactory().create(obo_path)
    annoTerms = None
    # create a csv file with children info
    if not os.path.exists('/home/d4gl0s/diploma/Experiment/Datasets/annotationTerms.csv'):
        annoTerms = icu.frequencyANDprobability(geneData , ont)
        annoTerms.to_csv('/home/d4gl0s/diploma/Experiment/Datasets/annotationTerms.csv')
    #endif
    annoTerms = pd.read_csv('/home/d4gl0s/diploma/Experiment/Datasets/annotationTerms.csv')
    annoTerms.columns=['terms','frequency']
    ic=None
    if os.path.exists('/home/d4gl0s/diploma/Experiment/Datasets/informationContent.csv'):
        ic = pd.read_csv('/home/d4gl0s/diploma/Experiment/Datasets/informationContent.csv')
    else:
        ic = icu.calcIC(annoTerms , ont)
        ic.to_csv('/home/d4gl0s/diploma/Experiment/Datasets/informationContent.csv')
    ic = pd.read_csv('/home/d4gl0s/diploma/Experiment/Datasets/informationContent.csv')
    ic.columns=['terms','frequency','probability','IC']
    P = None
    if os.path.exists('/home/d4gl0s/diploma/Experiment/Datasets/transitionProb.csv'):
        P = pd.read_csv('/home/d4gl0s/diploma/Experiment/Datasets/transitionProb.csv')
        P.set_index(ic['terms'].values.tolist(), inplace=True)
    else:
        # 1. get transition prob for each parent-child
        tp = hsu.getTransitionProb(ic, ont)
        # 2. initialize matrix P
        P = pd.DataFrame(0, index=ic['terms'].values.tolist(), columns=ic['terms'].values.tolist(), dtype=np.float64)
        print(P)
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
    print(P)
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
    print(W)
    print(W.loc[ic['terms'].tolist()][list(terms)].max())
    # 5. get leaves - parents matrix
    W_ll = W.loc[ic['terms'].tolist()][list(terms)]
    # 6. Get leaves hosted similarity matrix
    HSM = pd.DataFrame(0, index=ic['terms'].values.tolist(), columns=ic['terms'].values.tolist(), dtype=np.float64)
    tcounter=0
    for i in terms :
        for j in terms :
            tcounter+=1
            icu.progressBar(tcounter , len(terms)**2)
            anc , HSM.loc[i][j] = icu.simResnik(i,j,ont,ic)
        #endfor
    #endfor
    '''
        ic =
        ic = pd.concat([ic, -np.log(ic['probability'])], axis=1)
        new_columns = ['frequency', 'probability', 'IC']
        ic.columns = new_columns
        ic.to_csv('/home/d4gl0s/diploma/Experiment/Datasets/informationContent.csv')
    ic.columns = ['terms','frequency', 'probability', 'IC']
    print(f'{ic}')
    hsu.getTransitionProb(ic, ont)
    for i in range(250):
        t1 = random.choice(terms)
        t2 = random.choice(terms)
        icu.grasm(t1, t2 , ont , ic)
    # resnik
        anc , simRes = icu.simResnik(t1 , t2 , ont ,ic)
        if not anc==None:
            icu.similarityJiang(simRes , t1 , t2 , ic , anc)
            icu.similarityLin(simRes , t1 , t2 , ic , anc)
            icu.similarityRelevance(t1 , t2 , ic , anc)
            icu.similarityIC(t1 , t2 , ic , anc)
        else:
            print('*'*25)
            print(f'Term {t1} , {t2} similarity is : {0}')
            print('No common ancestor .')
            print('*'*25)
    '''

if __name__ == "__main__":
    main()
