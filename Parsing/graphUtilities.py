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
def allAncestors(term , root , ont):
    '''
    Input : id of a go term
            subontology root
            ontology graph
    - - - - -
    Process : 1. get parents
              2. add parents
              3. if parent root , stop
              4. if parent !root , step 1
    - - - - -
    Output : a dictionary with all ancestors and distance from them
    '''
    ret = {}
    dist = 1
    parents = ont.parents(term)# 1. get parents
    ret[dist] = parents# 2. store parents using distance
    while not len(parents)==0:
        temp = []
        for p in parents :
            temp = temp + ont.parents(p)
        #endfor
        dist+=1
        # * unique values *
        parents = list(set(temp))
        # 2. store parents using distance
        ret[dist] = parents
    #endwhile
    return ret

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
    '''anc = allAncestorsAllTerms(terms,ont)'''
    ic = icu.frequencyANDprobability(geneData , ont)
    ic = pd.concat([ic, -np.log(ic['probability'])], axis=1)
    new_columns = ['frequency', 'probability', 'IC']
    ic.columns = new_columns
    print(f'{ic}')
    hsu.getTransitionProb(ic, ont)
    '''
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
