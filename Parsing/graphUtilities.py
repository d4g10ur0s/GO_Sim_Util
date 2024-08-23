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
        children=list(set(tempChildren))
    #endwhile
    return dists , n
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
def allAncestors(t , ont , G=None):
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
    if G==None :
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
    print('Finding all ancestors for all terms .')
    anc = {}
    counter=0
    for t in terms :
        counter+=1
        icu.progressBar(counter, len(terms))
        pdist = {}
        dist = 1
        parents=list(ont.parents(t))
        while len(parents)>0:
            pdist[dist]=parents
            dist+=1
            tparents=[]
            for p in parents:
                tparents+=list(ont.parents(p))
            #endfor
            parents=list(set(tparents))
        #endwhile
        anc[t]=pdist
    #endfor
    return anc
#
#
#
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
    # 0. choose two random terms
    G = ont.get_graph()
    for i in range(250):
        icu.progressBar(i, 251)
        t1 = random.choice(terms)
        t2 = random.choice(terms)
        hsu.hybridRSS(t1, t2 , ont , ic , G)
    #endfor


if __name__ == "__main__":
    main()
