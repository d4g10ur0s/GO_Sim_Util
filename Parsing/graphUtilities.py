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
#from MachineLearningMethods.anc2vec import train as builder
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
def getAllLeaves(t , ont):
    tchildren = ont.children(t)
    children = tchildren
    leaves=[]
    while len(tchildren)>0:
        temp = []
        for p in tchildren :
            tc=ont.children(p)
            if len(tc)==0:# its a leaf node , no children
                leaves.append(p)
            #endif
            temp+=tc
        #endfor
        tchildren=temp
        children+=tchildren
    #endwhile
    return leaves
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
    dist = 1
    parents=list(ont.parents(t))
    while len(parents)>0:
        dists[dist]=parents
        dist+=1
        tparents=[]
        for p in parents:
            tparents+=list(ont.parents(p))
        #endfor
        parents=list(set(tparents))
    #endwhile
    return dists
#
#
#
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
        pdist = allAncestors(t , ont)
        anc[t]=pdist
    #endfor
    return anc
#
#
#
def main():
    output_file_path = 'interactions.csv'
    retrieve_and_store_interactions(output_file_path)

    if len(sys.argv) != 2:
        print("Usage: python script.py <obo_file_path> <ga_file_path>")
        sys.exit(1)
    '''
    # path for .obo
    obo_path = sys.argv[1]
    embeds = builder.fit(obo_path, embedding_sz=200, batch_sz=64, num_epochs=100)
    print(embeds['GO:0005737'])
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
    '''

if __name__ == "__main__":
    main()
