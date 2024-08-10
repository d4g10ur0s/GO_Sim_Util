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
#
#
#
def findRoot(t,ont ,namespaces,rootNodes):
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
    - - - - -
    Output : ancestors of all terms
    '''
    # 1. get roots
    print('Getting Root Nodes !')
    rootNodes={}
    namespace = []
    for r in ont.get_roots():
        rootNodes[str(ont.node(r)['label'])] = r
        namespace.append(str(ont.node(r)['label']))
    # 2. for each term get all ancestors
    anc = {}
    counter=0
    for t in terms :
        try :
            #print(str(t))
            node = ont.node(t)
            #print(f'Node\'s namespace : {node['meta']['basicPropertyValues'][0]['val']}')
            root = None
            tnamespace = None
            for i in node['meta']['basicPropertyValues']:
                if i['val'] in namespace :
                    root = rootNodes[i['val']]
                    tnamespace=i['val']
            anc[t] = (allAncestors(t , root , ont) , tnamespace)
            counter+=1
            if counter%10 == 0:
                print(f'Terms processed : {counter}')
            #print(str(anc[t]))
        except KeyError as e:
            print(str(e))
            print(str(ont.node(t)))
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
    anc = allAncestorsAllTerms(terms,ont)
    '''
    rada similarity
    simRada = ebm.simRada(genes , geneData , anc)
    print(str(simRada))
    '''
    # 1. get roots
    print('Getting Root Nodes !')
    rootNodes={}
    namespace = []
    for r in ont.get_roots():
        rootNodes[str(ont.node(r)['label'])] = r
        namespace.append(str(ont.node(r)['label']))
    # ---
    #
    # ---
    # 0. get two random nodes
    for i in range(100):
        t1 = random.choice(terms)
        t2 = random.choice(terms)
        root1 ,namespace1 = findRoot(t1, ont , namespace, rootNodes)
        root2 ,namespace2 = findRoot(t2, ont , namespace, rootNodes)
        tval1 = ebm.getTValues(t1, root1 , ont)
        tval2 = ebm.getTValues(t2, root2 , ont)
        anc1 = allAncestors(t1 , root1 , ont)
        anc2 = allAncestors(t2 , root2 , ont)
        mpath = ebm.findMinimumPath(t1, t2, (anc1 , root1) , (anc2 , root2) ,ont)
        sim = 0
        print(f'Shortest path : {mpath}')
        if not mpath==None :
            for node in mpath :
                if node in tval1.keys():
                    sim+=tval1[node]
                else:
                    sim+=tval2[node]
            #endfor
            print(f'Semantic Distance of {t1} , {t2} : {np.arctan(sim)/(np.pi/2)}')
            print(f'Semantic Similarity of {t1} , {t2} : {1 - np.arctan(sim)/(np.pi/2)}')
        else :
            print(f'Semantic Distance of {t1} , {t2} : {1}')
            print(f'Semantic Similarity of {t1} , {t2} : {0}')
    #d = findAllChildrenInGraph(tnew,ont)

if __name__ == "__main__":
    main()
