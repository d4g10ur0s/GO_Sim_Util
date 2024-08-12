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
    ic = icu.frequencyANDprobability(geneData , ont)
    ic = pd.concat([ic, -np.log(ic['probability'])], axis=1)
    new_columns = ['frequency', 'probability', 'IC']
    ic.columns = new_columns
    print(f'{ic}')
    # 0. get roots
    rootNodes={}
    namespace = []
    for r in ont.get_roots():
        rootNodes[str(ont.node(r)['label'])] = r
        namespace.append(str(ont.node(r)['label']))
    # 1. get all ancestors
    for i in range(250):
        t1 = random.choice(terms)
        root1 = findRoot(t1, ont , namespace, rootNodes)
        anc1 = allAncestors(t1 , root1[0] , ont)
        t2 = random.choice(terms)
        root2 = findRoot(t2, ont , namespace, rootNodes)
        anc2 = allAncestors(t2 , root2[0] , ont)
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
            if npaths in ca.keys():
                ca[npaths].append(abs(npaths1 - npaths2))
            else:
                ca[npaths] = [abs(npaths1 - npaths2)]
            #endif
        #endfor
        # 4. exclude multiple ancestors from each level
        dca = {}
        for i in ca :
            p = None
            ic = -1
            for a in ca[i]:
                if ic['IC'][a] > ic :
                    ic = ic['IC'][a]
                    p = a
            #endfor
            dca[i] = (p , ic)
        #endfor
        print(dca)
        sm = 0
        for i in dca :
            sm+=dca[i][1]
        if len(dca.keys())==0:
            print(f'DCA similarity of {t1} , {t2} is {sm}')
        else:
            sm /= len(dca.keys())
            print(f'DCA similarity of {t1} , {t2} is {sm}')

    '''
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
