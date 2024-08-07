import sys
sys.path.append('/home/d4gl0s/diploma/Experiment')
import json
import random
# data analysis tools
import pandas as pd
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
    t = list(anc.keys())
    for i in range(50):
        t1 = random.choice(t)
        t2 = random.choice(t)
        lcaDist , lca = ebm.lowest_common_ancestor(anc[t1],anc[t2])
        if not lca==None :
            print('-'*20)
            print(f'LCA : {lca[0]}')
            # find lca's maximum path to root
            node = ont.node(lca[0])
            root = None
            tnamespace = None
            for i in node['meta']['basicPropertyValues']:
                if i['val'] in namespace :
                    root = rootNodes[i['val']]
                    tnamespace=i['val']
            dist = allAncestors(lca[0], root , ont)
            print(f'All ancestors : {str(dist)}')
            print(f'Distance from root : {sum([.815,] + [.815**(i+1) for i in range(1,len(dist.keys()))])}')
            avgL = ebm.averageLengthToLCA(t1,t2,lca[0],ont)
            print(f'Average path length : {avgL}')
            normalizationFactor = (avgL + sum([.815,] + [.815**(i+1) for i in range(1,len(dist.keys()))]))/avgL
            print('-'*20)

if __name__ == "__main__":
    main()
