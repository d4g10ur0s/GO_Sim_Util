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
    # 1. get a random node
    for k in range(100):
        t1 = random.choice(terms)
        t2 = random.choice(terms)
        sval1 , semanticValue1 = ebm.getSvalue(t1 , ont , rootNodes , namespace)
        sval2 , semanticValue2 = ebm.getSvalue(t2 , ont , rootNodes , namespace)
        # Calculate semantic similarity of t1 , t2
        inTerms = set(sval1.keys()) & set(sval2.keys())
        svalPair = sum([sval1[i]+sval2[i] for i in inTerms])
        print(f'Similarity of {t1} , {t2} : {svalPair/(semanticValue1+semanticValue2)}')
    '''
    # 2. construct its graph
    node = ont.node(t1)
    root = None
    for i in node['meta']['basicPropertyValues']:
        if i['val'] in namespace :
            root = rootNodes[i['val']]
    graph1 = allAncestors(t1,root,ont)
    # 3. calculate all semantic values for the graph
    sval = {}
    weightFactor = .815
    # starting node's semantic value is 1
    sval[t1] = 1
    for i in graph1 :
        lterms = graph1[i]
        # semantic value depends on children
        for p in lterms :
            # if first layer , then term's child is the term under process
            if i==1 :
                sval[p] = 1 * weightFactor
            elif len(lterms)>0 :
                # 3.1 get node's children in ontology graph
                children = ont.children(p)
                # 3.2 get the intersection with all the previous layers
                player = graph1[i-1]
                childIntersection = set(sval.keys())&set(children)
                print(f'Last layer terms : {str(lterms)} , {str(i-1)}')
                print(f'Children of {p} Intersection : {str(childIntersection)}')
                sval[p]=0
                for c in childIntersection:
                    # 3.3 find the maximum svalue
                    if sval[p]<sval[c] * weightFactor:
                        sval[p] = sval[c] * weightFactor
                    #endif
                #endfor
            #endif
        #endfor
    #endfor
    print(str(sval))
    print(f'Term {t1} Semantic Value : {str(sum([sval[i] for i in sval.keys()]))}')
    '''

if __name__ == "__main__":
    main()
