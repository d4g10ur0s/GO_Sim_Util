import sys
import json
# data analysis tools
import pandas as pd
# graph analysis tools
import networkx as nx
# ontology analysis tools
import ontobio as ob
# parallel programming
from dask.distributed import Client,LocalCluster,progress
import dask

def save_to_json(gene_dict, file_path):
    with open(file_path, 'w') as f:
        json.dump(gene_dict, f, indent=4)

def parseGAF(ga_path):
    '''
    Input : Path for .gaf file .
    - - - - -
    Process : Create a dictionary of DB:id as keys with values [(GO:id , evidence), ... , (GO:id, evidence)]
    - - - - -
    Output : dict of genes with their corresponding go terms and evidence for annotation
    '''
    print('Reading GAF .')
    print(str(ga_path))
    geneDict = {}
    gap = ob.io.gafparser.GafParser()
    with open(ga_path , 'r+') as gaf :
        line = gaf.readline()
        counter = 1
        while line :
            try : # first lines do not have useful information
                line = gaf.readline()
                gline = gap.parse_line(line)
                g = gline.associations
                for asc in g :
                    if not (asc.subject.id.namespace+':'+asc.subject.id.identity) in geneDict.keys():
                        geneDict[asc.subject.id.namespace+':'+asc.subject.id.identity] = []
                    geneDict[asc.subject.id.namespace+':'+asc.subject.id.identity].append((asc.object.id.namespace+asc.object.id.identity , str(asc.evidence.gaf_evidence_code())))
            except (AttributeError,TypeError) as e :
                print(str(e))
            finally :
                counter+=1
                if counter%10000 == 0:
                    print(f'Terms processed : {counter}')
        #endwhile
        return geneDict

'''
def retPath(ont ,G , root , n ):
    paths = list(nx.all_simple_paths(G, root, n))
    longest = len(max(paths, key=lambda p: len(p)))
    return {'id':n,
            'label': ont.label(n),
            'pathcount': len(paths),
            'longest': longest}

def get_pathstats(ont ,G , root , nodes):
    """
    for any given node, return a table row with stats
    """
    # parallel programming client
    #cluster = LocalCluster(processes=False, threads_per_worker=2)
    #client = Client(cluster)
    items = []
    counter=0
    for n in nodes:
        counter+=1
        items.append(retPath(ont ,G , root , n))
        if counter%10 == 0:
            print(f'Terms processed : {counter}')
    #progress(futures)
    return items
'''

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <obo_file_path> <ga_file_path>")
        sys.exit(1)
    # path for .obo
    obo_path = sys.argv[1]
    # path for .gaf
    ga_path = sys.argv[2]
    # create ontology
    #print('Creating Ontology from file .')
    #ont = ob.OntologyFactory().create(obo_path)
    geneDict = parseGAF(ga_path)
    save_to_json(geneDict,"./Datasets/genes.json")
    exit(0)
    '''
    # get root nodes
    print('Getting Root Nodes !')
    rootNodes=[]
    for r in ont.get_roots():
        rootNodes.append(r)
    # get ontology graph
    print('Creating Ontology Graph !')
    G = ont.get_graph()
    print(str(len(G)))
    # get different subontologies
    print('Descendants for : ' + ont.label(rootNodes[1]))
    cc = list(ont.descendants(rootNodes[1]))
    cc_df = None
    print('Descendants for : ' + ont.label(rootNodes[2]))
    bp = list(ont.descendants(rootNodes[2]))
    bp_df = None
    print('Descendants for : ' + ont.label(rootNodes[3]))
    mf = list(ont.descendants(rootNodes[3]))
    mf_df = None
    # get all obsolete terms , DO NOT USE THEM
    obs = ont.all_obsoletes()
    # create the corresponding dataframes
    print('Stats for CC !')
    cc_df = get_pathstats(ont , G , rootNodes[1] , cc)
    print(cc_df)
    print('Stats for BP !')
    bp_df = get_pathstats(ont , G , rootNodes[2] , bp)
    print(bp_df)
    print('Stats for MF !')
    mf_df = get_pathstats(ont , G , rootNodes[3] , mf)
    print(mf_df)
    '''

if __name__ == "__main__":
    main()
