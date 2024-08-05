import sys
# data analysis tools
import pandas as pd
# graph analysis tools
import networkx as nx
# ontology analysis tools
import ontobio as ob
# parallel programming
from dask.distributed import Client,LocalCluster
import dask

def retPath(ont ,G , root , n,counter):
    counter+=1
    paths = list(nx.all_simple_paths(G, root, n))
    longest = len(max(paths, key=lambda p: len(p)))
    # some info
    if counter % 100 == 0:
        print(f"Processed {counter} terms...")
    return {'id':n,
            'label': ont.label(n),
            'pathcount': len(paths),
            'longest': longest}

def get_pathstats(ont ,G , root , nodes):
    """
    for any given node, return a table row with stats
    """
    # parallel programming client
    cluster = LocalCluster(processes=False, threads_per_worker=8)
    client = Client(cluster)
    items = []
    counter=0
    futures = [client.submit(retPath,ont ,G , root , n,counter) for n in nodes]
    progress(futures)
    return futures

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <obo_file_path>")
        sys.exit(1)
    # path for .obo
    obo_path = sys.argv[1]
    # create ontology
    print('Creating Ontology from file .')
    ont = ob.OntologyFactory().create(obo_path)
    # get root nodes
    print('Getting Root Nodes !')
    rootNodes=[]
    for r in ont.get_roots():
        rootNodes.append(r)
    # get ontology graph
    print('Creating Ontology Graph !')
    G = ont.get_graph()
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

if __name__ == "__main__":
    main()
