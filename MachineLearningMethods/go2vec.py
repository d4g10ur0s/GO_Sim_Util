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
from node2vec import Node2Vec
# ontology analysis tools
import ontobio as ob
# parallel programming
from dask.distributed import Client,LocalCluster,progress
import dask
# custom modules
from EdgeBasedSimilarity import edgeBasedMethods as ebm
from InformationContentSimilarity import informationContentUtility as icu
from HybridSimilarity import hybridSimilarityUtils as hsu

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <obo_file_path> <ga_file_path>")
        sys.exit(1)
    # path for .obo
    obo_path = sys.argv[1]
    # create ontology
    print('Creating Ontology from file .')
    ont = ob.OntologyFactory().create(obo_path)
    G = ont.get_graph().to_undirected()# graph must be undirected
    node2vec = Node2Vec(G, dimensions=64, walk_length=30, num_walks=200, workers=2)
    print('Training')
    model = node2vec.fit(window=10, min_count=1, batch_words=10)
    model.wv.most_similar('GO:1990462')
#
#
#
if __name__ == "__main__":
    main()
