# data analysis modules
import numpy as np
import pandas as pd
# ontology modules
import ontobio as ob
import networkx as nx
# custom modules
from Parsing import graphUtilities as gu

def parentFrequency(terms, ont):
    pf = {}
    # 1. for each term get its parents' frequency
    tKeys = terms.columns.tolist()
    counter = 0
    for t in tKeys :
        print(f'Number of terms processed : {counter}')
        counter+=1
        parents = ont.parents(t)
        while len(parents) > 0:
            tparents = []
            # 2. for each parent add one to parent frequency
            for p in parents:
                if p in pf.keys():
                    pf[p]+=1
                else:
                    pf[p]=1
                #endif
                # 3. get parents of parents
                tparents += ont.parents(p)# u need to have multiple times one parent , so not use of set
            #endfor
            parents = tparents
        #endwhile
    #endfor
    tf = pd.DataFrame.from_dict(pf,orient='index')
    return tf
