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
    for t in tKeys :
        parents = ont.parents(t)
        while len(parents) > 0:
            # 2. for each parent add one to parent frequency
            for p in parents:
                if p in pf.keys():
                    pf[p]+=1
                else:
                    pf[p]=1
                #endif
            #endfor
            # 3. get parents of parents
            parents = ont.parents(t)
        #endwhile
    #endfor
    print(pf)
    tf = pd.DataFrame.from_dict(pf,orient='index')
    print(str(tf))
    return tf
