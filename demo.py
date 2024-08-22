import os
import sys
import json
# data analysis tools
import pandas as pd
import numpy as np
# graph analysis tools
import networkx as nx
# ontology analysis tools
import ontobio as ob
# Custom Modules
from Parsing import parsing as pu
from Parsing import graphUtilities as gu
from EdgeBasedSimilarity import edgeBasedMethods as ebm

def edgeBasedMethodsMenu(geneData, ont):
    menu = '''
    ** Edge Based Methods **
    1. Simple Counting Method
    2. Simple Weighting Method
    3. Similarity GO
    '''
    choice = None
    while 1 :
        print(menu)
        try :
            choice = int(input('Select a type of methods : '))
            if choice < 1 or choice > 4:
                raise ValueError
            else:
                break
            #endif
        except ValueError :
            print('Give a valid integer !')
    #endwhile
    if choice==1 :
        ebm.simRada(geneData, ont)
    elif choice==2:
        ebm.simpleWeightedDistance(geneData , ont)
    elif choice==3:
        ebm.semanticValueSimilarity(geneData , ont)
    elif choice==4:
        ebm.shortestSemanticDifferentiationDistance(geneData , ont)
#
#
#
def main():
    '''
    1. Create appropriate format for GAF files and parse ontology
    1.1 Check if file exists
    1.2 Create/Read file
    '''
    '''
    pafn = input('Give processed annotation filename :')
    if '.json' in str(pafn):
        pass
    else:
        pafn+'.json'
    #endif
    geneData = None
    pafn = os.getcwd()+'/Datasets/'+pafn
    if os.path.exists(pafn):
        geneData = gu.read_json_file(pafn)
    else:# create annotation file
        geneData = pu.parseGAF(pafn)
    #endif
    #print(geneData)
    ontfile = input('Give ontology filename :')
    if '.obo' in str(ontfile):
        pass
    else:# u must fix this
        print('Error')
    #endif
    ontfile = os.getcwd()+'/Datasets/'+ontfile
    print('Creating Ontology from file .')
    ont = ob.OntologyFactory().create(ontfile)
    '''
    '''
    2. Main Menu
    2.1 Edge Based Methods
    2.2 Information Content Methods
    2.3 Hybrid Methods
    2.4 Machine Learning Methods
    '''
    pafn='genes.json'
    pafn=os.getcwd()+'/Datasets/'+pafn
    print('Reading genes from file .')
    geneData=gu.read_json_file(pafn)
    ontfile='go-basic.obo'
    ontfile=os.getcwd()+'/Datasets/'+ontfile
    print('Creating Ontology from file .')
    ont=ob.OntologyFactory().create(ontfile)
    menu = '''
    ** Main Menu **
    1. Edge Based Methods
    2. Information Content Methods
    3. Hybrid Methods
    4. Machine Learning Methods
    '''
    choice = None
    while 1 :
        print(menu)
        try :
            choice = int(input('Select a type of methods : '))
            if choice < 1 or choice > 4:
                raise ValueError
            #endif
            if choice==1:
                edgeBasedMethodsMenu(geneData , ont)
        except ValueError :
            print('Give a valid integer !')
    #endwhile

if __name__ == "__main__":
    main()
