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
from InformationContentSimilarity import informationContentUtility as icu
from HybridSimilarity import hybridSimilarityUtils as hsu
#
#
#
def hybridMethodsMenu(geneData, ont):
    menu = '''
    ** Information Content Based Methods **
    1. Integrated Similarity Measure
    2. Hybrid Relative Specificity Similarity
    3. Semantic Weights Hybrid Similarity
    '''
    # 0. get information content
    prob = icu.calculateInformationContent(geneData , ont)
    prob.columns = ['terms' , 'probability']
    # 0.1 get frequency
    tFrequency = pd.read_csv(os.getcwd()+'/Datasets/termFrequency.csv')
    tFrequency.columns=['terms' , 'frequency']
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
        hsu.integratedSM(tFrequency , prob , ont)
    elif choice==2:
        hsu.calculateHRSS(prob , ont)
    elif choice==3:
        hsu.calculateSemanticWeightsHybridSimilarity(prob , ont)
#
#
#
def icBasedMethodsMenu(geneData, ont):
    menu = '''
    ** Information Content Based Methods **
    1. Resnik Semantic Similarity
    2. Jiang Semantic Similarity
    3. Lin Semantic Similarity
    4. Relevance Semantic Similarity
    5. Information Coefficient Semantic Similarity
    6. Graph Information Content Similarity
    '''
    # 0. get information content
    prob = icu.calculateInformationContent(geneData , ont)
    prob.columns = ['terms' , 'probability']
    print(prob)
    choice = None
    while 1 :
        print(menu)
        try :
            choice = int(input('Select a type of methods : '))
            if choice < 1 or choice > 7:
                raise ValueError
            else:
                break
            #endif
        except ValueError :
            print('Give a valid integer !')
    #endwhile
    if choice==1 :
        icu.calculateSimResnik(prob, ont)
    elif choice==2:
        icu.calculateSimJiang(prob , ont)
    elif choice==3:
        icu.calculateSimLin(prob , ont)
    elif choice==4:
        icu.calculateSimRel(prob , ont)
    elif choice==5:
        icu.calculateSimInfoCoeff(prob , ont)
    elif choice==6:
        icu.calculateSimGIC(prob , ont)
    elif choice==7:
        icu.calculateDiShin(prob , ont)
#
#
#
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
        #try :
        choice = int(input('Select a type of methods : '))
        if choice < 1 or choice > 4:
            raise ValueError
        #endif
        if choice==1:
            edgeBasedMethodsMenu(geneData , ont)
        if choice==2:
            icBasedMethodsMenu(geneData, ont)
        if choice==3:
            hybridMethodsMenu(geneData , ont)
        #except ValueError :
        #    print('Give a valid integer !')
    #endwhile

if __name__ == "__main__":
    main()
