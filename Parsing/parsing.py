import sys
import json
sys.path.append('/home/d4gl0s/diploma/Experiment')
import os
import datetime as dt
import time
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
from InformationContentSimilarity import informationContentUtility as icu
from HybridSimilarity import hybridSimilarityUtils as hsu

#
#
#
def save_to_json(gene_dict, file_path):
    '''
    Just save a file in json format
    '''
    with open(file_path, 'w') as f:
        json.dump(gene_dict, f, indent=4)
#
#
#
def mapU2S(genesUniProt):
    dpath = 'protein_aliases.tsv'#input('Give dataset name :')
    start_time = time.time()
    if os.path.exists(os.getcwd()+'/Datasets/'+dpath):
        # 1. get proteins
        proteins = pd.read_csv(os.getcwd()+'/Datasets/'+dpath, sep='\t')
        # 2. extract names and corresponding go terms
        pnames = proteins['alias'].values.tolist()
        toSave = {}
        print(proteins['#string_protein_id'])
        counter=0
        for p in pnames:
            counter+=1
            icu.progressBar(counter , len(pnames) , start_time)
            temp = genesUniProt[genesUniProt['Entry Name']==p+'_HUMAN']
            if temp.empty or p in toSave.keys():
                pass
            else:
                go_terms = temp['Gene Ontology IDs'].values[0].split(';')
                toSave[p]=go_terms
        #endfor
        print(toSave)
        save_to_json(toSave, os.getcwd()+'/Datasets/termsU2S.json')
    else:
        return None
#
#
#
def parseUniProt():
    dpath = 'uniprotkb.tsv'#input('Give dataset name :')
    if os.path.exists(os.getcwd()+'/Datasets/genesUniProtKB.tsv'):
        toret = pd.read_csv(os.getcwd()+'/Datasets/genesUniProtKB.tsv')
        toret.drop(columns=['Unnamed: 0'])
        return toret
    #endif
    if os.path.exists(os.getcwd()+'/Datasets/'+dpath):
        # 1. get proteins
        proteins = pd.read_csv(os.getcwd()+'/Datasets/'+dpath, sep='\t')
        # 2. extract names and corresponding go terms
        rFiltered = proteins[proteins['Reviewed']=='reviewed']# Get only reviewed proteins
        rFiltered.dropna(inplace=True)# clear from NaN values from Gene Ontology IDs
        toSave = rFiltered[['Entry Name', 'Gene Ontology IDs']]
        # 3. create json file
        print(toSave)
        toSave.to_csv(os.getcwd()+'/Datasets/genesUniProtKB.tsv')
        return toSave
    else:
        return None
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
                    if not (asc.subject.label) in geneDict.keys():
                        geneDict[asc.subject.label] = []
                    geneDict[asc.subject.label].append((asc.object.id.namespace+':'+asc.object.id.identity , str(asc.evidence.gaf_evidence_code())))
            except (AttributeError,TypeError) as e :
                print(str(e))
            finally :
                counter+=1
                if counter%10000 == 0:
                    print(f'Terms processed : {counter}')
        #endwhile
        return geneDict
#
#
#
def getStringAssoc(genes_path , string_path):
    with open(genes_path, 'r') as f:# get genes' names
        data = json.load(f)
    # create a dictionary of the form gID:gName
    with open(string_path , 'r') as fstring :
        lines = fstring.readlines()
        stringGenes = {}
        for l in lines:
            gene = l.split()
            if gene[1] in data.keys():
                stringGenes[gene[0]]=gene[1]
            #endif
        #endfor
        save_to_json(stringGenes,"./Datasets/stringGenes.json")
        return stringGenes
#
#
#
def constuctProteinInteractions(stringAssoc):
    # 1. read file and split lines accordingly
    with open("./Datasets/STRING_interactions.txt", 'r') as f:
        pInter = f.readlines()
        pInter.pop(0)# first line is nothing useful
        temp=[]
        interSet = []
        for l in pInter:
            # split lines
            sl = l.split(' ')
            if sl[0] in stringAssoc.keys() and sl[1] in stringAssoc.keys():
                vInter = int(sl[2])/1000
                interSet.append((stringAssoc[sl[0]] ,stringAssoc[sl[1]] ,vInter))
            #endfor
        #endfor
        # save as a dataframe of values
        interSetDF = pd.DataFrame(data=interSet , columns=['p_name_1' , 'p_name_2' , 'score'])
        interSetDF.to_csv('./Datasets/protein_interactions.csv')
        return interSetDF
#
#
#
def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <obo_file_path> <ga_file_path> <STRING_protein_info>")
        sys.exit(1)
    # path for .obo
    obo_path = sys.argv[1]
    # path for .gaf
    ga_path = sys.argv[2]
    string_path = sys.argv[3]
    # create ontology
    pafn='genes.json'
    pafn=os.getcwd()+'/Datasets/'+pafn
    print('Reading genes from file .')
    # gene pairwise similarity
    with open("./Datasets/genes.json", 'r') as f:# get genes' names
        geneData = json.load(f)
    #endfor
    print('Creating Ontology from file .')
    ont = ob.OntologyFactory().create(obo_path)
    prob = icu.calculateInformationContent(geneData , ont)
    prob.columns = ['terms' , 'probability']
    print(prob)
    #geneDict = parseGAF(ga_path)
    #save_to_json(geneDict,"./Datasets/genes.json")
    # match STRING IDS to GENE NAMES
    stringAssoc=None
    if os.path.exists("./Datasets/stringGenes.json"):
        with open("./Datasets/stringGenes.json", 'r') as f:
          stringAssoc = json.load(f)
    else:
        stringAssoc = getStringAssoc("./Datasets/genes.json",string_path)
    # find interactions based on valid string ids
    stringInteractions = None
    if os.path.exists('./Datasets/protein_interactions.csv'):
        stringInteractions=pd.read_csv('./Datasets/protein_interactions.csv')
    else:
        stringInteractions=constuctProteinInteractions(stringAssoc)
    #endif
    stringInteractions.drop(columns=['Unnamed: 0'], inplace=True)
    # get Resnik Similarity
    simRes = icu.calculateSimRel(prob, ont)
    # get first two genes
    for inter in range(int(len(stringInteractions))):
        gene_1 = stringInteractions.loc[inter , 'p_name_1']
        gene_2 = stringInteractions.loc[inter , 'p_name_2']
        geneSimilarity = []
        for t1 in geneData[gene_1]:
            tsim=[]
            for t2 in geneData[gene_2]:
                tsim.append(simRes.loc[t1[0] , t2[0]])
            #endfor
            if len(tsim)>0:
                geneSimilarity.append(sum(tsim)/len(tsim))
            else:
                geneSimilarity.append(0)
        #endfor
        print(f'Similarity of Genes {gene_1} and {gene_2} : {sum(geneSimilarity)/len(geneSimilarity)}')
        print(f'Score of Genes {gene_1} and {gene_2} : {stringInteractions.loc[inter , 'score']}')


if __name__ == "__main__":
    main()
