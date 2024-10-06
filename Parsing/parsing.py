import sys
import json
sys.path.append('/home/d4gl0s/diploma/Experiment')
import os
import datetime as dt
import time
# data analysis tools
import pandas as pd
import numpy as np
# graph analysis tools
import networkx as nx
from node2vec import Node2Vec
import gensim
# ontology analysis tools
import ontobio as ob
# parallel programming
from dask.distributed import Client,LocalCluster,progress
import dask
# custom modules
from EdgeBasedSimilarity import edgeBasedMethods as ebm
from InformationContentSimilarity import informationContentUtility as icu
from HybridSimilarity import hybridSimilarityUtils as hsu
from Parsing import graphUtilities as gu
#
#
#
def my_cosine_similarity(arr1,arr2):
    dot = sum(a*b for a,b in zip(arr1,arr2) )
    norm_arr1 = sum(a*a for a in arr1) ** 0.5
    norm_arr2 = sum(b*b for b in arr2) ** 0.5
    csimilarity = dot/(norm_arr1*norm_arr2)
    return csimilarity
#
#
#
def vectorProteinVector(annotations , pairs , nm):
    start_time=time.time()
    model=gensim.models.Word2Vec.load(os.getcwd()+'/go2vec_model')
    # 1. for each pair , get it's annotations
    pVectors=[]
    counter=0
    for pair_index in range(len(pairs)):
        counter+=1
        icu.progressBar(counter , len(pairs), start_time)
        pair = pairs.iloc[pair_index]
        # 2. match with annotation
        annot_1 = [x[0] for x in annotations[pair['p1']]]
        annot_2 = [x[0] for x in annotations[pair['p2']]]
        # 3. get maximum simNormalized
        maxSim=[]
        # 3.1 find similarities for each pair of terms
        # check if t1 , t2 in same subontology
        for i in ['molecular_function' , 'cellular_component' , 'biological_process']:
            for t1 in annot_1:
                temp=[]
                if not(nm[t1]==i):
                    continue
                #endif
                for t2 in annot_2:
                    if nm[t1]==nm[t2]:
                        temp.append(my_cosine_similarity(model.wv[t1] , model.wv[t2]))
                    #endif
                #endfor
            #endfor
            if len(temp)==0:
                print('It was 0...')
                maxSim.append(0)
            else:
                maxSim.append(sum(temp)/len(temp))
        #endfor
        pVectors.append(maxSim)
    #endfor
    # Find the maximum length of the arrays
    max_len = max(len(array) for array in pVectors)
    # Create a list of lists with 0 padding for missing values
    padded_arrays = [[value for value in array] + [0] * (max_len - len(array)) for array in pVectors]
    # Create a DataFrame from the padded arrays
    df = pd.DataFrame(padded_arrays)
    #df.replace(0,df.mean(axis=0),inplace=True)
    df['svm']=pairs['svm']
    print(df)
    return df
#
#
#
def proteinVectors(annotations , pairs , similarity , nm):
    start_time=time.time()
    counter=0
    max_value = similarity.max().max()
    # Divide all values by the maximum value
    simNormalized = similarity.div(max_value, axis='index')
    # 1. for each pair , get it's annotations
    pVectors=[]
    for pair_index in range(len(pairs)):
        counter+=1
        icu.progressBar(counter , len(pairs), start_time)
        pair = pairs.iloc[pair_index]
        # 2. match with annotation
        annot_1 = [x[0] for x in annotations[pair['p1']]]
        annot_2 = [x[0] for x in annotations[pair['p2']]]
        # 3. get maximum simNormalized
        maxSim=[]
        # 3.1 find similarities for each pair of terms
        # check if t1 , t2 in same subontology
        for i in ['molecular_function' , 'cellular_component' , 'biological_process']:
            for t1 in annot_1:
                temp=[]
                if not(nm[t1]==i):
                    continue
                #endif
                for t2 in annot_2:
                    if nm[t1]==nm[t2]:
                        temp.append(simNormalized.loc[t1 , t2])
                    #endif
                #endfor
            #endfor
            if len(temp)==0:
                print('It was 0...')
                maxSim.append(0)
            else:
                maxSim.append(sum(temp)/len(temp))
        #endfor
        pVectors.append(maxSim)
    #endfor
    # Find the maximum length of the arrays
    max_len = max(len(array) for array in pVectors)
    # Create a list of lists with 0 padding for missing values
    padded_arrays = [[value for value in array] + [0] * (max_len - len(array)) for array in pVectors]
    # Create a DataFrame from the padded arrays
    df = pd.DataFrame(padded_arrays)
    #df.replace(0,df.mean(axis=0),inplace=True)
    df['svm']=pairs['svm']
    print(df)
    return df
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
def getProteinScore(geneData):
    dpath = 'protein_aliases.tsv'
    start_time = time.time()
    pairScore = pd.DataFrame(0, index=geneData.keys(), columns=geneData.keys(), dtype=np.float64)
    # 1. read interactions file
    with open(os.getcwd()+'/Datasets/STRING_interactions.txt','r') as f:
        lines = f.readlines()
    #endwith
    proteins = pd.read_csv(os.getcwd()+'/Datasets/'+dpath, sep='\t')
    counter=0
    for l in lines[1:]:
        counter+=1
        icu.progressBar(counter,len(lines),start_time)
        temp=l.split(' ')
        pNames1 = proteins[proteins['#string_protein_id']==temp[0]]['alias']
        pNames2 = proteins[proteins['#string_protein_id']==temp[1]]['alias']
        p1=None
        for p in pNames1:
            if p in geneData.keys():
                p1=p
                break
            #endif
        #endfor
        p2=None
        for p in pNames2:
            if p in geneData.keys():
                p2=p
                break
            #endif
        #endfor
        if p1==None or p2==None:
            continue
        #endif
        try :
            pairScore.loc[p1,p2]=float(int(temp[2])/1000)# divide by 1000
        except KeyError:
            print(f'Protein {p1} or {p2} do not exist .')
    #endfor
    pairScore.to_csv(os.getcwd()+'/Datasets/protein_interactions_extended.csv')
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
            try :
                temp = genesUniProt[genesUniProt['Entry Name']==p+'_HUMAN']
                if temp.empty or p in toSave.keys():
                    pass
                else:
                    go_terms = temp['Gene Ontology IDs'].values[0].split(';')
                    toSave[p]=[gt.strip() for gt in go_terms]
            except :
                pass
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
def mapInteractions(prot):
    interPath = '/home/d4gl0s/diploma/Experiment/Datasets/STRING_interactions.txt'
    pInter = pd.read_csv(interPath, sep=' ')
    interFiltered = pInter[pInter['protein1'].isin(prot)]
    interFiltered=interFiltered[interFiltered['protein2'].isin(prot)]
    return interFiltered
#
#
#
def mapGAF(geneData):
    if os.path.exists(os.getcwd()+'/Datasets/pair_scoring.csv'):
        toSave = pd.read_csv(os.getcwd()+'/Datasets/pair_scoring.csv')
        toSave.drop(columns=['Unnamed: 0'],inplace=True)
        return toSave
    #endif
    start_time=time.time()
    geneNames = geneData.keys()
    dpath = 'protein_aliases.tsv'#input('Give dataset name :')
    start_time = time.time()
    proteins = pd.read_csv(os.getcwd()+'/Datasets/'+dpath, sep='\t')
    validProteins = proteins[proteins['alias'].isin(geneNames)].reset_index(drop=True)
    filteredProt = validProteins['#string_protein_id'].unique()
    interFiltered=mapInteractions(filteredProt)
    # construct dataframe having p1 , p2 , score
    toSave = pd.DataFrame(columns=['p1', 'p2' , 'score'])
    counter=0
    for i in range(len(interFiltered)):
        counter+=1
        icu.progressBar(counter ,len(interFiltered), start_time)
        pair = interFiltered.iloc[i]
        name_1 = validProteins[validProteins['#string_protein_id']==pair['protein1']]['alias'].unique()
        name_2 = validProteins[validProteins['#string_protein_id']==pair['protein2']]['alias'].unique()
        temp = (name_1[0] , name_2[0] , float(pair['combined_score'])/1000)
        toSave.loc[len(toSave.index)] = temp
    #endfor
    toSave.to_csv(os.getcwd()+'/Datasets/pair_scoring.csv')
    return toSave

#
#
#
def extractTerms(genes):
    terms=[]
    for g in genes:
        geneTerms=genes[g]
        for t in geneTerms:
            terms.append(t[0])
        #endfor
    #endfor
    return list(set(terms))
#
#
#
def extractTermsFromGenes(genes):
    terms = []
    for g in genes:
        geneTerms = genes[g]
        for t in geneTerms:
            terms.append(t.strip())
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
