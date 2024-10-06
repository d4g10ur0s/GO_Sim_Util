import sys
import os
sys.path.append('/home/d4gl0s/diploma/Experiment')
import json
import random
import time
# data analysis tools
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
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
from MachineLearningMethods import anc2vec

def trainModel(ont):
    G = ont.get_graph().to_undirected()# graph must be undirected
    node2vec = Node2Vec(G, dimensions=100, walk_length=100, num_walks=20, workers=2)
    print('Training')
    model = node2vec.fit(window=10, min_count=1, batch_words=10)
    model.save(os.getcwd()+'/go2vec_model')
    return model
#
#
#
def createProteinMatrices(geneData , model):
    # Create a file for each protein
    start_time=time.time()
    counter=0
    for g in geneData.keys():
        gVec = []
        counter+=1
        icu.progressBar(counter, len(geneData.keys()), start_time)
        for term in geneData[g]:
            gVec.append(model.wv[term.strip()])
        #endfor
        # save with protein name to a csv
        proteinDataFrame = pd.DataFrame(data=gVec)
        proteinDataFrame.to_csv(os.getcwd()+'/Datasets/'+g+'.csv')
    #endfor
#
#
#
def executeModel(geneData,ont):
    ontfile='go-basic.obo'
    ontfile=os.getcwd()+'/Datasets/'+ontfile
    embeds = anc2vec.train.fit(ontfile, embedding_sz=200, batch_sz=64, num_epochs=100)
    print(embeds)
    '''
    model=None
    if not os.path.exists(os.getcwd()+'/go2vec_model'):
        model=trainModel(ont)
    else:
        model=gensim.models.Word2Vec.load(os.getcwd()+'/go2vec_model')
        #createProteinMatrices(geneData , model)
        pVectors=[]
        for g in geneData:
            pVec=pd.read_csv(os.getcwd()+'/Datasets/'+g+'.csv')
            pVectors.append(pVec.mean(axis=1))
        #endfor
        proteinDataFrame=pd.DataFrame(data=pVectors)
        centroids ,classes = my_kmeans(proteinDataFrame, 4 , 1000)
        for i in classes:
            print(i)
        #Plotting the results
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(class_1[2] , class_1[3] , class_1[4] , color = 'red')
        ax.scatter(class_2[2] , class_2[3] , class_2[4] , color = 'black')
        ax.set_xlabel(ks[2])
        ax.set_ylabel(ks[3])
        ax.set_zlabel(ks[4])
        plt.savefig("Graphs\\kMeansGraphs\\"+str(i)+"-Means-"+ks[2]+" "+ks[3]+" "+ks[4]+".png")
        plt.clf()
        '''
'''
JUST A DEMO

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
    node2vec = Node2Vec(G, dimensions=50, walk_length=40, num_walks=100, workers=2)
    print('Training')
    model = node2vec.fit(window=10, min_count=1, batch_words=10)
    model.wv.most_similar('GO:1990462')
    model.save('go2vec_model')
#
#
#
if __name__ == "__main__":
    main()
'''
def my_kmeans(df,k = 2,maxIterations=1000):
    #arxikopoihsh kentrweidwn
    always_centroid = []
    c1 = None
    c2 = None
    choose = np.random.randint(df.shape[1] , size=k)
    my_centroids = []
    for i in range(0,k):
        my_centroids.append( df[str(choose[i])].values.tolist() )
        always_centroid.append( df[str(choose[i])] )
    #ta exw kanei lista
    i = 0
    to_centroid = []
    for i in range(0,df.shape[1]):
        if i in choose:
            pass
        else:
            similarities = []
            for j in range(0,len(my_centroids)):
                #vazw tis omoiothtes se lista ke pairnw thn megaluterh apoluth timh
                similarities.append( my_cosine_similarity(np.squeeze( np.asarray(my_centroids[j] ) ) ,np.squeeze( np.asarray(df[str(i)].values.tolist() ) ) ) )
            #dialegw to megalutero similarity
            best = 0
            for j in range(0,len(similarities)):
                if abs(similarities[j]) > best:
                    best = similarities[j]
                    #prepei na kanw ke ena pop
                    if len(to_centroid)-1 == i:#to plh8os twn stoixeiwn einai iso me to i panta!1 kentroeides gia ka8e perilhpsh
                        to_centroid.pop(len(to_centroid) -1)
                    #to dianusma 8a paei sto kentroeides tade
                    to_centroid.append(j)
    iterations = -1
    while iterations < maxIterations:
        c1 = always_centroid#prin allaksei to kentroeides
        iterations+=1
        kappa = 0
        #update centroids
        for i in range(0,len(my_centroids)):#gia ka8e kedroeides
            for j in range(0,len(to_centroid)):
                #an eimai sto katallhlo kanw summ
                if to_centroid[j] == i:
                    #kane sum
                    always_centroid[i] = always_centroid[i]+df[str(j)]
                else:
                    pass
            #sto telos pollaplasiazw ola ta stoixeia
            always_centroid[i] = always_centroid[i]*(1/len(always_centroid[i]))
        #ksanakanw thn diadikasia ?
        my_centroids = []
        for i in range(0,k):
            my_centroids.append( always_centroid[i].values.tolist() )
        #ta exw kanei lista
        i = 0
        to_centroid = []
        for i in range(0,df.shape[1]):
            if i in choose:
                pass
            else:
                similarities = []
                for j in range(0,len(my_centroids)):
                    #vazw tis omoiothtes se lista ke pairnw thn megaluterh apoluth timh
                    similarities.append( my_cosine_similarity(np.squeeze( np.asarray(my_centroids[j] ) ) ,np.squeeze( np.asarray(df[str(i)].values.tolist() ) ) ) )
                #dialegw to megalutero similarity
                best = 0
                for j in range(0,len(similarities)):
                    if abs(similarities[j]) > best:
                        best = similarities[j]
                        if len(to_centroid)-1 == i:
                            to_centroid.pop(len(to_centroid) - 1)
                        to_centroid.append(j)
        c2 = my_centroids
        #an ta kedroeidh idia tote break
        p = True
        for i in range(0,k):
            print("Finished in : "+ str(iterations) +" iterations .")
            if c1[i].equals(c2[i]):
                pass
            else:
                p = False

    return (choose, to_centroid)
#
#
#
def my_cosine_similarity(arr1,arr2):
    dot = sum(a*b for a,b in zip(arr1,arr2) )
    norm_arr1 = sum(a*a for a in arr1) ** 0.5
    norm_arr2 = sum(b*b for b in arr2) ** 0.5
    csimilarity = dot/(norm_arr1*norm_arr2)
    return csimilarity
