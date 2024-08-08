import ontobio as ont
import networkx as nx

from Parsing import graphUtilities as gu

def averageLengthToLCA(t1 , t2 , lca , ont):
    G=ont.get_graph()
    path1 = list(nx.all_simple_paths(G, lca, t1))
    path2 = list(nx.all_simple_paths(G, lca, t2))
    return sum([len(i) for i in path1 + path2])/len(path1+path2)
#
#
#
def lowest_common_ancestor(t1, t2):
    '''
    1. Given two terms t1 , t2 find their LCA .
    2. Find maximum distance of LCA from root .
    '''
    if not(t1[1] == t2[1]):
        return [0,None]
    dist_1 = 0
    dist_2 = 0
    for i in t1[0]:
        for j in t2[0]:
            intersection=list(set(t1[0][i])&set(t2[0][j]))
            if len(intersection)>0:
                return [i+j , intersection]
        #endfor
    #endfor
#
#
#
def simpleWeightedDistance(t1 , t2 , root , anc1 , anc2 , ont):
    # 1. find lowest common ancestor
    lcaDist , lca = lowest_common_ancestor(anc1,anc2)
    # 2. if they are in the same namespace
    if not lca==None :
        print('-'*20)
        print(f'LCA : {lca[0]}')
        # 3. find lca's maximum path to root
        dist = gu.allAncestors(lca[0], root , ont)
        print(f'All ancestors : {str(dist)}')
        print(f'Distance from root : {sum([.815,] + [.815**(i+1) for i in range(1,len(dist.keys()))])}')
        avgL = averageLengthToLCA(t1,t2,lca[0],ont)
        print(f'Average path length : {avgL}')
        normalizationFactor = (avgL + sum([.815,] + [.815**(i+1) for i in range(1,len(dist.keys()))]))/avgL
        print('-'*20)
#
#
#
def minimumPath(t1 , t2):
    '''
    1.
    '''
    if not(t1[1] == t2[1]):
        return 0
    dist_1 = 0
    dist_2 = 0
    # 1. for each layer of t1
    for i in t1[0]:
        # 2. for each layer of t2
        for j in t2[0]:
            intersection=list(set(t1[0][i])&set(t2[0][j]))
            if len(intersection)>0:
                dist_1 = i
                dist_2 = j
                break
            else:
                continue
        #endfor
    #endfor
    return dist_1 + dist_2

def simRada(genesID, genesData , ancestors):
    '''
    1. for every possible combination of nodes find the minimum path between them .
    2. sum all the distances and divide by product of number of elements of each set .
    '''
    terms = list(ancestors.keys())
    simRada = {}
    for i in genesID:
        if not (i in list(simRada.keys())):
            simRada[i] = {}
        for j in genesID:
            if not(j in list(simRada.keys())):
                simRada[j] = {}
            if not(j in list(simRada[i].keys())):
                simRada[i][j] = 0
            else :
                #similarity is already in
                continue
            print('-'*20)
            # get genes
            g1=genesData[i]
            tset1 = [t[0] for t in g1]
            g2=genesData[j]
            tset2 = [t[0] for t in g2]
            dists = []
            normSet = []
            for t1 in tset1:
                for t2 in tset2:
                    normSet = [ancestors[t1][0][i] for i in ancestors[t1][0].keys()] + [ancestors[t2][0][i] for i in ancestors[t2][0].keys()]
                    print(str(normSet))
                    mPath=minimumPath(ancestors[t1],ancestors[t2])
                    dists.append(mPath)
                #endfor
            #endfor
            score = sum(dists)/(len(tset1)*len(tset2))
            print(f'Overall score of {i} and {j} : {score}')
            simRada[i][j] = score
            simRada[j][i] = score
            print('-'*20)
        #endfor
    #endfor
    return simRada
