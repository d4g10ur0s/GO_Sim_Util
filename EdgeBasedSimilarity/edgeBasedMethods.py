import ontobio as ont
import networkx as nx

from Parsing import graphUtilities as gu

def getSvalue(t,ont,rootNodes,namespace):
    # 2. construct its graph
    node = ont.node(t)
    root = None
    for i in node['meta']['basicPropertyValues']:
        if i['val'] in namespace :
            root = rootNodes[i['val']]
    graph = gu.allAncestors(t,root,ont)
    # 3. calculate all semantic values for the graph
    sval = {}
    weightFactor = .815
    # starting node's semantic value is 1
    sval[t] = 1
    for i in graph :
        lterms = graph[i]
        # semantic value depends on children
        for p in lterms :
            # if first layer , then term's child is the term under process
            if i==1 :
                sval[p] = 1 * weightFactor
            elif len(lterms)>0 :
                # 3.1 get node's children in ontology graph
                children = ont.children(p)
                # 3.2 get the intersection with all the previous layers
                player = graph[i-1]
                childIntersection = set(sval.keys())&set(children)
                sval[p]=0
                for c in childIntersection:
                    # 3.3 find the maximum svalue
                    if sval[p]<sval[c] * weightFactor:
                        sval[p] = sval[c] * weightFactor
                    #endif
                #endfor
            #endif
        #endfor
    #endfor
    #print(f'Term {t} Semantic Value : {str(sum([sval[i] for i in sval.keys()]))}')
    return sval , sum([sval[i] for i in sval.keys()])
#
#
#
def pairSemanticValue(t1, t2 , ont):
    rootNodes={}
    namespace = []
    for r in ont.get_roots():
        rootNodes[str(ont.node(r)['label'])] = r
        namespace.append(str(ont.node(r)['label']))
    sval1 , semanticValue1 = getSvalue(t1 , ont , rootNodes , namespace)
    sval2 , semanticValue2 = getSvalue(t2 , ont , rootNodes , namespace)
    # Calculate semantic similarity of t1 , t2
    inTerms = set(sval1.keys()) & set(sval2.keys())
    svalPair = sum([sval1[i]+sval2[i] for i in inTerms])
    pairSimilarity = svalPair/(semanticValue1+semanticValue2)
    print(f'Similarity of {t1} , {t2} : {pairSimilarity}')
    return pairSimilarity
#
#
#
def geneSemanticValue(g1 , g2 , gene1 , gene2 , ont):
    # 1. for each possible term pair get the maximum similarity
    sim = 0
    for term1 in gene1 :
        tsim = 0
        for term2 in gene2 :
            psv=ebm.pairSemanticValue(t1, t2, ont)
            if tsim < psv:
                tsim=psv
            #endif
        #endfor
        sim+=tsim
    #endfor
    for term2 in gene2 :
        tsim = 0
        for term1 in gene1 :
            psv=ebm.pairSemanticValue(t1, t2, ont)
            if tsim < psv:
                tsim=psv
            #endif
        #endfor
        sim+=tsim
    #endfor
    # use the formula to calculate gene similarity
    print(f'Similarity between genes {g1} , {g2} : {sim/(len(gene1) + len(gene2))}')
    return sim/(len(gene1) + len(gene2))
#
#
#
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
