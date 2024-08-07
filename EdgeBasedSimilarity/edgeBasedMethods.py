import ontobio as ont

'''
Method no.1 : Given two terms t1 ,t2 find their Lowest Common Ancestor .
              Count the edges with weight 1 .
'''
def minimumPath(t1 , t2):
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
'''
Method no.2 : Given two terms t1 , t2 find their LCA .
              Find distance of LCA from root .

'''
