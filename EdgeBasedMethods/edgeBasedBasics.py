import pronto

import json
import random

'''
Utility no.1 : Find the lowest common ancestor of two terms .
               Returns the length of the path to LCA and the id of the LCA .
'''
def findMinimumPath(term_1, term_2):
    '''
    0. if terms belong to different subontologies , then their distance is 0
    '''
    if not(term_1.namespace == term_2.namespace) :
        return 0 , None
    '''
    1. get parents of term_1 and parents of term_2
    2. if they have common anscestor , then return distance 1
    '''
    distance = 1
    parents_1 = set(term_1.superclasses(distance = distance))
    parents_2 = set(term_2.superclasses(distance = distance))
    '''
    3. get the grandparents of term_1 and term_2
    4. update distance by 1
    5. if these sets have common terms then return updated distance else :
    6. go to step 3.

    Intuition :
        1. We travel the ontology graph in a bottom - up fashion
        2. For each edge we traverse add one to distance from node_i if the set of parents have more parents than before
        3. if the set has not more parent than before , then we have reached the root , no more addition , wait for the other node to find the root
    '''
    pastLen_1 = 1
    pastLen_2 = 1
    distance_1 = 0
    distance_2 = 0
    while len(parents_1.intersection(parents_2)) < 1 :
        distance+=1
        parents_1 = set(term_1.superclasses(distance = distance))
        parents_2 = set(term_2.superclasses(distance = distance))
        if len(parents_1) > pastLen_1 :
            pastLen_1 = len(parents_1)
            distance_1+=1
        if len(parents_2) > pastLen_2 :
            pastLen_2 = len(parents_2)
            distance_2+=1
    '''
    return number_of_edges , lowest_common_ancestors
    '''
    return distance_1 + distance_2 , list(parents_1.intersection(parents_2))

'''
Utility no.2 : Use distance of LCA from root .
               Returns least number of edges to root .
               There is an example for two random chosen terms .
'''
def getDistanceFromRoot(term):
    '''
    1. get LCA namespace
    2. travel the graph bottom - up
    3. stop the first time you meet an ancestor with term name as LCA's namespace
    '''
    rootTermName = term.namespace
    distance = 1
    parents = [x.name for x in term.superclasses(distance=distance)]
    while not term.namespace in parents :
        distance+=1
        parents = [x.name for x in term.superclasses(distance=distance)]
    return distance
'''
This is just an example
'''
def distanceFromRootExample(go):
    arr = [x.id for x in go.terms()]
    a = go.get_term(random.choice(arr))
    b = go.get_term(random.choice(arr))
    '''
    if terms are obsolete , do not use them !
    '''
    if b.obsolete or a.obsolete :
        print('Obsolete Terms !')
    else:
        dist , commonParents = findMinimumPath(a,b)
        print('Distance from root : ' + str(dist))
        print('Node weight : ' + str(sum([0.75**d for d in range(dist)])))
        print(str(commonParents))
# minimum distance from root method
def distanceFromRoot(go):
    arr = [x for x in go.terms()]
    termDict = {}
    counter=0
    for t1 in arr :
        for t2 in arr :
            counter+=1
            '''
            1. if one of terms is obsolete just continue .
            '''
            if t1.obsolete or t2.obsolete :
                continue
            else :
                '''
                2. find the minimum path .
                3. if it has been found earlier continue .
                '''
                if t1.id in termDict.keys():
                    if t2.id in termDict[t1.id].keys():
                        print('Has been checked .')
                        continue
                    else :
                        dist , commonParents = findMinimumPath(t1,t2) # if it is necessary , it is has high computational cost
                        if not commonParents==None:
                            parentDistFromRoot = []
                            for i in list(commonParents):
                                parentDistFromRoot.append(getDistanceFromRoot(i))
                            termDict[t1.id][t2.id] = {'dist' : dist , 'LCA' : [x.id for x in list(commonParents)],'LCA_dist' : parentDistFromRoot}
                            '''
                            Add second term info
                            '''
                            if t2.id in termDict.keys():
                                if t1.id in termDict[t2.id].keys():
                                    print('Has been checked .')
                                    continue
                                else :
                                    termDict[t2.id][t1.id] = {'dist' : dist , 'LCA' : [x.id for x in list(commonParents)],'LCA_dist' : parentDistFromRoot}
                            else :
                                termDict[t2.id] = {}
                                termDict[t2.id][t1.id] = {'dist' : dist , 'LCA' : [x.id for x in list(commonParents)],'LCA_dist' : parentDistFromRoot}
                else :
                    dist , commonParents = findMinimumPath(t1,t2) # if it is necessary , it is has high computational cost
                    if not commonParents==None :
                        termDict[t1.id] = {}
                        parentDistFromRoot = []
                        for i in list(commonParents):
                            parentDistFromRoot.append(getDistanceFromRoot(i))
                        termDict[t1.id][t2.id] = {'dist' : dist , 'LCA' : [x.id for x in list(commonParents)],'LCA_dist' : parentDistFromRoot}
                        '''
                        Add second term info
                        '''
                        if t2.id in termDict.keys():
                            if t1.id in termDict[t2.id].keys():
                                print('Has been checked .')
                                continue
                            else :
                                termDict[t2.id][t1.id] = {'dist' : dist , 'LCA' : [x.id for x in list(commonParents)],'LCA_dist' : parentDistFromRoot}
                        else :
                            termDict[t2.id] = {}
                            termDict[t2.id][t1.id] = {'dist' : dist , 'LCA' : [x.id for x in list(commonParents)],'LCA_dist' : parentDistFromRoot}
            # some info
            if counter % 1000 == 0:
                print(f"Processed {counter} term pairs...")
        # endfor
    #endfor
    return termDict
