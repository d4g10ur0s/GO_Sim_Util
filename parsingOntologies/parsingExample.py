import pronto
import sys
import random
import json

# Given two terms get their closest common ancestor
def findMinimumPath(term_1, term_2):
    '''
    0. if terms belong to different subontologies , then their distance is 0
    1. get parents of term_1 and parents of term_2
    2. if they have common anscestor , then return distance 1 else :
        3. get the grandparents of term_1 and term_2
        4. update distance by 1
        5. if these sets have common terms then return updated distance else :
        6. go to step 3.
    '''
    if not(term_1.namespace == term_2.namespace) :
        return 0

    distance = 1
    parents_1 = set(term_1.superclasses(distance = distance))
    parents_2 = set(term_2.superclasses(distance = distance))
    '''
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
    # add one more edge for each starting term
    return distance_1 + distance_2

def getAllMinimumPaths(go):
    dists = {}
    '''
    Input : an ontology
    Output : every subontologie's pair of terms distance
    '''
    counter = 0
    for t_1 in go.terms():
        for t_2 in go.terms():
            counter+=1
            '''
            if it is the same term continue
            if one of the terms is obsolete continue
            if distance is already captured
            '''
            if t_1.id == t_2.id or t_1.obsolete or t_2.obsolete :
                continue
            '''
            else :
                1. find the minimum path
                2. add minimum path in dictionary format
            '''
            dist = None
            if t_1.id in dists.keys():
                if not t_2.id in dists[t_1.id].keys():
                    dist=findMinimumPath(t_1,t_2)
                    dists[t_1.id][t_2.id] = dist
                    if t_2.id in dists.keys():
                        if not t_1.id in dists[t_2.id].keys():
                            dists[t_2.id][t_1.id] = dist
                    else:
                        dists[t_2.id] = {}
                        dists[t_2.id][t_1.id] = dist
            else:
                dist=findMinimumPath(t_1,t_2)
                dists[t_1.id] = {}
                dists[t_1.id][t_2.id] = dist
                if t_2.id in dists.keys():
                    if not t_1.id in dists[t_2.id].keys():
                        dists[t_2.id][t_1.id] = dist
                else:
                    dists[t_2.id] = {}
                    dists[t_2.id][t_1.id] = dist
            # some info
            if counter % 1000 == 0:
                print(f"Processed {counter} term pairs...")
    return dists

def parseOntology(path):
    go = pronto.Ontology(path)
    terms = go.terms()
    relationships = go.relationships()
    return go,terms,relationships

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <obo_file_path>")
        sys.exit(1)

    obo_path = sys.argv[1]
    # get ontologie's terms and relationships
    go,terms,relationships = parseOntology(obo_path)
    count = 0
    arr = []
    for t in terms :
        arr.append(t.id)
        print(str(t.name))
    # get two random terms
    print('Minimum distance between : ' + str(go.get_term(random.choice(arr)).name) + ' and ' + str(go.get_term(random.choice(arr)).name) + ' is ' + str(findMinimumPath(go.get_term(arr[0]),go.get_term(arr[1]))))
    paths = getAllMinimumPaths(go)
    with open('min_paths.json', 'w+') as outfile:
        json.dump(data, outfile, indent=2)


if __name__ == "__main__":
    main()
