import pronto

import json
import random

import sys

from EdgeBasedMethods import edgeBasedBasics as eb

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
                    dist=eb.findMinimumPath(t_1,t_2)
                    dists[t_1.id][t_2.id] = dist
                    if t_2.id in dists.keys():
                        if not t_1.id in dists[t_2.id].keys():
                            dists[t_2.id][t_1.id] = dist
                    else:
                        dists[t_2.id] = {}
                        dists[t_2.id][t_1.id] = dist
            else:
                dist=eb.findMinimumPath(t_1,t_2)
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
    #print('Minimum distance between : ' + str(go.get_term(random.choice(arr)).name) + ' and ' + str(go.get_term(random.choice(arr)).name) + ' is ' + str(eb.findMinimumPath(go.get_term(arr[0]),go.get_term(arr[1]))))
    print('''
    ** Menu **
    1. Edge - Based Methods
    2. Information Content Methods
    3. Hybrid Methods
    4. Machine Learning

    Give -1 to exit
    ''')
    while 1 :
        try :
            choice = int(input('Choose Type of Method : '))
            if choice == 1:
                while 1 :
                    print('''
                    ** Edge - Based Methods **
                    1. Minimum Distance
                    2. Maximum Distance from Root
                    ''')
                    methodChoice = int(input('Choose Method : '))
                    if methodChoice == 1:
                        paths = getAllMinimumPaths(go)
                        with open('min_paths.json', 'w+') as outfile:
                            json.dump(data, outfile, indent=2)
                    if methodChoice == 2:
                        eb.distanceFromRootExample(go)
                        eb.distanceFromRoot(go)
                        with open('dist_from_root.json', 'w+') as outfile:
                            json.dump(data, outfile, indent=2)
                    else :
                        print('Give a valid number .')
            elif choice == -1 :
                print('Bye !')
                break
            else :
                print('Give a valid number .')
        except ValueError :
            print('Insert a valid number')

if __name__ == "__main__":
    main()
