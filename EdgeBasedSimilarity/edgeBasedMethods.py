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
                print(intersection)
                break
            else:
                continue
        #endfor
    #endfor
    return dist_1 + dist_2

'''
Method no.2 : Given two terms t1 , t2 find their LCA .
              Find distance of LCA from root .

'''
