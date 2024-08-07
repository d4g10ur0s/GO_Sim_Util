# General

This is a package for GO term similarity comparison . There are multiple methods from differect categories of bibliography concerning GO term similarity .

# Edge Based Methods 

## 1. Simple Edge Counting Method 

- Source : [Development and application of a metric on semantic nets](https://ieeexplore.ieee.org/abstract/document/24528/)
- Implementation : ./EdgeBasedSimilarity/edgeBasedMethods.py

#### 1.1 Between Term Similarity

- Find the minimum lenght path .
- Each transition costs 1 .

#### 1.2 Between Entities

- For each combination of terms from each entity find the minimum path of terms .
- Sum all the distances and divide by product of number of elements of each set . 

## 2. Weight Edges using the minimum distance from LCA to the root

- Source : [A Knowledge-Based Clustering Algorithm Driven by
Gene Ontology](https://www.tandfonline.com/doi/abs/10.1081/BIP-200025659)
- Implementation : ./EdgeBasedSimilarity/edgeBasedMethods.py

#### 2.1 Between Term Similarity

- Find the LCA of terms .
- Find the maximum distance of LCA from the root .
- Use the formula from :

$$
W_{path(t_1,t_2)} = \sum_{i=0}^{n}{(w)^i}
$$

- w is a hyperparameter for weight . (.815)

#### 2.2 Between Entities
