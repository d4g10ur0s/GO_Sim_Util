## General

This is a package for GO term similarity comparison . There are multiple methods from differect categories of bibliography concerning GO term similarity .

## Edge Based Methods 

#### 1. Simple Edge Counting Method

- Find the minimum lenght path for each pair of terms .
- The number of edges is normalized by the maximum distance .
- This metric can be used as similarity .

#### 2. Weight Edges using the minimum distance from LCA to the root

- Find the LCA for each pair of terms .
- Find the minimum distance to the root .
- Use the formula from :

$$
W_{path(t_1,t_2)} = \sum_{i=0}^{n}{(w)^i}
$$

- w is a hyperparameter for weight .
