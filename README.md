# A. General

This is a package for GO term similarity comparison . There are multiple methods from differect categories of bibliography concerning GO term similarity .

# B. Set Up

## 1. Java - 11

- Install : `sudo apt-get install openjdk-11-jdk`
- Add to path : `export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64`

## 2. OWL Tools : 

- Source : [owltools](https://github.com/owlcollab/owltools)

1. Get from GitHub : `git clone https://github.com/owlcollab/owltools.git`
2. Move to `OWLTools-Parent` directory .
3. Install using maven : Run `mvn clean install` .
4. Add owltools to path : `export PATH=$PATH:directory/to/where/you/downloaded/owltools/OWLTools-Runner/target`

## 3. Python3 Modules

- ontobio : `pip install ontobio`
- networkx : `pip install networkx[default,extra]`

# C. Methods

## Edge Based Methods 

### 1. Simple Edge Counting Method 

- Source : [Development and application of a metric on semantic nets](https://ieeexplore.ieee.org/abstract/document/24528/)
- Implementation : ./EdgeBasedSimilarity/edgeBasedMethods.py

#### 1.1 Between Term Similarity

- Find the minimum lenght path .
- Each transition costs 1 .

##### 1.2 Between Entities

- For each combination of terms from each entity find the minimum path of terms .
- Sum all the distances and divide by product of number of elements of each set . 

### 2. Weight Edges using the maximum distance from LCA to the root

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

### 3. Semantic Value Method

- Source : [A new method to measure the semantic similarity of GO terms](https://academic.oup.com/bioinformatics/article/23/10/1274/197095)
- Implementation : ./EdgeBasedSimilarity/edgeBasedMethods.py

#### 3.1 Between Term Similarity

- Construct the subgraph for each term .
- Calculate svalue using the formula :
  
$$
S_{t_i}(t_i) = \begin{cases}S_{t_i}(t_i) = 1 \\ S_t(t_i) = \max{\{w_e\cdot S_{t_i}(t_i')|t_i' \in childrenof(t_i)\}} \ if \ t \ \neq \ t_i\end{cases}
$$

- Calculate term similarity using the formulas :

$$
SV(t) = \sum_{t_i \in N_t}{S_t(t_i)}
$$

$$
S_{GO}(t_1, t_2) = \frac{\sum_{t_i \in N_{t_1}\cap N_{t_2}}{(S_{t_1}(t_i) + S_{t_2}(t_i)}}{SV(t_1) + SV(t_2)}
$$

### 4. Shortest Semantic Differentiation Distance - SSDD

- Source : [A novel insight into Gene Ontology semantic similarity](https://www.sciencedirect.com/science/article/pii/S0888754313000876)
- Implementation : ./EdgeBasedSimilarity/edgeBasedMethods.py

#### 4.1 Between Term Similarity

---

## Information Content Methods 

### 1. Resnik Semantic Similarity

### 2. Jiang & Conrath Semantic Similarity

### 3. Lin Semantic Similarity

### 4. Relevance Semantic Similarity

### 5. SimIC

---

## Hybrid Methods 

### 1. Integrated Semantic Similarity


### 2. Hybrid Relative Specificity Similarity

- Source : [Improving the Measurement of Semantic Similarity between Gene Ontology Terms and Gene Products: Insights from an Edge- and IC-Based Hybrid Method](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0066745)
- 

#### 2.1 Between Term Similarity

- Find Most Informative Common Ancestor - MICA .
- Calculate Information Content Distance from root .

$$
\alpha_{IC} = -\ln{(p_{IC}(MICA)}
$$

- Find Most Informative Leaf for each term .
- Calculate Average Information Content Distance from each MIL and term .

$$
\beta_{IC} = \frac{-\ln{(p_{IC}(t_1))} + \ln{(p_{IC}(MIL_1))} - \ln{(p_{IC}(t_2))} + \ln{(p_{IC}(MIL_2))} }{2}
$$

- Find Minimum Path from each term to MICA .
- 

---

## Graph Embedding Methods 

### 1. go2vec
