# Cliquify

Cliquify - Robust representation of molecular graphs to trees structures is an extension to the work from [Junction Tree Variational Autoencoder (JTVAE)](https://github.com/wengong-jin/icml18-jtnn).

This work aims to improve the tree representation of the molecular graph by introducing a variation of Hugin's Algorithm through the formation of chordal graphs.

![image](https://user-images.githubusercontent.com/69520909/206056421-ec8e5e94-5449-41a2-9d6e-4c63e0c79dbf.png)

## Vocab Generalization
- We define the importance of tree molecular **vocabulary** through its ability of **representing more and diverse molecules**.
* JTVAE **ring vocabulary** constraints the number of molecules generated due to its poor generalizability.
+ Our solution fixes the problem by using more **generalizable** triangular cliques as vocabulary.
- Generalizable vocabulary helps in generative model, eg. VAE or GAN, to generate more diverse molecules without the need of redefining vocabulary based on new dataset.


### JTVAE (random sampled 2 vocabs)
- Vocabulary used

![image](https://user-images.githubusercontent.com/69520909/206058746-ad83ac75-d832-4f98-894c-63d58c54afb0.png)
- Samples generated (Pruned)

![image](https://user-images.githubusercontent.com/69520909/206060185-334a4a93-3c0f-428d-928b-2f82b8961d69.png)

### Cliquify (random sampled 2 vocabs)
- Vocabulary used

![image](https://user-images.githubusercontent.com/69520909/206060392-a9a72663-03bc-47b9-98a0-ae2c50392ae8.png)

- Samples generated (Pruned)

![image](https://user-images.githubusercontent.com/69520909/206077036-8c31488c-b7d8-4707-8fe6-6de4e5a141e3.png)

As you can from the comparison above, by using the more generalizable vocabulary from Cliquify, after random sampling, cliquify can produce molecules with diversified components.

## Candidate Shrinkage
Cliquify uses triangular clique decomposition, which helps in 
- reducing number of candidates per node generation (**candidate generation explosion** when involving large rings mentioned in [hgraph2graph](https://github.com/wengong-jin/hgraph2graph)),
* control the number and characteristics of candidates being generated for each fragments.

### Given sample molecule
![image](https://user-images.githubusercontent.com/69520909/206062937-f93d5c10-10bd-4ba3-9be0-20d7d122f041.png)

### JTVAE | Cliquify
--------
![image](https://user-images.githubusercontent.com/69520909/206063323-e26da9d7-ffcd-44b4-ab91-5209dbe4d20a.png)
- The diagram above shows how Cliquify reduces the possibility of candidates generation per node. 

![image](https://user-images.githubusercontent.com/69520909/206076233-2b2be06d-58b8-49e7-9c84-3b2112b4ea55.png)

* The diagram above shows the average candidate generation per tree node, from molecules which has 6 membered rings and above
- Cliquify has low fluctuation of average numbers of candidates generated as compared to JTVAE


## Tree Similarity 
* JTVAE junction tree (ring vocabulary) is not deterministic since there are potentially many molecules that correspond to the same junction tree. -
- Using Cliquify, using the triangulation clique method
  - increase deterministic properties of junction tree
  - allow lesser one to many relationship between junction tree and corresponding molecule
  ![image](https://user-images.githubusercontent.com/69520909/206075265-3c51bc64-2c0c-44a5-9c30-6984ab94a484.png)

+ We quantify the tree similarity between molecules using Graph Edit Distance (GED) from Networkx Library
  - GED based on tree nodes
  
  ![image](https://user-images.githubusercontent.com/69520909/206075722-26481893-5a13-40c5-b2da-3cd0518f2e51.png)
  
  - GED based on tree nodes and edges
  
  ![image](https://user-images.githubusercontent.com/69520909/206075820-57194205-654b-4c06-9988-736558e6a733.png)
  
* Based on the two diagrams above, we can infer that cliquify produces more unique trees as compared to JTVAE, 
making the tree structure more determistic for decoding, encourages more one to one relationship between molecules and tree structure representation.

## Descendant Orientation Awareness
- JTVAE – due to its neighborhood to neighborhood decoding process, it does not consider the orientation of the existing decoded molecule
+ Cliquify eliminate this possibility by restricting the location of possible attachment, reducing/eliminating the possibility of orientation identification error.
  - It does that through prioritizing Non Ring Bonds attachment during graph to tree decomposition, reducing the possible triangular cliques attached to the Non Ring Bond.

### Original molecule (JTVAE)
![image](https://user-images.githubusercontent.com/69520909/206077460-2192207a-87a2-4415-91e2-6fb7816094fe.png)
### Decoded molecule (JTVAE)
![image](https://user-images.githubusercontent.com/69520909/206077540-e511b7b5-a427-49e2-8d88-7e70dde2ea24.png)

## Honeycomb Problem
- Honeycomb structure is prevalent in large organic molecules. JTVAE fails to capture such formation
* Honeycomb formation requires recursive build, thus the more complicated the neighboring molecules, the larger the candidate count would be.
+ This would like result in possibility of candidate explosion.

![image](https://user-images.githubusercontent.com/69520909/206078608-769f7129-c8da-47b8-a658-6a4b40ca5e55.png)

### JTVAE 
- Due to its inherent tree structure decoding, JTVAE fails to capture how multiple children of the same parent are being connected to one another.

  ![image](https://user-images.githubusercontent.com/69520909/206078973-e7223a1c-25b9-46c7-bfd7-d4fc0a6030e0.png)

### Cliquify
- Cliquify is able to decompose the honeycomb structure, reducing the possibility of candidate explosion through pruning.

  ![image](https://user-images.githubusercontent.com/69520909/206079002-d5fe4567-4716-4412-9fb9-6a5bf9c9c233.png)


