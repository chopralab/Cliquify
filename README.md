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

![image](https://user-images.githubusercontent.com/69520909/206062058-07dd3dcf-ac30-4a43-8815-69ad775c4cd7.png)

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


* This diagram shows the average candidate generation per tree node, from molecules which has 6 membered rings and above
- Cliquify has low fluctuation of numbers of candidates generated as compared to JTVAE
![image](https://user-images.githubusercontent.com/69520909/206063589-9514c14d-bbcb-463b-8c6f-b59e246ec816.png)



