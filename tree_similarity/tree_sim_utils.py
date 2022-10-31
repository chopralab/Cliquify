
import matplotlib.pyplot as plt
import numpy as np

def generate_index(n_target, n_sample):
    # 240_000 - 250K ZINC dataset
    target = np.random.randint(240_000, size=n_target)
    sample = np.random.randint(240_000, size=n_sample)
    
    # with open("index.txt") as f:
    #     for val in target: f.writelines("{},".format(val))
    #     for val in sample: f.writelines("{},".format(val))


def plot_ged():
    with open("../a_JTVAE/tree_similarity/ged_mean_values2.txt") as f:
        jtvae_ged_list = f.readlines()
        jtvae_ged = [float(ged.split(",")[1]) for ged in jtvae_ged_list if len(ged.split(",")) == 2]

    with open("ged_mean_values.txt") as f:
        cliquify_ged_list = f.readlines()
        cliquify_ged = [float(ged.split(",")[1]) for ged in cliquify_ged_list if len(ged.split(",")) == 2]

    plt.bar(x=[i for i, val in enumerate(cliquify_ged)], color="b", label="cliquify", height=cliquify_ged)
    plt.bar(x=[i for i, val in enumerate(jtvae_ged)], color="g", label="jtvae", height=jtvae_ged)
    # plt.xticks(rotation = 90)
    plt.title("Graph Edit Distance(GED) on 5000 smiles")
    plt.xlabel('smiles idx')
    plt.ylabel('ged score (nodes)')
    plt.legend()
    plt.savefig("tree_ged_score_compiled_trees.png")

plot_ged()