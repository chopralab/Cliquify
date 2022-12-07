
import matplotlib.pyplot as plt
import numpy as np

def plot_cand_shrinkage():

    jtvae_avg = []
    with open("../a_JTVAE/candidate_shrinkage/candidate_count.txt") as f:
        jtvae_cand_count_list = f.readlines()
        for cand_count_list in jtvae_cand_count_list:
            avg, cand_count = cand_count_list.split("| ")
            avg, cand_count = float(avg), eval(cand_count)
            jtvae_avg.append(avg)


    cliquify_avg = []
    with open("candidate_count.txt") as f:
        cliquify_cand_count_list = f.readlines()
        for cand_count_list in cliquify_cand_count_list:
            avg, cand_count = cand_count_list.split("| ")
            avg, cand_count = float(avg), eval(cand_count)
            cliquify_avg.append(avg)
    
    jtvae_avg = jtvae_avg[2000:]
    cliquify_avg = cliquify_avg[2000:]

    x = np.arange(len(cliquify_avg))

    plt.plot(x, jtvae_avg, color="g", label="jtvae")
    plt.plot(x, cliquify_avg, color="b", label="cliquify")

    plt.title("Avg Candidate Generation Count on {} smiles".format(len(cliquify_avg)))
    plt.xlabel('smiles idx')
    plt.ylabel('avg candidate count')
    plt.legend()
    plt.savefig("avg_candidate_generation_count.png")

    

plot_cand_shrinkage()
    