import sys
sys.path.append('../')
from rdkit import Chem
import concurrent.futures
import networkx as nx
from multiprocessing import Pool

from enumerate3 import dfs_assemble, node_labelling, remove_edges_reset_idx, reconstruction_evaluation
from tree_decomposition2 import tree_decomp
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso, draw_mol


import timeit
with open("../zdata/zinc/all.txt") as f:
    smiles_list = f.readlines()



# filename = "fail_mol_list.txt"
# filename = "fail_mol_list_enum_reduced(multithreaded).txt"
filename = "fail_mol_list_enum_reduced(multiprocessing).txt"

def decompose_reconstruct(idx):
    chosen_smiles = smiles_list[idx]

    mol = Chem.MolFromSmiles(chosen_smiles)
    cliques, molTreeEdges, triangulated_graph = tree_decomp(mol)

    try:
        root, cur_graph, global_amap = node_labelling(mol, cliques, molTreeEdges, triangulated_graph)
    except:
        gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False)
        return "{},{},Decompose\n".format(gold_smiles, idx)
    
    dfs_assemble(cur_graph, global_amap, [], root, None, print_out=False)

    final_graph = remove_edges_reset_idx(cur_graph)
    # draw_mol(final_graph, 7778, folder="../space_time_benchmarking")

    gold_smiles, dec_smiles, graph_match  = reconstruction_evaluation(chosen_smiles, final_graph)

    if gold_smiles != dec_smiles or not graph_match:
        print(gold_smiles, dec_smiles, idx, "fail_reconstruct")
        return "{},{},Reconstruct\n".format(gold_smiles, idx)
    
    return 

def multiprocess_decode():
    pool = Pool()
    fail_list = pool.map(decompose_reconstruct, range(len(smiles_list)))

    return fail_list

def write_file(smiles_list):
    smiles_list = [smiles for smiles in smiles_list if smiles]
    with open(filename, "a") as myfile:
        myfile.writelines(smiles_list)


# for idx, smiles in enumerate(smiles_list):
#     decompose_reconstruct(idx)

print(timeit.Timer(multiprocess_decode).timeit(number=1))
