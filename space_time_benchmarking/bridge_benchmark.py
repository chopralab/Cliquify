import sys
sys.path.append('../')
from rdkit import Chem
import concurrent.futures
import networkx as nx
import multiprocessing

from enumerate3 import dfs_assemble, node_labelling, remove_edges_reset_idx, reconstruction_evaluation
from tree_decomposition2 import tree_decomp
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso, draw_mol


with open("../candidate_shrinkage/bridge_rings.txt") as f:
    smiles_list = f.readlines()


filename = "fail_mol_list_bridge_decomp.txt"

def decompose_reconstruct(idx):
    chosen_smiles = smiles_list[idx]
    # print(chosen_smiles)

    mol = Chem.MolFromSmiles(chosen_smiles)
    cliques, molTreeEdges, triangulated_graph = tree_decomp(mol)

    # for k, v in molTreeEdges:
    #     print(cliques[k], cliques[v])
    # raise

    try:
        root, cur_graph, global_amap = node_labelling(mol, cliques, molTreeEdges, triangulated_graph)
    except:
        gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False)
        print(gold_smiles, idx, "fail_deconstruct")
        with open(filename, "a") as myfile:
            myfile.writelines("{},{},Decompose\n".format(gold_smiles, idx))
        return
    
    dfs_assemble(cur_graph, global_amap, [], root, None, print_out=False)

    final_graph = remove_edges_reset_idx(cur_graph)
    # draw_mol(final_graph, 7778, folder="../space_time_benchmarking")

    gold_smiles, dec_smiles, graph_match  = reconstruction_evaluation(chosen_smiles, final_graph)
    
    # if chosen_smiles != dec_smiles:
    #     print(idx, dec_smiles, "chirality")

    if gold_smiles != dec_smiles or not graph_match:
        print(gold_smiles, dec_smiles, idx, "fail_reconstruct")

        # with open(filename, "a") as myfile:
        #     myfile.writelines("{},{},Reconstruct\n".format(gold_smiles, idx))
    
    return

for i in range(len(smiles_list)):
    if i % 100 == 0: print(i)
    decompose_reconstruct(i)

    # if i == 3:
    #     decompose_reconstruct(i)
    #     raise

# with concurrent.futures.ThreadPoolExecutor() as executor:
#     executor.map(decompose_reconstruct, [i for i in range(len(smiles_list))])