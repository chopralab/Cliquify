import sys
sys.path.append('../')
from rdkit import Chem
import concurrent.futures
import networkx as nx
import multiprocessing

from Cliquify.enumerate3 import dfs_assemble, node_labelling, remove_edges_reset_idx, reconstruction_evaluation
from Cliquify.tree_decomposition2 import tree_decomp
from Cliquify.utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso

# with open("C:\\Users\\fongm\\Downloads\\icml18-jtnn\\data\\zinc\\all.txt") as f:
#     smiles_list = f.readlines()
with open("zinc\\all.txt") as f:
    smiles_list = f.readlines()


# smiles_list = smiles_list[14755:]


def decompose_reconstruct(idx):
    chosen_smiles = smiles_list[idx]
    mol = Chem.MolFromSmiles(chosen_smiles)
    cliques, molTreeEdges, triangulated_graph = tree_decomp(mol)
    try:
        root, cur_graph, global_amap = node_labelling(mol, cliques, molTreeEdges, triangulated_graph)
    except:
        gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False)
        print(gold_smiles, "fail_deconstruct")
        with open("space_time_benchmarking\\fail_mol_list.txt", "a") as myfile:
            myfile.writelines("{},{},Decompose\n".format(gold_smiles, idx))
        return
    
    dfs_assemble(cur_graph, global_amap, [], root, None, print_out=False)

    final_graph = remove_edges_reset_idx(cur_graph)
    gold_smiles, dec_smiles, graph_match  = reconstruction_evaluation(chosen_smiles, final_graph)
    
    if gold_smiles != dec_smiles or not graph_match:
        print(gold_smiles, "fail_reconstruct")
        with open("space_time_benchmarking\\fail_mol_list.txt", "a") as myfile:
            myfile.writelines("{},{},Reconstruct\n".format(gold_smiles, idx))
    
    return

with concurrent.futures.ThreadPoolExecutor() as executor:
    executor.map(decompose_reconstruct, [i for i in range(len(smiles_list))])