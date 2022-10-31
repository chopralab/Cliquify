import sys
sys.path.append('../')
from rdkit import Chem
import concurrent.futures
import networkx as nx
import multiprocessing
import networkx.algorithms.isomorphism as iso

from enumerate3 import dfs_assemble, node_labelling, remove_edges_reset_idx, reconstruction_evaluation
from tree_decomposition2 import tree_decomp
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso

with open("../zinc/all.txt") as f:
    original_smiles_list = f.readlines()

with open("honeycomb_structure.txt") as f:
    smiles_list = f.readlines()

honeycomb_solve = {
    "solved": 0,
    "partial_solved": 0,
    "unsolved": 0,
    "deconstruct": 0,
}

def solve_honeycomb(idx):
    smiles_info = smiles_list[idx]
    _, idx = smiles_info.split(",")
    chosen_smiles = original_smiles_list[int(idx)]

    mol = Chem.MolFromSmiles(chosen_smiles)
    cliques, molTreeEdges, triangulated_graph = tree_decomp(mol)
    try:
        root, cur_graph, global_amap = node_labelling(mol, cliques, molTreeEdges, triangulated_graph)
    except:
        gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False)
        print(gold_smiles, "fail_deconstruct")
        honeycomb_solve["deconstruct"] += 1
        return
    
    dfs_assemble(cur_graph, global_amap, [], root, None, print_out=False)

    final_graph = remove_edges_reset_idx(cur_graph)
    gold_smiles, dec_smiles, graph_match  = reconstruction_evaluation(chosen_smiles, final_graph)
    
    if gold_smiles != dec_smiles or not graph_match:

        graph1 = mol_to_nx(Chem.MolFromSmiles(gold_smiles))
        graph2 = mol_to_nx(Chem.MolFromSmiles(dec_smiles))
        GM = iso.GraphMatcher(graph1, graph2, node_match=node_equal_iso2)


        if GM.is_isomorphic():
            honeycomb_solve["partial_solved"] += 1
        else:
            # print("unsolved", idx)
            # print(gold_smiles)
            # print(dec_smiles)
            honeycomb_solve["unsolved"] += 1
    else:
        print("solved", idx)
        honeycomb_solve["solved"] += 1
        
    return

# for i in range(len(smiles_list)):
#     solve_honeycomb(i)
with concurrent.futures.ThreadPoolExecutor() as executor:
    executor.map(solve_honeycomb, [i for i in range(len(smiles_list))])

print(honeycomb_solve)