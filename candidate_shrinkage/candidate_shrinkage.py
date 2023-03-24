import sys
sys.path.append('../')
from rdkit import Chem
import concurrent.futures
import networkx as nx
import multiprocessing

from enumerate_shrinkage import dfs_assemble, node_labelling, remove_edges_reset_idx, reconstruction_evaluation
from tree_decomposition2 import tree_decomp
from utils import nx_to_mol, mol_to_nx, node_equal_iso2, ring_edge_equal_iso

with open("many_large_rings_increment.txt") as f:
    smiles_list = f.readlines()


def candidate_enumeration_count(idx):
    chosen_smiles = smiles_list[idx].split(",")[0]
    mol = Chem.MolFromSmiles(chosen_smiles)
    cliques, molTreeEdges, triangulated_graph = tree_decomp(mol)

    try:
        root, cur_graph, global_amap = node_labelling(mol, cliques, molTreeEdges, triangulated_graph)
    except:
        return
    
    enumerate_cand_per_node = []
    dfs_assemble(cur_graph, global_amap, [], root, None, enumerate_cand_per_node, print_out=False)

    # print(enumerate_cand_per_node)
    return enumerate_cand_per_node

with open("candidate_count.txt", "w") as myfile:
    for i in range(len(smiles_list)):
        enumerate_cand_per_node = candidate_enumeration_count(i)

        avg = sum(enumerate_cand_per_node) / len(enumerate_cand_per_node)
        print(avg, enumerate_cand_per_node)

        myfile.writelines("{}| {}\n".format(avg, enumerate_cand_per_node))


# with concurrent.futures.ThreadPoolExecutor() as executor:
#     executor.map(candidate_enumeration_count, [i for i in range(len(smiles_list))])