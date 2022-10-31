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
    smiles_list = f.readlines()


def count_honeycomb(idx):
    chosen_smiles = smiles_list[idx]
    gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False)

    graph = mol_to_nx(Chem.MolFromSmiles(gold_smiles))
    honeycomb_graph = mol_to_nx(Chem.MolFromSmiles("C1CC2CCCC3CCCC(C1)C23"))
    honeycomb_graph2 = mol_to_nx(Chem.MolFromSmiles("C1CC2CCCC3CCC(C1)C23"))
    honeycomb_graph3 = mol_to_nx(Chem.MolFromSmiles("C1CC2CCC3CCCC4CCC(C1)C2C34"))
    honeycomb_graph4 = mol_to_nx(Chem.MolFromSmiles("C1CC2CCC3CCC4CCCC5CC(C1)C2C3C45"))
    honeycomb_graph4 = mol_to_nx(Chem.MolFromSmiles("C1CC2CCC3CCC4CCCC5CC(C1)C2C3C45"))
    honeycomb_graph5 = mol_to_nx(Chem.MolFromSmiles("C1CC2CCCC3CCC(C1)C23"))
    honeycomb_graph6 = mol_to_nx(Chem.MolFromSmiles("C1CC2CCC3CCC4CCC(C1)C2C34"))


    GM = iso.GraphMatcher(graph, honeycomb_graph)
    GM2 = iso.GraphMatcher(graph, honeycomb_graph2)
    GM3 = iso.GraphMatcher(graph, honeycomb_graph3)
    GM4 = iso.GraphMatcher(graph, honeycomb_graph4)
    GM5 = iso.GraphMatcher(graph, honeycomb_graph5)
    GM6 = iso.GraphMatcher(graph, honeycomb_graph5)
    if GM.subgraph_is_isomorphic() or GM2.subgraph_is_isomorphic() or GM3.subgraph_is_isomorphic() or \
        GM4.subgraph_is_isomorphic() or GM5.subgraph_is_isomorphic() or GM6.subgraph_is_isomorphic():
        with open("honeycomb_structure2.txt", "a") as myfile:
            print(gold_smiles, idx)
            myfile.writelines("{},{}\n".format(gold_smiles, idx))

    return

# with concurrent.futures.ThreadPoolExecutor() as executor:
#     executor.map(count_honeycomb, [i for i in range(len(smiles_list))])
for i in range(len(smiles_list)):
    count_honeycomb(i)
