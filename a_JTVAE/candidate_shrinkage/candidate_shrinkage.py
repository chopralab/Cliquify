import sys
sys.path.append('../')
from rdkit import Chem
import concurrent.futures
import networkx as nx
import multiprocessing

from jtnn.chemutils import *
from jtnn.mol_tree import MolTree

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

# with open("../../zdata/zinc/all.txt") as f:
#     smiles_list = f.readlines()

with open("../../candidate_shrinkage/many_large_rings_increment.txt") as f:
    smiles_list = f.readlines()


def candidate_enumeration_count(idx):
    chosen_smiles = smiles_list[idx].split(",")[0]
    mol = Chem.MolFromSmiles(chosen_smiles)


    tree = MolTree(chosen_smiles)
    tree.recover()

    label_idx = 0
    cur_mol = copy_edit_mol(tree.nodes[label_idx].mol)
    global_amap = [{}] + [{} for node in tree.nodes]
    global_amap[label_idx + 1] = {atom.GetIdx():atom.GetIdx() for atom in cur_mol.GetAtoms()}

    
    enumerate_cand_per_node = []
    dfs_assemble_shrinkage(cur_mol, global_amap, [], tree.nodes[label_idx], None, enumerate_cand_per_node)


    # print(enumerate_cand_per_node)
    return enumerate_cand_per_node

with open("candidate_count.txt", "a") as myfile:
    for i in range(len(smiles_list)):
        enumerate_cand_per_node = candidate_enumeration_count(i)

        avg = sum(enumerate_cand_per_node) / len(enumerate_cand_per_node)
        print(avg, enumerate_cand_per_node)

        myfile.writelines("{}| {}\n".format(avg, enumerate_cand_per_node))


# with concurrent.futures.ThreadPoolExecutor() as executor:
#     executor.map(candidate_enumeration_count, [i for i in range(len(smiles_list))])