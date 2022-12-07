import sys
sys.path.append('../')
from rdkit import Chem
import numpy as np

from jtnn.chemutils import *
from jtnn.mol_tree import MolTreeNode

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


with open("vocab.txt") as f:
    smiles_list = f.readlines()

fragments = ["C1CCCCC1", "C1CCCC1"]
# fragments = ["C1CCCCC1", "CC"]

np.random.seed(42)

MAX_COUNT = 3

def random_sample_tree(root, depth):
    global count
    global edges_list
    global nodes_dict

    if count > MAX_COUNT:
        return
        
    # neighbors_count = np.random.randint(0,MAX_COUNT5-depth) # max of 4 neigbors
    neighbors_count = np.random.randint(1, 4) # max of 4 neigbors


    for i in range(neighbors_count):

        if count > MAX_COUNT: return

        count += i + 1
        topology = np.random.rand()

        if topology > 0.5:
            idx = np.random.randint(0, len(fragments))
            neigh = MolTreeNode(fragments[idx])
            neigh.idx = count

            edges_list.append((root.idx, neigh.idx))
            nodes_dict[count] = neigh

            random_sample_tree(neigh, depth+1)

        else:
            neigh = MolTreeNode(fragments[1])
            neigh.idx = count

            edges_list.append((root.idx, neigh.idx))
            nodes_dict[count] = neigh



gen_smiles_list = []
for _ in range(1000):

    #--------------------------------------------------
    count = 0
    root = MolTreeNode(fragments[1])
    root.idx = 0
    nodes_dict = {root.idx: root}
    edges_list = []

    random_sample_tree(root, 0)

    for x,y in edges_list:
        nodes_dict[x].add_neighbor(nodes_dict[y])
        nodes_dict[y].add_neighbor(nodes_dict[x])

    
    list_of_nodes = list(nodes_dict.values())

    for i,node in enumerate(list_of_nodes):
        node.nid = i + 1
        node.is_leaf = (len(node.neighbors) == 1)
        if len(node.neighbors) > 1:
            set_atommap(node.mol, node.nid)

    #---------------------------------------------------

    root_idx = 0
    root = list_of_nodes[root_idx]

    cur_mol = copy_edit_mol(root.mol)
    global_amap = [{}] + [{} for node in list_of_nodes]
    global_amap[root_idx + 1] = {atom.GetIdx():atom.GetIdx() for atom in cur_mol.GetAtoms()}

    dfs_random_assemble(cur_mol, global_amap, [], root, None)


    dec_smiles = Chem.MolToSmiles(cur_mol, isomericSmiles=False)
    dec_mol  = Chem.MolFromSmiles(dec_smiles)

    cur_mol = cur_mol.GetMol()
    cur_mol = Chem.MolFromSmiles(Chem.MolToSmiles(cur_mol))

    if cur_mol:
        set_atommap(cur_mol)
        dec_smiles = Chem.MolToSmiles(cur_mol, isomericSmiles=False)
        print(dec_smiles)

        if dec_smiles not in gen_smiles_list:
            gen_smiles_list.append(dec_smiles)

    else:
        print(None)

with open("gen_mol_count5.txt", "a") as myfile:
    for gen_smiles in gen_smiles_list:
        myfile.writelines("{}\n".format(gen_smiles))
    # raise

