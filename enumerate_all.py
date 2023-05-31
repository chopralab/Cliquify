
from rdkit import Chem
import itertools
import networkx as nx
from utils import *
import networkx.algorithms.isomorphism as iso
from debug_script import *
from rdkit import RDLogger
import numpy as np
from tree_decomposition2 import tree_decomp
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
import copy

np.random.seed(42)
KEY = "ghost"

def copy_atom(atom):
    new_atom = Chem.Atom(atom.GetSymbol())
    new_atom.SetFormalCharge(atom.GetFormalCharge())
    new_atom.SetAtomMapNum(atom.GetAtomMapNum())
    new_atom.SetChiralTag(atom.GetChiralTag())
    new_atom.SetHybridization(atom.GetHybridization())
    new_atom.SetNumExplicitHs(atom.GetNumExplicitHs())
    new_atom.SetIsAromatic(atom.GetIsAromatic())
    return new_atom

def set_atommap(mol, num=0):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(num)

def set_atommap_graph(G, num=0):
    for node in G.nodes():
        G.nodes[node]["map_num"] = num

def get_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: 
        return None
    Chem.Kekulize(mol)
    return mol

def get_triangulated_graph(mol):
    if mol and mol.GetNumAtoms() == 3 and mol.GetNumBonds() <= 2:
        triangulated_mol = Chem.RWMol(Chem.MolFromSmiles(''))
        for atom in mol.GetAtoms():
            new_atom = copy_atom(atom)
            triangulated_mol.AddAtom(new_atom)

        subset_possible_bonds = list(itertools.combinations(mol.GetAtoms(), 2))
        subset = [(bond[0].GetIdx(), bond[1].GetIdx()) for bond in subset_possible_bonds]

        # G = nx.Graph()
        for bond in subset:
            a1, a2 = bond[0], bond[1]
            bond_obj = mol.GetBondBetweenAtoms(a1, a2)
            if bond_obj:
                triangulated_mol.AddBond(a1, a2, order=bond_obj.GetBondType())
                new_bond = triangulated_mol.GetBondBetweenAtoms(a1, a2)
                new_bond.SetBoolProp(KEY, False)
                # G.add_edge(a1, a2, order=bond_obj.GetBondTypeAsDouble())
            else:
                triangulated_mol.AddBond(a1, a2, order=Chem.BondType.SINGLE)
                new_bond = triangulated_mol.GetBondBetweenAtoms(a1, a2)
                new_bond.SetBoolProp(KEY, True)
                # G.add_edge(a1, a2, order=0.0)        

        mol = triangulated_mol.GetMol()
        mol.UpdatePropertyCache() # getNumImplicitHs() called without preceding call to calcImplicitValence()
        G = mol_to_nx(mol)

        return mol, G
    else:
        for bond in mol.GetBonds():
            bond.SetBoolProp(KEY, False)
        G = mol_to_nx(mol)
        return mol, G

def get_fragments(mol, atoms):
    try: return Chem.MolFragmentToSmiles(mol, atoms, kekuleSmiles=True)
    except: return Chem.MolFragmentToSmiles(mol, atoms)

def get_fragments2(mol, atoms):
    try: return Chem.MolFragmentToSmiles(mol, atoms, kekuleSmiles=True, allHsExplicit=True)
    except: 
        return Chem.MolFragmentToSmiles(mol, atoms, allHsExplicit=True)

def get_smarts_fragments(mol, atoms):
    return Chem.MolFragmentToSmarts(mol, atoms)

def get_mol2(mol, clique):
    new_mol = Chem.RWMol()

    node_to_idx = {}

    for idx in clique:
        symbol=mol.GetAtomWithIdx(idx).GetSymbol()
        formal_charge=mol.GetAtomWithIdx(idx).GetFormalCharge()
        chiral_tag=mol.GetAtomWithIdx(idx).GetChiralTag()
        hybridization=mol.GetAtomWithIdx(idx).GetHybridization()
        num_explicit_hs=mol.GetAtomWithIdx(idx).GetNumExplicitHs()
        is_aromatic=mol.GetAtomWithIdx(idx).GetIsAromatic()

        a=Chem.Atom(symbol)
        a.SetChiralTag(chiral_tag)
        a.SetFormalCharge(formal_charge)
        a.SetIsAromatic(is_aromatic)
        a.SetHybridization(hybridization)
        a.SetNumExplicitHs(num_explicit_hs)
        out_idx = new_mol.AddAtom(a)
        node_to_idx[idx] = out_idx

    subset_possible_bonds = list(itertools.combinations(clique, 2))
    for sub in subset_possible_bonds:
        bond = mol.GetBondBetweenAtoms(*sub)
        if bond:
            ifirst = node_to_idx[sub[0]]
            isecond = node_to_idx[sub[1]]
            bond_type = bond.GetBondType()
            new_mol.AddBond(ifirst, isecond, bond_type)

    new_mol2 = new_mol.GetMol()
    new_mol2.UpdatePropertyCache()

    return new_mol2


class MolTreeNode(object):

    def __init__(self, mol, clique=[], smiles=None):
        self.mol, self.tri_mol = mol, mol
        self.graph = mol_to_nx(mol)
        self.neighbors = []
               
    def add_neighbor(self, nei_node):
        self.neighbors.append(nei_node)

def copy_edit_mol(mol):
    new_mol = Chem.RWMol(Chem.MolFromSmiles(''))
    for atom in mol.GetAtoms():
        new_atom = copy_atom(atom)
        new_mol.AddAtom(new_atom)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        bt = bond.GetBondType()
        new_mol.AddBond(a1, a2, bt)

        # added
        b_ghs = bond.GetBoolProp(KEY)
        new_bond = new_mol.GetBondBetweenAtoms(a1, a2)
        new_bond.SetBoolProp(KEY, b_ghs)
    return new_mol

#---------------------------------------------------------------------------

def atom_equal(a1, a2):
    return a1.GetSymbol() == a2.GetSymbol() and a1.GetFormalCharge() == a2.GetFormalCharge()

def bond_prop_equal(b1, b2):
    return b1.GetBoolProp(KEY) == b2.GetBoolProp(KEY)

#Bond type not considered because all aromatic (so SINGLE matches DOUBLE)
def ring_bond_equal(b1, b2, reverse=False):
    bond_prop = bond_prop_equal(b1, b2)
    b1 = (b1.GetBeginAtom(), b1.GetEndAtom())
    if reverse:
        b2 = (b2.GetEndAtom(), b2.GetBeginAtom())
    else:
        b2 = (b2.GetBeginAtom(), b2.GetEndAtom())
    return atom_equal(b1[0], b2[0]) and atom_equal(b1[1], b2[1]) and bond_prop

def enclosed_tri_clique(b1, b2, enclosed):
    if enclosed: return b1.GetBoolProp(KEY) is True and b2.GetBoolProp(KEY) is True # only allow ghost bond attachment
    else: return True # allow all attachment

def attach_mols(ctr_mol, neighbors, prev_nodes, nei_amap):
    prev_nids = [node.nid for node in prev_nodes]

    for nei_node in prev_nodes + neighbors:
        nei_id,nei_mol = nei_node.nid,nei_node.tri_mol # mol -> tri_mol
        amap = nei_amap[nei_id]
        for atom in nei_mol.GetAtoms():
            if atom.GetIdx() not in amap:
                new_atom = copy_atom(atom)
                amap[atom.GetIdx()] = ctr_mol.AddAtom(new_atom) # 2(nei) : 3(ctr)

        if nei_mol.GetNumBonds() == 0:
            nei_atom = nei_mol.GetAtomWithIdx(0)
            ctr_atom = ctr_mol.GetAtomWithIdx(amap[0])
            ctr_atom.SetAtomMapNum(nei_atom.GetAtomMapNum())
        else:
            for bond in nei_mol.GetBonds():
                a1 = amap[bond.GetBeginAtom().GetIdx()]
                a2 = amap[bond.GetEndAtom().GetIdx()] # get ctr_idx 
                # print(a1, a2)
                if ctr_mol.GetBondBetweenAtoms(a1, a2) is None:
                    ctr_mol.AddBond(a1, a2, bond.GetBondType())
                    #-----------------------------------------------
                    new_bond = ctr_mol.GetBondBetweenAtoms(a1, a2)
                    new_bond.SetBoolProp(KEY, bond.GetBoolProp(KEY))
                elif nei_id in prev_nids: #father node overrides
                    ctr_mol.RemoveBond(a1, a2)
                    ctr_mol.AddBond(a1, a2, bond.GetBondType())
                    #------------------------------------------------
                    new_bond = ctr_mol.GetBondBetweenAtoms(a1, a2)
                    new_bond.SetBoolProp(KEY, bond.GetBoolProp(KEY))
    return ctr_mol


def local_attach(ctr_mol, neighbors, prev_nodes, amap_list):
    ctr_mol = copy_edit_mol(ctr_mol)
    nei_amap = {nei.nid:{} for nei in prev_nodes + neighbors}

    for nei_id,ctr_atom,nei_atom in amap_list:
        nei_amap[nei_id][nei_atom] = ctr_atom

    ctr_mol = attach_mols(ctr_mol, neighbors, prev_nodes, nei_amap)
    return ctr_mol.GetMol()


MAX_NCAND = 6000

#This version records idx mapping between ctr_mol and nei_mol
#Keep attaching 
def enum_attach(ctr_node, nei_node, amap, singletons):
    ctr_mol = ctr_node.tri_mol
    nei_mol,nei_idx = nei_node.tri_mol,nei_node.nid # mol -> tri_mol


    att_confs = []
    black_list = [atom_idx for nei_id,atom_idx,_ in amap if nei_id in singletons]
    ctr_atoms = [atom for atom in ctr_mol.GetAtoms() if atom.GetIdx() not in black_list]
    ctr_bonds = [bond for bond in ctr_mol.GetBonds()]

    #---------------------------------------------------------------
    count_ctr_ghost = sum([bond.GetBoolProp(KEY) for bond in ctr_bonds])
    count_nei_ghost = sum([bond.GetBoolProp(KEY) for bond in nei_mol.GetBonds()])
    count_neigh_of_nei = sum([1 for nei in nei_node.neighbors if nei.nid != ctr_node.nid])
    # print(count_ctr_ghost, count_nei_ghost, count_neigh_of_nei)
    #---------------------------------------------------------------

    if nei_mol.GetNumBonds() == 0: #neighbor singleton
        nei_atom = nei_mol.GetAtomWithIdx(0)
        used_list = [atom_idx for _,atom_idx,_ in amap]
        for atom in ctr_atoms:
            if atom_equal(atom, nei_atom) and atom.GetIdx() not in used_list:
                new_amap = amap + [(nei_idx, atom.GetIdx(), 0)]
                att_confs.append( new_amap )

        if not att_confs:
            for atom in ctr_atoms:
                if atom_equal(atom, nei_atom):
                    new_amap = amap + [(nei_idx, atom.GetIdx(), 0)]
                    att_confs.append( new_amap )
   
    elif nei_mol.GetNumBonds() == 1: #neighbor is a bond
        bond = nei_mol.GetBondWithIdx(0)
        bond_val = int(bond.GetBondTypeAsDouble())
        b1,b2 = bond.GetBeginAtom(), bond.GetEndAtom()

        for atom in ctr_atoms: 
            #Optimize if atom is carbon (other atoms may change valence)
            # if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() < bond_val:
            #     continue
            if atom_equal(atom, b1):
                new_amap = amap + [(nei_idx, atom.GetIdx(), b1.GetIdx())]
                att_confs.append( new_amap )
            elif atom_equal(atom, b2):
                new_amap = amap + [(nei_idx, atom.GetIdx(), b2.GetIdx())]
                att_confs.append( new_amap )
    else:
        enclosed = True if nei_node.is_leaf and count_nei_ghost == 1 else False # enclosed leaf triangular clique

        #intersection is an bond
        if ctr_mol.GetNumBonds() > 1:
            for b1 in ctr_bonds:
                for b2 in nei_mol.GetBonds():
                    if ring_bond_equal(b1, b2) and enclosed_tri_clique(b1, b2, enclosed):
                        new_amap = amap + [(nei_idx, b1.GetBeginAtom().GetIdx(), b2.GetBeginAtom().GetIdx()), (nei_idx, b1.GetEndAtom().GetIdx(), b2.GetEndAtom().GetIdx())]
                        att_confs.append( new_amap )

                    if ring_bond_equal(b1, b2, reverse=True) and enclosed_tri_clique(b1, b2, enclosed):
                        new_amap = amap + [(nei_idx, b1.GetBeginAtom().GetIdx(), b2.GetEndAtom().GetIdx()), (nei_idx, b1.GetEndAtom().GetIdx(), b2.GetBeginAtom().GetIdx())]
                        att_confs.append( new_amap )

    return att_confs

def enum_assemble(node, neighbors, prev_nodes=[], prev_amap=[], print_out=False):
    all_attach_confs = []
    singletons = [nei_node.nid for nei_node in neighbors + prev_nodes if nei_node.graph.number_of_nodes() == 1]

    if print_out: print(len(neighbors))

    def search(cur_amap, depth):

        if len(all_attach_confs) > MAX_NCAND:
            return
        if depth == len(neighbors):
            all_attach_confs.append(cur_amap)
            return

        nei_node = neighbors[depth]

        cand_amap = enum_attach(node, nei_node, cur_amap, singletons) # mol -> tri_mol
        if print_out: print("cand_amap", cand_amap, nei_node.nid)
        true_cand_graphs = []
        candidates = []
        for i, amap in enumerate(cand_amap):
            cand_mol = local_attach(node.tri_mol, neighbors[:depth+1], prev_nodes, amap) # mol -> tri_mol
            cand_graph = mol_to_nx(cand_mol)

            # draw_mol(cand_graph, count * 1000 + i, folder="../subgraph")

            smiles = Chem.MolToSmiles(cand_mol)
            if print_out: print(smiles)

            duplicate = len([1 for G in true_cand_graphs if nx.is_isomorphic(G, cand_graph, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)]) # more candidate due to more specific
            if duplicate: continue
            true_cand_graphs.append(cand_graph)
            
            candidates.append(amap)

        if len(candidates) == 0:
            return

        for new_amap in candidates:
            search(new_amap, depth + 1)

    search(prev_amap, 0)

    cand_smiles = set()
    candidates = []
    candidates_G = []
    for i, amap in enumerate(all_attach_confs):
        cand_mol = local_attach(node.tri_mol, neighbors, prev_nodes, amap)
        # cand_mol = Chem.MolFromSmiles(Chem.MolToSmiles(cand_mol))
        cand_G = mol_to_nx(cand_mol)
        smiles = Chem.MolToSmiles(cand_mol)
        cand_smiles.add(smiles)
        try: Chem.Kekulize(cand_mol)
        except: pass
        candidates.append( (smiles,cand_mol,amap) )
        candidates_G.append( (smiles,cand_G,amap) )

    return candidates, candidates_G


def dfs_all_assemble(cur_graph, global_amap, fa_amap, cur_node, fa_node, possible_attachment):

    fa_nid = fa_node.nid if fa_node is not None else -1
    prev_nodes = [fa_node] if fa_node is not None else []

    children = [nei for nei in cur_node.neighbors if nei.nid != fa_nid]

    neighbors = [nei for nei in children if nei.graph.number_of_nodes() > 1]
    neighbors = sorted(neighbors, key=lambda x:x.graph.number_of_nodes(), reverse=True)
    singletons = [nei for nei in children if nei.graph.number_of_nodes() == 1]
    neighbors = singletons + neighbors

    cur_amap = [(fa_nid,a2,a1) for nid,a1,a2 in fa_amap if nid == cur_node.nid] # check if there is any atommap has been occupied previously when attaching with parent
    
    cands, cands_G = enum_assemble(cur_node, neighbors, prev_nodes, cur_amap)
    
    # add all possible attachments
    possible_attachment[cur_node.nid] = len(cands)
    
    if (not cands) and (not cands_G):
        return
    
    cand_smiles, cand_Gs, cand_amap = zip(*cands_G)
        
    if not cur_node.cand_idx:
        label_idx = np.random.randint(0, len(cands_G))
        label_amap = cand_amap[label_idx]
    else:
        label_amap = cand_amap[cur_node.cand_idx]

    for nei_id,ctr_atom,nei_atom in label_amap:
        if nei_id == fa_nid:
            continue
        global_amap[nei_id][nei_atom] = global_amap[cur_node.nid][ctr_atom]

    cur_graph = attach_graphs(cur_graph, children, [], global_amap) #father is already attached

    for nei_node in children:
        if not nei_node.is_leaf:
            dfs_all_assemble(cur_graph, global_amap, label_amap, nei_node, cur_node, possible_attachment)

def node_labelling(mol, cliques, molTreeEdges, triangulated_graph):

    list_of_nodes = []
    for i, clique in enumerate(cliques):
        if isinstance(clique, tuple):
            c = list(clique)
            m = MolTreeNode(mol, c)
        list_of_nodes.append(m)

    for x,y in molTreeEdges:
        list_of_nodes[x].add_neighbor(list_of_nodes[y])
        list_of_nodes[y].add_neighbor(list_of_nodes[x])

    for i,node in enumerate(list_of_nodes):
        node.nid = i + 1
        if len(node.neighbors) > 1: #Leaf node mol is not marked
            set_atommap(node.mol, node.nid)
            set_atommap(node.tri_mol, node.nid)
            set_atommap_graph(node.graph, node.nid)
        node.is_leaf = (len(node.neighbors) == 1)

    for node in list_of_nodes:
        node.recover_G(triangulated_graph.copy())

    root_idx = 0
    root = list_of_nodes[root_idx]
    cur_graph = root.graph.copy()
    global_amap = [{}] + [{} for node in list_of_nodes] # nid starts with 1, thus the first dict is not used, just as offset
    global_amap[root_idx+1] = {node:node for node in cur_graph.nodes()}

    return root, cur_graph, global_amap


def remove_edges_reset_idx(cur_graph):
    final_graph = cur_graph.copy()
    for a1, a2, data in cur_graph.edges(data=True):
        if data.get(KEY): final_graph.remove_edge(a1, a2)

    set_atommap_graph(final_graph)

    return final_graph
