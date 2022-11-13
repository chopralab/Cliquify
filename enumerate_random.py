
from rdkit import Chem
import itertools
import networkx as nx
# from Cliquify.utils import *
from utils import *
import networkx.algorithms.isomorphism as iso
# from Cliquify.debug_script import *
from debug_script import *
from rdkit import RDLogger
import numpy as np
# from Cliquify.tree_decomposition2 import tree_decomp
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
        if mol and smiles: 
            self.smiles = smiles
            self.mol = mol

            self.tri_mol, self.graph = get_triangulated_graph(self.mol)
            self.clique = [x for x in clique] #copy
            self.neighbors = []
        else:
            self.smiles = get_fragments(mol , clique)
            self.smarts = get_smarts_fragments(mol , clique)

            self.mol = get_mol2(mol, clique) # use mol object directly
            # self.mol = Chem.MolFromSmiles(self.smiles) # use smiles fragment string
            # self.mol = Chem.MolFromSmarts(self.smarts) # use smarts

            self.tri_mol, self.graph = get_triangulated_graph(self.mol)

            self.clique = [x for x in clique] #copy
            self.neighbors = []
        
    def add_neighbor(self, nei_node):
        self.neighbors.append(nei_node)

    def recover_G(self, original_graph):
        clique = []
        clique.append(self.clique)

        if not self.is_leaf:
            for cidx in self.clique:
                original_graph.nodes[cidx]["map_num"] = self.nid
        
        # all_edges = list(self.graph.edges())
        for nei_node in self.neighbors:
            clique.append(nei_node.clique)
            # all_edges.extend(list(self.graph.edges()))
            if nei_node.is_leaf: #Leaf node, no need to mark 
                continue
            for cidx in nei_node.clique:
                #allow singleton node override the atom mapping
                if cidx not in self.clique or len(nei_node.clique) == 1: # neighboring node atom will only override non ctr clique atom
                    node = original_graph.nodes[cidx]
                    node["map_num"] = nei_node.nid
                    # print(cidx, nei_node.nid, 'current', self.nid)
        
        valid_bonds = []
        for cliq in clique:
            subset_possible_bonds = list(itertools.combinations(cliq, 2))
            for bond in subset_possible_bonds:
                if original_graph.has_edge(*bond): valid_bonds.append(set(bond))

        clique = list(set([node for cliq in clique for node in cliq]))
        temp_graph = original_graph.subgraph(clique).copy()
        self.label_G = original_graph.subgraph(clique).copy()

        for edge in list(temp_graph.edges()):
            if set(edge) not in valid_bonds:
                self.label_G.remove_edge(*edge)

        # if self.nid == 1:
        #     draw_mol(self.label_G, self.nid + 10_000, ['map_num', 'bond_type', 'color'])
        
        return self.label_G
            
def get_smiles(mol):
    return Chem.MolToSmiles(mol, kekuleSmiles=True)

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


def remapping(mol):
    triangulated_mol = Chem.RWMol(Chem.MolFromSmiles(''))
    for atom in mol.GetAtoms():
        new_atom = copy_atom(atom)
        triangulated_mol.AddAtom(new_atom)

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()

        try:
            if not bond.GetBoolProp(KEY):
                triangulated_mol.AddBond(a1, a2, order=bond.GetBondType())
            else:
                triangulated_mol.AddBond(a1, a2, order=Chem.BondType.DOUBLE)
        except:
            triangulated_mol.AddBond(a1, a2, order=bond.GetBondType())

    return triangulated_mol.GetMol()

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


def local_attach2(ctr_graph, neighbors, prev_nodes, amap_list):
    inside_graph = ctr_graph.copy()
    nei_amap = {nei.nid:{} for nei in prev_nodes + neighbors}

    for nei_id,ctr_atom,nei_atom in amap_list:
        nei_amap[nei_id][nei_atom] = ctr_atom

    inside_graph = attach_graphs(inside_graph, neighbors, prev_nodes, nei_amap)
    return inside_graph


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
        # if not nei_node.is_leaf and False:

        #intersection is an atom
        if count_nei_ghost == 2 and count_neigh_of_nei >= 2:
            for a1 in ctr_atoms:
                for a2 in nei_mol.GetAtoms():
                    if atom_equal(a1, a2):
                        #Optimize if atom is carbon (other atoms may change valence)
                        # if a1.GetAtomicNum() == 6 and a1.GetTotalNumHs() + a2.GetTotalNumHs() < 4:
                        #     continue
                        new_amap = amap + [(nei_idx, a1.GetIdx(), a2.GetIdx())]
                        att_confs.append( new_amap )

        #intersection is an bond
        if ctr_mol.GetNumBonds() > 1:
            for b1 in ctr_bonds:
                for b2 in nei_mol.GetBonds():
                    if ring_bond_equal(b1, b2):
                        new_amap = amap + [(nei_idx, b1.GetBeginAtom().GetIdx(), b2.GetBeginAtom().GetIdx()), (nei_idx, b1.GetEndAtom().GetIdx(), b2.GetEndAtom().GetIdx())]
                        att_confs.append( new_amap )

                    if ring_bond_equal(b1, b2, reverse=True):
                        new_amap = amap + [(nei_idx, b1.GetBeginAtom().GetIdx(), b2.GetEndAtom().GetIdx()), (nei_idx, b1.GetEndAtom().GetIdx(), b2.GetBeginAtom().GetIdx())]
                        att_confs.append( new_amap )

    return att_confs

def enum_assemble(node, neighbors, prev_nodes=[], prev_amap=[], print_out=False):
    all_attach_confs = []
    singletons = [nei_node.nid for nei_node in neighbors + prev_nodes if nei_node.graph.number_of_nodes() == 1]

    if print_out: print(len(neighbors))

    # count = 0
    def search(cur_amap, depth):
        # global count
        # count += 1

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
            # print(amap)

            cand_mol = local_attach(node.tri_mol, neighbors[:depth+1], prev_nodes, amap) # mol -> tri_mol
            cand_graph = mol_to_nx(cand_mol)

            # draw_mol(cand_graph, count * 1000 + i, folder="../subgraph")

            smiles = Chem.MolToSmiles(cand_mol)
            if print_out: print(smiles)

            # duplicate = len([1 for G in true_cand_graphs if nx.is_isomorphic(G, cand_graph)]) # less candidate due to less specific
            duplicate = len([1 for G in true_cand_graphs if nx.is_isomorphic(G, cand_graph, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)]) # more candidate due to more specific
            if duplicate: continue
            true_cand_graphs.append(cand_graph)
            
            candidates.append(amap)

        if len(candidates) == 0:
            return

        for new_amap in candidates:
            search(new_amap, depth + 1)

    search(prev_amap, 0)

    # if print_out: print("cands_id", cands_id), print(len(all_attach_confs))

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


def dfs_random_assemble(cur_graph, global_amap, fa_amap, cur_node, fa_node, print_out=False):
    global count

    fa_nid = fa_node.nid if fa_node is not None else -1
    prev_nodes = [fa_node] if fa_node is not None else []

    children = [nei for nei in cur_node.neighbors if nei.nid != fa_nid]

    neighbors = [nei for nei in children if nei.graph.number_of_nodes() > 1]
    neighbors = sorted(neighbors, key=lambda x:x.graph.number_of_nodes(), reverse=True)
    singletons = [nei for nei in children if nei.graph.number_of_nodes() == 1]
    neighbors = singletons + neighbors

    cur_amap = [(fa_nid,a2,a1) for nid,a1,a2 in fa_amap if nid == cur_node.nid] # check if there is any atommap has been occupied previously when attaching with parent
    
    cands, cands_G = enum_assemble(cur_node, neighbors, prev_nodes, cur_amap)
    
    if (not cands) and (not cands_G):
        return
    
    cand_smiles, cand_Gs, cand_amap = zip(*cands_G)

    # heuristics on cand_Gs can be done here

    # random sample between plausible labels
    label_idx = np.random.randint(0, len(cands_G))

    label_amap = cand_amap[label_idx]

    for nei_id,ctr_atom,nei_atom in label_amap:
        if nei_id == fa_nid:
            continue
        global_amap[nei_id][nei_atom] = global_amap[cur_node.nid][ctr_atom]

    cur_graph = attach_graphs(cur_graph, children, [], global_amap) #father is already attached

    for nei_node in children:
        if not nei_node.is_leaf:
            dfs_random_assemble(cur_graph, global_amap, label_amap, nei_node, cur_node)

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
    
def reconstruction_evaluation(chosen_smiles, final_graph):
    cur_mol = nx_to_mol(mol_to_nx(nx_to_mol(final_graph)))

    chosen_graph = mol_to_nx(Chem.MolFromSmiles(chosen_smiles))
    graph_match = nx.is_isomorphic(final_graph, chosen_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)


    gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False)
    dec_smiles = Chem.MolToSmiles(cur_mol, isomericSmiles=False)

    return gold_smiles, dec_smiles, graph_match
    
def add_ghost_edges(G, ghost_edges):

    for edge in ghost_edges:
        G.add_edge(edge[0],
                edge[1],
                bond_type=Chem.BondType.SINGLE,
                ghost=True,
                color='r',
                )
    return G
