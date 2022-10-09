from rdkit import Chem
import itertools
import networkx as nx
from utils import *
import networkx.algorithms.isomorphism as iso
from debug_script import *
from rdkit import RDLogger
from tree_decomposition2 import tree_decomp
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)
import copy

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
    try:
        return Chem.MolFragmentToSmiles(mol, atoms, kekuleSmiles=True)
    except:
        return Chem.MolFragmentToSmiles(mol, atoms)

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

    def __init__(self, mol, clique=[]):
        self.smiles = get_fragments(mol , clique)
        # self.mol = get_mol(self.smiles)
        # print()
        # print(Chem.MolToSmiles(self.mol))

        self.mol = get_mol2(mol, clique)
        # print(Chem.MolToSmiles(self.mol))
        # print("\n\n")

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

    # print()
    # print(nei_amap)
    # print("nei_id : {nei_atom : ctr_atom }")

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


MAX_NCAND = 9000

#This version records idx mapping between ctr_mol and nei_mol
#Keep attaching 
def enum_attach(ctr_mol, nei_node, amap, singletons):
    nei_mol,nei_idx = nei_node.tri_mol,nei_node.nid # mol -> tri_mol

    att_confs = []
    black_list = [atom_idx for nei_id,atom_idx,_ in amap if nei_id in singletons]
    ctr_atoms = [atom for atom in ctr_mol.GetAtoms() if atom.GetIdx() not in black_list]
    ctr_bonds = [bond for bond in ctr_mol.GetBonds()]

    if nei_mol.GetNumBonds() == 0: #neighbor singleton
        nei_atom = nei_mol.GetAtomWithIdx(0)
        used_list = [atom_idx for _,atom_idx,_ in amap]
        for atom in ctr_atoms:
            # if nei_idx == 29:
            #     print(atom.GetIdx(), "used_list", used_list)
            #     print(atom_equal(atom, nei_atom), atom.GetIdx() not in used_list)
            #     print()
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
        #intersection is an atom
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

import json
class Graph:
    def __init__(self, graph) -> None:
        self.graph = graph
    def __eq__(self, other):
        if type(other) is type(self):
            GM = iso.GraphMatcher(self.graph, other.graph, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)
            return GM.is_isomorphic()
        else:
            return False
    def __hash__(self) -> int:
        # return hash(nx.weisfeiler_lehman_graph_hash(self.graph,
        #     edge_attr='color',
        #     node_attr='map_num'))
        # return hash(nx.weisfeiler_lehman_graph_hash(self.graph))

        edges_data = []
        for n1, n2, data in self.graph.edges.data():
            data = json.dumps(data)
            # edges_data.append((n1, n2, data))
            edges_data.append((n1, n2))

        nodes_data = []
        for n, data in self.graph.nodes.data():
            data = json.dumps(data)
            nodes_data.append((n))

        return hash((tuple(edges_data), tuple(nodes_data)))

count_i = 0
def enum_assemble_singleton_tri(node, prev_nodes, neighbors, prev_amap): # neighbors exclude prev_nodes

    # neighbors = node.neighbors
    possible_cands_G = []
    # cur_graph = neighbors[0].graph.copy()
    cur_graph = node.graph.copy()
    # print(neighbors[0].clique)
    # draw_mol(node.graph, 0, ["map_num", "bond_type", "color"])
    # draw_mol(cur_graph, 1, ["map_num", "bond_type", "color"])
    # nei_print = [nei.clique for nei in node.neighbors]
    nei_print = [nei.clique for nei in neighbors]
    print(nei_print)
    # neighbors = neighbors[:2]

    def search(cur_graph, depth):
        # print('current depth', depth)

        if depth == len(neighbors):
            possible_cands_G.append(cur_graph)
            return

        nei_node = neighbors[depth]

        # print('depth', depth)
        try: prev_node = neighbors[depth - 1]
        except: prev_node = []
        # cands_G = enum_attach_single_bond(cur_graph, nei_node, prev_node)
        cands_G = enum_attach_single_bond(cur_graph, nei_node)

        if node.tri_mol.GetNumAtoms() == 1 and depth >= 2: # atom as singleton
            cands_G.extend(enum_attach_double_bond(cur_graph, nei_node.graph))

        if node.tri_mol.GetNumAtoms() == 2 and depth >= 3: # bond as singleton
            cands_G.extend(enum_attach_double_bond(cur_graph, nei_node.graph))


        filtered_graph = []
        for i, cand_G in enumerate(cands_G): # candidate tree pruning
            node_adj = [True for node in cand_G.nodes() if len(list(cand_G.neighbors(node))) > 4]
            if node_adj: continue

            node_adj2 = [True for node in cand_G.nodes() if len(list(cand_G.neighbors(node))) == 2 and depth >= 3]
            if len(node_adj2) > 2: continue

            cand_ghost_edge = sum(nx.get_edge_attributes(cand_G, "ghost").values()) # sum of ghost edge
            cur_ghost_edge = sum(nx.get_edge_attributes(cur_graph, "ghost").values()) # sum of ghost edge
            nei_ghost_edge = sum(nx.get_edge_attributes(nei_node.graph, "ghost").values()) # sum of ghost edge
            if (cand_ghost_edge < (cur_ghost_edge + nei_ghost_edge - 1)) or (cand_ghost_edge > (cur_ghost_edge + nei_ghost_edge)): continue

            # GM = iso.GraphMatcher(node.label_G, cand_G)
            # if not GM.subgraph_is_isomorphic(): continue
                
            filtered_graph.append(cand_G)
        
        if len(filtered_graph) == 0:
            return
        
        for i, cand_G in enumerate(filtered_graph):
            global count_i
            # if depth == 0:
            #     draw_mol(cand_G, 1000 + count_i, ["map_num", "bond_type", "color"])
            #     count_i += 1
            # if depth == 1:
            #     draw_mol(cand_G, 2000 + count_i, ["map_num", "bond_type", "color"])
            #     count_i += 1
            search(cand_G, depth + 1)

    search(cur_graph, 0)

    print("candidates:", len(possible_cands_G))
    print("\n")
    # raise

    temp_label = node.label_G.copy()
    temp_label.remove_edge(5, 9)
    # temp_label.remove_edge(2, 12)
    cand_graphs = set()
    # draw_mol(temp_label, 100000, ['map_num', 'bond_type', 'color'])
    for i, cand_G in enumerate(possible_cands_G):
        if Graph(cand_G) in cand_graphs: continue
        cand_graphs.add(Graph(cand_G))

        # GM = iso.GraphMatcher(temp_label, cand_G, edge_match=ring_edge_equal_iso)
        GM = iso.GraphMatcher(temp_label, cand_G, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        graph_match = GM.is_isomorphic()
        if graph_match:
            print(200000 + i, "graph_match", graph_match)
            print()
            draw_mol(cand_G, 200000 + i, ['map_num', 'bond_type', 'color'])

    raise

    cand_graphs = set()
    # draw_mol(node.label_G, 100000, ['map_num', 'bond_type', 'color'])
    for i, cand_G in enumerate(possible_cands_G):
        try: mol = nx_to_mol(cand_G)
        except: continue

        if Graph(cand_G) in cand_graphs: continue
        cand_graphs.add(Graph(cand_G))

        cand_smiles = Chem.MolToSmiles(mol)
        # graph_match = nx.is_isomorphic(node.label_G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)
        # GM = iso.GraphMatcher(node.label_G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)
        GM = iso.GraphMatcher(node.label_G, cand_G, edge_match=ring_edge_equal_iso)
        graph_match = GM.is_isomorphic()
        # smiles_match = node.label == cand_smiles
        # if smiles_match and graph_match:
        if graph_match:
            print(cand_smiles)
            # print("smiles_match", smiles_match)
            print("graph_match", graph_match)
            print()
            draw_mol(cand_G, 100000 + i, ['map_num', 'bond_type', 'color'])
    raise

    return possible_cands

def enum_assemble(node, neighbors, prev_nodes=[], prev_amap=[], print_out=False):
    all_attach_confs = []
    singletons = [nei_node.nid for nei_node in neighbors + prev_nodes if nei_node.graph.number_of_nodes() == 1]

    # #------------------------ GET TO THIS AFTER ALL GRAPH CONV-----------------------------
    # if node.graph.number_of_nodes() <= 2 and len(node.neighbors) >= 3 and sum(tri_mol_count) >= 3 * len(tri_mol_count):
    #     candidates = enum_assemble_singleton_tri(node, prev_nodes, neighbors, prev_amap)
    #     return candidates

    if print_out: print(len(neighbors))
    # cands_id = []
    def search(cur_amap, depth):
        if len(all_attach_confs) > MAX_NCAND:
            return
        if depth == len(neighbors):
            all_attach_confs.append(cur_amap)
            return

        nei_node = neighbors[depth]

        cand_amap = enum_attach(node.tri_mol, nei_node, cur_amap, singletons) # mol -> tri_mol
        if print_out: print("cand_amap", cand_amap, nei_node.nid)
        cand_smiles = set()
        cand_graphs = set()
        true_cand_graphs = []
        candidates = []
        for i, amap in enumerate(cand_amap):
            cand_mol = local_attach(node.tri_mol, neighbors[:depth+1], prev_nodes, amap) # mol -> tri_mol
            cand_graph = mol_to_nx(cand_mol)
            smiles = Chem.MolToSmiles(cand_mol)
            if print_out: print(smiles)

            # cand_graph = local_attach2(node.graph, neighbors[:depth+1], prev_nodes, amap) # graph alternative

            # smiles = Chem.MolToSmiles(cand_mol)
            # if smiles in cand_smiles:
            #     continue
            # cand_smiles.add(smiles)

            # if Graph(cand_graph) in cand_graphs:
            #     continue
            # cand_graphs.add(Graph(cand_graph))

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
        # if smiles in cand_smiles:
        #     continue
        cand_smiles.add(smiles)
        try: Chem.Kekulize(cand_mol)
        except: pass
        candidates.append( (smiles,cand_mol,amap) )
        candidates_G.append( (smiles,cand_G,amap) )

    return candidates, candidates_G

def enum_assemble_singleton_tri2(cur_graph, cur_node, neighbors, prev_nodes, temp_global_amap): # neighbors exclude prev_nodes
    # cur_graph and cur_node.graph is a different thing, cur_graph indicates the dfs assembled up to this point
    nid_graph = {nei.nid: nei.graph.copy() for nei in neighbors}
    
    if prev_nodes: cur_nid, fa_nid = cur_node.nid, prev_nodes[0].nid
    else: cur_nid, fa_nid = cur_node.nid, 0

    cur_clique = list(temp_global_amap[cur_nid].values())
    fa_clique = list(temp_global_amap[fa_nid].values())
    label_starting_clique = set(cur_clique + fa_clique)

    starting_graph = cur_graph.subgraph(label_starting_clique).copy()
    possible_cands_G_set = set()
    possible_cands_G = []
    all_attach_confs_set = set()
    all_attach_confs = []

    draw_mol(cur_graph, 5999, ["map_num", "bond_type", "color"], folder="extension")
    draw_mol(starting_graph, 6000, ["map_num", "bond_type", "color"], folder="extension")
    # # draw_mol(cur_node.label_G, 6001, ["map_num", "bond_type", "color"], folder="extension")
    draw_mol(cur_node.graph, 6001, ["map_num", "bond_type", "color"], folder="extension")
    # print('cur_node.nid', cur_node.nid)
    print('prev_nodes clique', prev_nodes[0].clique)
    print('temp_global_amap', temp_global_amap)

    print()
    for nei in neighbors:
        print(nei.clique)

    def search(starting_graph, cur_amap, temp_global_amap, depth):

        if depth == len(neighbors):
            # if Graph(starting_graph) not in possible_cands_G_set and tuple(cur_amap) not in all_attach_confs_set:
            #     possible_cands_G_set.add(Graph(starting_graph))
            #     possible_cands_G.append(starting_graph)

            #     all_attach_confs_set.add(tuple(cur_amap))
            #     all_attach_confs.append(cur_amap)

            possible_cands_G.append(starting_graph)
            all_attach_confs.append(cur_amap)
            return

        nei_node = neighbors[depth]
        
        cands_G, cands_G_amap = [], []
        # cands_G, cands_G_amap = enum_attach_double_bond2(starting_graph, nei_node, copy.deepcopy(temp_global_amap))
        # cands2_G, cands2_G_amap = enum_attach_single_bond(starting_graph, nei_node, copy.deepcopy(temp_global_amap))
        try: cands_G, cands_G_amap = enum_attach_double_bond2(starting_graph, nei_node, copy.deepcopy(temp_global_amap))
        except: pass
        try: cands2_G, cands2_G_amap = enum_attach_single_bond(starting_graph, nei_node, copy.deepcopy(temp_global_amap))
        except: return
        cands_G.extend(cands2_G), cands_G_amap.extend(cands2_G_amap)
    
        # print('depth', depth, temp_global_amap)
        # # draw_mol(nei_node.graph, 6000 + nei_node.nid, ["map_num", "bond_type", "color"], folder="extension")
        # print(len(cands_G_amap), 'cands produced')

        candidate_graphs = [starting_graph.copy() for _ in cands_G_amap]
        temp_global_amaps = [copy.deepcopy(temp_global_amap) for _ in cands_G_amap]


        match_new_amaps = []
        match_new_temp_global_amaps = []
        match_new_cand_graphs = []

        mismatch_new_amaps = []
        mismatch_new_temp_global_amaps = []
        mismatch_new_cand_graphs = []
        
        possible_cands_subgraph = []

        for i, label_amap in enumerate(cands_G_amap):
            cand_graph = candidate_graphs[i]
            temp_global_amap = temp_global_amaps[i]

            # -----------------------subgraph pruning is done here-----------------------------
            cand_G = cands_G[i]
            node_adj = [True for node in cand_G.nodes() if len(list(cand_G.neighbors(node))) > 4]
            if node_adj: continue

            node_adj2 = [True for node in cand_G.nodes() if len(list(cand_G.neighbors(node))) == 2 and depth >= 3]
            if len(node_adj2) > 2: continue

            cand_ghost_edge = sum(nx.get_edge_attributes(cand_G, "ghost").values()) # sum of ghost edge
            cur_ghost_edge = sum(nx.get_edge_attributes(starting_graph, "ghost").values()) # sum of ghost edge
            nei_ghost_edge = sum(nx.get_edge_attributes(nei_node.graph, "ghost").values()) # sum of ghost edge
            if (cand_ghost_edge < (cur_ghost_edge + nei_ghost_edge - 1)) or (cand_ghost_edge > (cur_ghost_edge + nei_ghost_edge)): continue

            if Graph(cand_G) in possible_cands_G_set: continue
            duplicate = len([1 for G in possible_cands_subgraph if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)]) # more candidate due to more specific
            if duplicate: continue
            possible_cands_subgraph.append(cand_G)
            # -------------------------------------------------------------------

            for nei_id, ctr_atom,nei_atom, ctr_id in label_amap: # nei nid, atom idx, atom idx
                if nei_id == fa_nid: continue

                try: temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]
                except:
                    # add artificial nodes
                    temp_global_amap[ctr_id][ctr_atom] = len(cand_graph.nodes)
                    temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]

                    node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
                    cand_graph.add_node(len(cand_graph.nodes), **node_attr)
            
            cand_graph = attach_graphs(cand_graph, [nei_node], prev_nodes, temp_global_amap, print_out=True) #father is already attached
            

            new_amap = cur_amap + [tuple(label_amap)]
            new_temp_global_amap = copy.deepcopy(temp_global_amap)
            new_cand_graph = cand_graph.copy()

            GM = iso.GraphMatcher(cur_node.label_G, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
            GM2 = iso.GraphMatcher(starting_graph, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
            GM3 = iso.GraphMatcher(cur_node.label_G, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
            if GM.subgraph_is_isomorphic() and not GM2.is_isomorphic() or GM3.is_isomorphic():
                print((depth + 1) * 100 + i, True, label_amap)

                match_new_amaps.append(new_amap)
                match_new_temp_global_amaps.append(new_temp_global_amap)
                match_new_cand_graphs.append(new_cand_graph)
            else:
                mismatch_new_amaps.append(new_amap)
                mismatch_new_temp_global_amaps.append(new_temp_global_amap)
                mismatch_new_cand_graphs.append(new_cand_graph)

        # if len(match_new_amaps) == 0:
        #     draw_mol(nei_node.graph, 7000 + depth * 100 +  i , ["map_num", "bond_type", "color"], folder="extension")
        #     for i, cand_G in enumerate(mismatch_new_cand_graphs):
        #         draw_mol(cand_G, 7000 + depth * 100 +  i + 1 , ["map_num", "bond_type", "color"], folder="extension")

        if len(match_new_amaps) > 0:
            for i, new_amap in enumerate(match_new_amaps):
                search(match_new_cand_graphs[i], new_amap, match_new_temp_global_amaps[i], depth + 1)
        else:
            for i, new_amap in enumerate(mismatch_new_amaps):
                search(mismatch_new_cand_graphs[i], new_amap, mismatch_new_temp_global_amaps[i], depth + 1)

    search(starting_graph, [], temp_global_amap, 0)

    # for confs in list(all_attach_confs):
    #     print(confs)

    # for i, cand_G in enumerate(possible_cands_G):
    #     GM = iso.GraphMatcher(cur_node.label_G, cand_G, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
    #     graph_match = GM.is_isomorphic()
    #     if graph_match:
    #         print(i, graph_match)
    #         draw_mol(cand_G, 6100, ["map_num", "bond_type", "color"], folder="extension")
    #         draw_mol(cur_node.label_G, 6101, ["map_num", "bond_type", "color"], folder="extension")

    print("num of cands inside:", len(all_attach_confs))
    # raise

    return possible_cands_G, all_attach_confs

count = 0
#Only used for debugging purpose
def dfs_assemble(cur_graph, global_amap, fa_amap, cur_node, fa_node):
    global count

    fa_nid = fa_node.nid if fa_node is not None else -1
    prev_nodes = [fa_node] if fa_node is not None else []

    children = [nei for nei in cur_node.neighbors if nei.nid != fa_nid]

    neighbors = [nei for nei in children if nei.graph.number_of_nodes() > 1]
    neighbors = sorted(neighbors, key=lambda x:x.graph.number_of_nodes(), reverse=True)
    singletons = [nei for nei in children if nei.graph.number_of_nodes() == 1]
    neighbors = singletons + neighbors

    cur_amap = [(fa_nid,a2,a1) for nid,a1,a2 in fa_amap if nid == cur_node.nid] # check if there is any atommap has been occupied previously when attaching with parent
    
    tri_mol_count = [nei.tri_mol.GetNumAtoms() for nei in cur_node.neighbors]
    # if cur_node.graph.number_of_nodes() <= 2 and len(cur_node.neighbors) >= 3 and sum(tri_mol_count) == 3 * len(tri_mol_count):
    if cur_node.graph.number_of_nodes() <= 2 and len(cur_node.neighbors) >= 3 and sum(tri_mol_count) == 3 * len(tri_mol_count) and not min(tri_mol_count) < 3:
        # print('in1', cur_node.nid)
        print("tri_mol_count", tri_mol_count)
        cands, cands_G = None, [(nx_to_mol(cur_node.label_G), cur_node.label_G, [])]
    else:
        # print('in2', cur_node.nid)
        cands, cands_G = enum_assemble(cur_node, neighbors, prev_nodes, cur_amap)

    if (not cands) and (not cands_G):

        # cands, cands_G = enum_assemble(cur_node, neighbors, prev_nodes, cur_amap, print_out=True)
        print('exit', cur_node.nid)
        print('nei', [nei.nid for nei in cur_node.neighbors])
        print('children', [nei.nid for nei in neighbors])
        return

    # cand_smiles, cand_mols, cand_amap = zip(*cands)
    # label_idx = cand_smiles.index(cur_node.label)
    # label_amap = cand_amap[label_idx]

    label_idx = None
    cand_smiles, cand_Gs, cand_amap = zip(*cands_G)
    for i, cand_G in enumerate(cand_Gs):
        GM = iso.GraphMatcher(cur_node.label_G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)
        if GM.is_isomorphic():
            label_idx = i
            break

    print('num of cands:', len(cand_Gs), "cur_id", cur_node.nid, "cur_clique", cur_node.clique, "nei", [nei.clique for nei in cur_node.neighbors])

    try:
        label_amap = cand_amap[label_idx]
        # if cur_node.nid == 25: raise
    except:
        draw_mol(cur_graph, 1, ["map_num", "bond_type", "color"], folder="next_cand", label='cur_graph')
        
        draw_mol(cur_node.label_G, 999, ["map_num", "bond_type", "color"], folder="next_cand", label="label_G")
        for i, cand_G in enumerate(cand_Gs):
            draw_mol(cand_G, 1000 + i, ["map_num", "bond_type", "color"], folder="next_cand")
            # GM = iso.GraphMatcher(cur_node.label_G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)
            GM = iso.GraphMatcher(cur_node.label_G, cand_G, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso) # exclude hss calculation
            if GM.is_isomorphic():
                print(GM.is_isomorphic(), 1000 + i)
                label_idx = i
                break
        # raise

    try:
        label_amap = cand_amap[label_idx]
    except:
        print('cur_node.nid2', cur_node.nid)
        draw_mol(cur_graph, 6000 + count, ["map_num", "bond_type", "color"], folder="cur_graph")

        draw_mol(cur_node.label_G, 999, ["map_num", "bond_type", "color"], folder="next_cand", label="label_G")
        gold_smiles = Chem.MolToSmiles(nx_to_mol(cur_node.label_G), isomericSmiles=False)
        for i, smiles in enumerate(cand_smiles):
            if smiles == gold_smiles:
                label_idx = i
                print('here1', i)
                break
            print(smiles)
            print(gold_smiles)
            print()
                
        raise
    
    try:
        label_amap = cand_amap[label_idx]
    except:
        for i, cand_G in enumerate(cand_Gs):
            GM = iso.GraphMatcher(cur_node.label_G, cand_G, node_match=node_equal_iso)
            if GM.is_isomorphic():
                print('here2', i)
                label_idx = i
                break

        label_amap = cand_amap[label_idx]



    for nei_id,ctr_atom,nei_atom in label_amap:
        if nei_id == fa_nid:
            continue
        global_amap[nei_id][nei_atom] = global_amap[cur_node.nid][ctr_atom]

    if len(label_amap) == 0:

        # debug_singleton_3nei(cur_graph, cur_node, neighbors, prev_nodes, fa_nid, global_amap)
        # raise

        nid_graph = {child.nid: child.graph.copy() for child in children}
        # temp_global_amap = [amap if i == cur_node.nid or i == fa_nid else {} for i, amap in enumerate(global_amap) ] # filter all those amap which are not part of cur_node and prev_node
        temp_global_amap = [amap for i, amap in enumerate(global_amap)] # filter all those amap which are not part of cur_node and prev_node
        cands_G_label, all_attach_confs = enum_assemble_singleton_tri2(cur_graph, cur_node, neighbors, prev_nodes, temp_global_amap)

        chosen_i = None
        for i, total_label_amap in enumerate(all_attach_confs):
            GM = iso.GraphMatcher(cur_node.label_G, cands_G_label[i], node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
            if GM.is_isomorphic():
                print(i, True)
                chosen_i  = i

        total_label_amap = all_attach_confs[chosen_i]
        concat_label_amap = [label_amap for label_amap_list in total_label_amap for label_amap in label_amap_list]

        # print('total_label_amap', total_label_amap)
        # print('concat_label_amap', concat_label_amap)
        # print("global_amap", global_amap)
        # global_amap[children[0].nid] = {node:node for node in children[0].graph.nodes()} # allow global amap for loop to feed value in from ctr_node/cur_node
        for nei_id,ctr_atom,nei_atom, ctr_id in concat_label_amap: # nei nid, atom idx, atom idx
            if nei_id == fa_nid: continue

            try: global_amap[nei_id][nei_atom] = global_amap[ctr_id][ctr_atom]
            except:
                # add artificial nodes
                global_amap[ctr_id][ctr_atom] = len(cur_graph.nodes)
                global_amap[nei_id][nei_atom] = global_amap[ctr_id][ctr_atom]

                node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
                cur_graph.add_node(len(cur_graph.nodes), **node_attr)
        
        # print("global_amap", global_amap)
        cur_graph = attach_graphs(cur_graph, children, prev_nodes, global_amap)
        # draw_mol(cur_graph, 6000 + count, ["map_num", "bond_type", "color"], folder="cur_graph")

        for i, nei_node in enumerate(children):
            if not nei_node.is_leaf:
                # return
                if i == len(children) - 1:
                    print(total_label_amap[i])
                    label_amap = [amap[:3] for amap in total_label_amap[i]]
                else:
                    label_amap = [amap[:3] for amap in total_label_amap[i] if amap[3] == cur_node.nid]
                dfs_assemble(cur_graph, global_amap, label_amap, nei_node, cur_node)

    else:
        # cur_mol = attach_mols(cur_mol, children, [], global_amap) #father is already attached
        cur_graph = attach_graphs(cur_graph, children, [], global_amap) #father is already attached

        for nei_node in children:
            if not nei_node.is_leaf:
                dfs_assemble(cur_graph, global_amap, label_amap, nei_node, cur_node)

def add_ghost_edges(G, ghost_edges):

    for edge in ghost_edges:
        G.add_edge(edge[0],
                edge[1],
                bond_type=Chem.BondType.SINGLE,
                ghost=True,
                color='r',
                )
    return G

def gen_mol1():
    mol = Chem.MolFromSmiles("C1CC2CCC3CCCC4CCC(C1)C2C34")
    for bond in mol.GetBonds(): bond.SetBoolProp(KEY, False)
    # for bond in mol.GetBonds(): print(mol.GetBoolProp(KEY))
    # cliques = [
    #     (9, 15, 5), (9, 15, 14), (9, 12, 14),
    #     (12, 14, 2), (14, 2, 15), (5, 15, 2), (15, 14), (8, 9, 5)
    # ] # ensure that neighboring cliques have two atoms in common
    # cliques = [
    #     (5, 14, 2), (9, 15, 5), (9, 15, 14), (9, 12, 14),
    #     (15, 5, 14), (12, 14, 2), (15, 14), (8, 9, 5)
    # ]

    # cliques = [
    #     (9, 15, 5), (9, 15, 14), (9, 12, 14), (14, 2, 12),
    #     (15, 2, 14), (15, 5, 2), (15, 14)
    # ]
    # molTreeEdges = [
    #     (0, 6), (1, 6), (2, 6), (3, 6), (4, 6), (5, 6)
    # ]

    # complete
    cliques = [
        (9, 15, 5), (9, 15, 14), (9, 12, 14), (14, 2, 12),
        (15, 2, 14), (15, 5, 2), (15, 14), (1, 2, 12), 
        (2, 4, 5), (5, 9, 6), (9, 11, 12), (7, 9, 6)
    ]
    molTreeEdges = [
        (0, 6), (1, 6), (2, 6), (3, 6), (4, 6), (5, 6), 
        (3, 7), (5, 8), (0, 9), (2, 10), (9, 11)
    ]

    # original_graph = get_subgraph(mol, cliques)
    original_graph = mol_to_nx(mol)
    ghost_edges = [(9, 5), (9, 6), (9, 7), (9, 14), (9, 12), (9, 11), (12, 2), (12, 1), (12, 0), (2, 15), (2, 5), (2, 4)]
    triangulated_graph = add_ghost_edges(original_graph, ghost_edges)

    return mol, cliques, molTreeEdges, triangulated_graph

def gen_mol2():
    mol = Chem.MolFromSmiles("C1CCC2(CC1)CCCC1CCCCC21")
    for bond in mol.GetBonds(): bond.SetBoolProp(KEY, False)

    # complete
    cliques = [
        (3, 2, 1), (3, 1, 0), (3, 0, 5), (3, 4, 5),
        (3, 7, 6), (3, 7, 8), (3, 8, 9), (3, 9, 14),
        (9, 10, 14)
    ]
    molTreeEdges = [
        (0, 1), (1, 2), (2, 3), (1, 5), (4, 5), (5, 6), (6, 7), (7, 8) 
    ]

    # original_graph = get_subgraph(mol, cliques)
    original_graph = mol_to_nx(mol)
    ghost_edges = [(3, 0), (3, 1), (3, 5), (3, 7), (8, 3), (9, 3), (14, 10), (14, 11), (14, 12)]
    triangulated_graph = add_ghost_edges(original_graph, ghost_edges)

    return mol, cliques, molTreeEdges, triangulated_graph

def gen_mol3():
    mol = Chem.MolFromSmiles("C1CCC(CC1)C1CC2CCCC3CCCC(C1)C23")
    for bond in mol.GetBonds(): bond.SetBoolProp(KEY, False)

    # complete
    cliques = [
        (3, 6), (6, 7, 8), (6, 8, 17), (17, 8, 16), 
        (8, 16, 18), (18, 12, 8), (18, 12, 16), (18,), 
        (8, 12, 11), (8, 11, 10), (8, 10, 9), 
        (16, 12, 15), (15, 14, 12), (14, 12, 13)
    ]
    molTreeEdges = [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 7), (5, 7), (6, 7), (5, 8), (8, 9), (9, 10), (6, 11), (11, 12), (12, 13)
    ]

    # original_graph = get_subgraph(mol, cliques)
    original_graph = mol_to_nx(mol)
    ghost_edges = [(8, 6), (17, 8), (16, 8), (8, 12), (8, 11), (8, 10), (12, 16), (12, 15), (12, 14)]
    triangulated_graph = add_ghost_edges(original_graph, ghost_edges)

    return mol, cliques, molTreeEdges, triangulated_graph

def main():
    # ------------------------------------------------------------------------
    # mol, cliques, molTreeEdges, triangulated_graph = gen_mol1()
    # mol, cliques, molTreeEdges, triangulated_graph = gen_mol2()
    # mol, cliques, molTreeEdges, triangulated_graph = gen_mol3()

    with open("C:\\Users\\fongm\\Downloads\\icml18-jtnn\\data\\zinc\\all.txt") as f:
        smiles = f.readlines()
    # with open("/mnt/c/Users/fongm/Downloads/icml18-jtnn/data/zinc/all.txt") as f:
    #     smiles = f.readlines()

    # mol = Chem.MolFromSmiles("CCN1C(=O)C(=O)NC2=CC=CC=C12")
    # draw_mol(mol_to_nx(mol), 1, ['symbol', 'bond_type', 'color'])
    # mol_f = (nx_to_mol(mol_to_nx(mol)))
    # smiles = Chem.MolToSmiles(mol_f, kekuleSmiles=True)
    # print(Chem.MolToSmiles(mol, kekuleSmiles=True))

    # print(smiles)
    # raise


    # idx 3, 6, second honeycomb 
    # smiles = ["O=C1[C@@H]2C=C[C@@H](C=CC2)C1(c1ccccc1)c1ccccc1","O=C([O-])CC[C@@]12CCCC[C@]1(O)OC(=O)CC2", "ON=C1C[C@H]2CC3(C[C@@H](C1)c1ccccc12)OCCO3", "C[C@H]1CC(=O)[C@H]2[C@@]3(O)C(=O)c4cccc(O)c4[C@@H]4O[C@@]43[C@@H](O)C[C@]2(O)C1", 'Cc1cc(NC(=O)CSc2nnc3c4ccccc4n(C)c3n2)ccc1Br', 'CC(C)(C)c1ccc(C(=O)N[C@H]2CCN3CCCc4cccc2c43)cc1', "O=c1c2ccc3c(=O)n(-c4nccs4)c(=O)c4ccc(c(=O)n1-c1nccs1)c2c34", "O=C(N1CCc2c(F)ccc(F)c2C1)C1(O)Cc2ccccc2C1"]
    # chosen_smiles = smiles[5]

    chosen_smiles = smiles[56898]

    # chosen_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False)
    # chosen_smiles = "CCNC(=C)C12CCC(C)(C3=C1N=C1C=CC=CC1=N3)C2(C)C"
    # chosen_smiles = "CC1(C)C2(C)CCC1(C(N)=C)C1=C2N=C2C=CC=CC2=N1"
    print('ch sm', chosen_smiles)
    # chosen_smiles = smiles[379]
    mol = Chem.MolFromSmiles(chosen_smiles)
    cliques, molTreeEdges, triangulated_graph = tree_decomp(mol)

    # for i, cliq in enumerate(cliques): print(i + 1, cliq)
    # for c1, c2 in molTreeEdges:
    #     print(cliques[c1], cliques[c2])
    # raise

    list_of_nodes = []
    for i, clique in enumerate(cliques):
        if isinstance(clique, tuple):
            # print(i, clique, mol, get_fragments(mol ,list(clique)))
            c = list(clique)
            # for cw in list(clique):
            #     val = mol.GetAtomWithIdx(cw).GetSymbol()
            #     val2 = mol.GetAtomWithIdx(cw).GetNumExplicitHs()
            #     print(val, val2)

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

        # set_atommap(node.mol, node.nid)
        # set_atommap(node.tri_mol, node.nid)
        # set_atommap_graph(node.graph, node.nid)
        # node.is_leaf = (len(node.neighbors) == 1)

    for node in list_of_nodes:
        node.recover_G(triangulated_graph.copy())


    # for node in list_of_nodes:
    #     print(node.nid, [nei.nid for nei in node.neighbors])
    # raise

    #---------------------------------------------------------------------------
    
    # testing
    # root = list_of_nodes[-2]
    # root_idx = 9
    # root_idx = 11
    root_idx = 0
    # root_idx = 4
    # root_idx = 5
    root = list_of_nodes[root_idx]
    # root = list_of_nodes[6]
    cur_graph = root.graph.copy()
    global_amap = [{}] + [{} for node in list_of_nodes] # nid starts with 1, thus the first dict is not used, just as offset
    global_amap[root_idx+1] = {node:node for node in cur_graph.nodes()}

    dfs_assemble(cur_graph, global_amap, [], root, None)

    final_graph = cur_graph.copy()
    for a1, a2, data in cur_graph.edges(data=True):
        if data.get(KEY): final_graph.remove_edge(a1, a2)

    draw_mol(final_graph, 1, ['symbol', 'bond_type', 'color'])

    print()
    set_atommap_graph(final_graph)

    cur_mol = nx_to_mol(mol_to_nx(nx_to_mol(final_graph)))

    chosen_graph = mol_to_nx(Chem.MolFromSmiles(chosen_smiles))
    graph_match = nx.is_isomorphic(final_graph, chosen_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
    print("graph_match", graph_match)

    # -------------num_explicit_hss-----------
    # print('okisTT')
    # mol_t = nx_to_mol(mol_to_nx(mol))
    # # print('okis')
    # # for node in final_graph.nodes(data=True):
    # #     print(node)
    # mol_g = nx_to_mol(final_graph)
    # gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False, kekuleSmiles=True)
    # print(gold_smiles)
    # print(Chem.MolToSmiles(mol_t, isomericSmiles=False, kekuleSmiles=True))
    # print(Chem.MolToSmiles(mol_g, isomericSmiles=False, kekuleSmiles=True))
    # -----------------------------------------

    gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False)
    dec_smiles = Chem.MolToSmiles(cur_mol, isomericSmiles=False)
    
    print(chosen_smiles)
    print(gold_smiles)
    print(dec_smiles)
    print(gold_smiles == dec_smiles)

    # print(get_smiles(set_atommap(cur_mol)))
    # draw_mol(triangulated_graph, 700001, ['map_num', 'bond_type', 'color'])


main()


