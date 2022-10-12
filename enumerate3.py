
from rdkit import Chem
import itertools
import networkx as nx
from Cliquify.utils import *
import networkx.algorithms.isomorphism as iso
from Cliquify.debug_script import *
from rdkit import RDLogger
from Cliquify.tree_decomposition2 import tree_decomp
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
    try: return Chem.MolFragmentToSmiles(mol, atoms, kekuleSmiles=True)
    except: return Chem.MolFragmentToSmiles(mol, atoms)

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

        cand_amap = enum_attach(node.tri_mol, nei_node, cur_amap, singletons) # mol -> tri_mol
        if print_out: print("cand_amap", cand_amap, nei_node.nid)
        true_cand_graphs = []
        candidates = []
        for i, amap in enumerate(cand_amap):
            cand_mol = local_attach(node.tri_mol, neighbors[:depth+1], prev_nodes, amap) # mol -> tri_mol
            cand_graph = mol_to_nx(cand_mol)
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

def enum_assemble_singleton_tri2(cur_graph, cur_node, neighbors, prev_nodes, temp_global_amap, print_out=False): # neighbors exclude prev_nodes
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

    if print_out:
        draw_mol(cur_graph, 5999, ["map_num", "bond_type", "color"], folder="extension")
        draw_mol(starting_graph, 6000, ["map_num", "bond_type", "color"], folder="extension")
        # # draw_mol(cur_node.label_G, 6001, ["map_num", "bond_type", "color"], folder="extension")
        draw_mol(cur_node.graph, 6001, ["map_num", "bond_type", "color"], folder="extension")
        print('prev_nodes clique', prev_nodes[0].clique)
        print('cur_node clique', cur_node.clique)
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
            
            cand_graph = attach_graphs(cand_graph, [nei_node], prev_nodes, temp_global_amap, print_out=False) #father is already attached
            
            new_amap = cur_amap + [tuple(label_amap)]
            new_temp_global_amap = copy.deepcopy(temp_global_amap)
            new_cand_graph = cand_graph.copy()

            GM = iso.GraphMatcher(cur_node.label_G, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
            GM2 = iso.GraphMatcher(starting_graph, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
            GM3 = iso.GraphMatcher(cur_node.label_G, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
            if GM.subgraph_is_isomorphic() and not GM2.is_isomorphic() or GM3.is_isomorphic():
                if print_out: print((depth + 1) * 100 + i, True, label_amap)

                match_new_amaps.append(new_amap)
                match_new_temp_global_amaps.append(new_temp_global_amap)
                match_new_cand_graphs.append(new_cand_graph)
            else:
                mismatch_new_amaps.append(new_amap)
                mismatch_new_temp_global_amaps.append(new_temp_global_amap)
                mismatch_new_cand_graphs.append(new_cand_graph)

        if len(match_new_amaps) > 0:
            for i, new_amap in enumerate(match_new_amaps):
                search(match_new_cand_graphs[i], new_amap, match_new_temp_global_amaps[i], depth + 1)
        else:
            for i, new_amap in enumerate(mismatch_new_amaps):
                search(mismatch_new_cand_graphs[i], new_amap, mismatch_new_temp_global_amaps[i], depth + 1)

    search(starting_graph, [], temp_global_amap, 0)

    if print_out: print("num of cands inside:", len(all_attach_confs))

    return possible_cands_G, all_attach_confs

count = 0
#Only used for debugging purpose
def dfs_assemble(cur_graph, global_amap, fa_amap, cur_node, fa_node, print_out=False):
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
        cands, cands_G = None, [(nx_to_mol(cur_node.label_G), cur_node.label_G, [])]
    else:
        cands, cands_G = enum_assemble(cur_node, neighbors, prev_nodes, cur_amap)

    if (not cands) and (not cands_G):
        return

    label_idx = None
    cand_smiles, cand_Gs, cand_amap = zip(*cands_G)
    for i, cand_G in enumerate(cand_Gs):
        GM = iso.GraphMatcher(cur_node.label_G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)
        if GM.is_isomorphic():
            label_idx = i
            break

    # candidate printing
    if len(cand_Gs) > 1 and print_out: print('num of cands:', len(cand_Gs), "cur_id", cur_node.nid)
    
    # if not label_idx: return
    try:
        label_amap = cand_amap[label_idx]
    except:
        # draw_mol(cur_graph, 1, ["map_num", "bond_type", "color"], folder="next_cand", label='cur_graph')
        
        # draw_mol(cur_node.label_G, 999, ["map_num", "bond_type", "color"], folder="next_cand", label="label_G")
        for i, cand_G in enumerate(cand_Gs):
            # draw_mol(cand_G, 1000 + i, ["map_num", "bond_type", "color"], folder="next_cand")
            # GM = iso.GraphMatcher(cur_node.label_G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)
            GM = iso.GraphMatcher(cur_node.label_G, cand_G, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso) # exclude hss calculation
            # print(GM.is_isomorphic(), 1000 + i)
            if GM.is_isomorphic():
                label_idx = i
                break
        # raise
        # if not label_idx: return

    try:
        label_amap = cand_amap[label_idx]
    except:
        # print('cur_node.nid2', cur_node.nid)
        # draw_mol(cur_graph, 6000 + count, ["map_num", "bond_type", "color"], folder="cur_graph")

        # draw_mol(cur_node.label_G, 999, ["map_num", "bond_type", "color"], folder="next_cand", label="label_G")
        gold_smiles = Chem.MolToSmiles(nx_to_mol(cur_node.label_G), isomericSmiles=False)
        for i, smiles in enumerate(cand_smiles):
            if smiles == gold_smiles:
                label_idx = i
                # print('here1', i)
                break
            # print(smiles)
            # print(gold_smiles)
            # print()     
        # raise     
    
    try:
        label_amap = cand_amap[label_idx]
    except:
        for i, cand_G in enumerate(cand_Gs):
            GM = iso.GraphMatcher(cur_node.label_G, cand_G, node_match=node_equal_iso)
            if GM.is_isomorphic():
                # print('here2', i)
                label_idx = i
                break
        if not label_idx: return
        label_amap = cand_amap[label_idx]


    for nei_id,ctr_atom,nei_atom in label_amap:
        if nei_id == fa_nid:
            continue
        global_amap[nei_id][nei_atom] = global_amap[cur_node.nid][ctr_atom]

    if len(label_amap) == 0:

        nid_graph = {child.nid: child.graph.copy() for child in children}
        temp_global_amap = [amap for i, amap in enumerate(global_amap)] # filter all those amap which are not part of cur_node and prev_node
        cands_G_label, all_attach_confs = enum_assemble_singleton_tri2(cur_graph, cur_node, neighbors, prev_nodes, temp_global_amap, print_out=False)

        if print_out: print('num of cands:', len(all_attach_confs), "cur_id", cur_node.nid)
        chosen_i = None
        for i, total_label_amap in enumerate(all_attach_confs):
            GM = iso.GraphMatcher(cur_node.label_G, cands_G_label[i], node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
            if GM.is_isomorphic():
                # print(i, True)
                chosen_i  = i

        # if not chosen_i: return # honeycomb not working because of this

        total_label_amap = all_attach_confs[chosen_i]
        concat_label_amap = [label_amap for label_amap_list in total_label_amap for label_amap in label_amap_list]

        for nei_id,ctr_atom,nei_atom, ctr_id in concat_label_amap: # nei nid, atom idx, atom idx
            if nei_id == fa_nid: continue

            try: global_amap[nei_id][nei_atom] = global_amap[ctr_id][ctr_atom]
            except:
                # add artificial nodes
                global_amap[ctr_id][ctr_atom] = len(cur_graph.nodes)
                global_amap[nei_id][nei_atom] = global_amap[ctr_id][ctr_atom]

                node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
                cur_graph.add_node(len(cur_graph.nodes), **node_attr)
        
        cur_graph = attach_graphs(cur_graph, children, prev_nodes, global_amap)

        for i, nei_node in enumerate(children):
            if not nei_node.is_leaf:
                # return
                if i == len(children) - 1:
                    # print(total_label_amap[i])
                    label_amap = [amap[:3] for amap in total_label_amap[i]]
                else:
                    label_amap = [amap[:3] for amap in total_label_amap[i] if amap[3] == cur_node.nid]
                dfs_assemble(cur_graph, global_amap, label_amap, nei_node, cur_node)

    else:
        cur_graph = attach_graphs(cur_graph, children, [], global_amap) #father is already attached

        for nei_node in children:
            if not nei_node.is_leaf:
                dfs_assemble(cur_graph, global_amap, label_amap, nei_node, cur_node)

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


def main():

    with open("C:\\Users\\fongm\\Downloads\\icml18-jtnn\\data\\zinc\\all.txt") as f:
        smiles = f.readlines()


    for count, chosen_smiles in enumerate(smiles):
        # if count < 22400: continue
        # # if count % 200 == 0: print("over: ", count)
        chosen_smiles = smiles[425]

        # chosen_smiles = smiles[600]

        mol = Chem.MolFromSmiles(chosen_smiles)
        cliques, molTreeEdges, triangulated_graph = tree_decomp(mol)

        try:
            list_of_nodes = []
            for i, clique in enumerate(cliques):
                if isinstance(clique, tuple):
                    c = list(clique)
                    m = MolTreeNode(mol, c)
                list_of_nodes.append(m)
        except:
            with open("fail_mol_list.txt", "a") as myfile:
                gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False)
                myfile.writelines("{},{},Decompose\n".format(gold_smiles, count))
            continue

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

        dfs_assemble(cur_graph, global_amap, [], root, None, print_out=False)

        final_graph = cur_graph.copy()
        for a1, a2, data in cur_graph.edges(data=True):
            if data.get(KEY): final_graph.remove_edge(a1, a2)

        # draw_mol(final_graph, 1, ['symbol', 'bond_type', 'color'])

        # print()
        set_atommap_graph(final_graph)

        cur_mol = nx_to_mol(mol_to_nx(nx_to_mol(final_graph)))

        chosen_graph = mol_to_nx(Chem.MolFromSmiles(chosen_smiles))
        graph_match = nx.is_isomorphic(final_graph, chosen_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)


        gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(chosen_smiles), isomericSmiles=False)
        dec_smiles = Chem.MolToSmiles(cur_mol, isomericSmiles=False)

        print(gold_smiles)
        print(dec_smiles)
        break      

        if gold_smiles != dec_smiles or not graph_match:
            print("count", count)
            # print(chosen_smiles, end="")
            # print(gold_smiles)
            # print(dec_smiles)
            # print(gold_smiles == dec_smiles, graph_match)
            # print(cliques)

            # with open("debug_fail_mol.txt", "a") as myfile:
            #     myfile.writelines("count: {}\n".format(count))
            #     myfile.writelines(chosen_smiles)
            #     myfile.writelines(gold_smiles + "\n")
            #     myfile.writelines(dec_smiles + "\n")
            #     myfile.writelines("{} {}\n".format(gold_smiles == dec_smiles, graph_match))
            #     myfile.writelines(str(cliques) + "\n")
            #     myfile.writelines("\n")

            with open("fail_mol_list.txt", "a") as myfile:
                myfile.writelines("{},{},Reconstruct\n".format(gold_smiles, count))

            # print(get_smiles(set_atommap(cur_mol)))
            # draw_mol(triangulated_graph, 700001, ['map_num', 'bond_type', 'color'])

if __name__ == "__main__":
    main()


