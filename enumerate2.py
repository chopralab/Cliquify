from rdkit import Chem
import itertools
import networkx as nx
from utils import *
import networkx.algorithms.isomorphism as iso
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

KEY = "gh"

def copy_atom(atom):
    new_atom = Chem.Atom(atom.GetSymbol())
    new_atom.SetFormalCharge(atom.GetFormalCharge())
    new_atom.SetAtomMapNum(atom.GetAtomMapNum())
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
    return Chem.MolFragmentToSmiles(mol, atoms, kekuleSmiles=True)



class MolTreeNode(object):

    def __init__(self, smiles, clique=[]):
        self.smiles = smiles
        self.mol = get_mol(self.smiles)
        self.tri_mol, self.graph = get_triangulated_graph(self.mol)

        self.clique = [x for x in clique] #copy
        self.neighbors = []
        
    def add_neighbor(self, nei_node):
        self.neighbors.append(nei_node)

    def recover_G(self, original_graph):
        clique = []
        clique.extend(self.clique)

        if not self.is_leaf:
            for cidx in self.clique:
                original_graph.nodes[cidx]["map_num"] = self.nid
        
        for nei_node in self.neighbors:
            clique.extend(nei_node.clique)
            if nei_node.is_leaf: #Leaf node, no need to mark 
                continue
            for cidx in nei_node.clique:
                #allow singleton node override the atom mapping
                if cidx not in self.clique or len(nei_node.clique) == 1: # neighboring node atom will only override non ctr clique atom
                    node = original_graph.nodes[cidx]
                    node["map_num"] = nei_node.nid
                    # print(cidx, nei_node.nid, 'current', self.nid)
        
        clique = list(set(clique))
        self.label_G = original_graph.subgraph(clique)
        # if self.nid == 7:
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
            if bond.GetBoolProp(KEY):
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

def attach_graphs(ctr_graph, neighbors, prev_nodes, nei_amap):
    prev_nids = [node.nid for node in prev_nodes]

    # print()
    # print(nei_amap)
    # print("nei_id : {nei_atom : ctr_atom }")

    for nei_node in prev_nodes + neighbors:
        nei_id,nei_graph = nei_node.nid, nei_node.graph # mol -> tri_mol
        amap = nei_amap[nei_id]
        for node in nei_graph.nodes():
            if node not in amap:
                node_attr = copy_node_attr(nei_graph, node)
                amap[node] = len(ctr_graph.nodes)
                # print(node_attr)
                ctr_graph.add_node(len(ctr_graph.nodes), **node_attr)

        if nei_graph.number_of_edges() == 0:
            nei_atom = nei_graph.nodes[0]
            ctr_atom = ctr_graph.nodes[amap[0]]
            ctr_graph.nodes[ctr_atom]["map_num"] = nei_graph.nodes[nei_atom]["map_num"]
        else:
            for node1, node2, data in nei_graph.edges(data=True):
                a1 = amap[node1]
                a2 = amap[node2]
                if not ctr_graph.has_edge(a1, a2):
                    ctr_graph.add_edge(a1, a2, **data)
                elif nei_id in prev_nids: #father node overrides
                    ctr_graph.remove_edge(a1, a2)
                    ctr_graph.add_edge(a1, a2, **data)

    return ctr_graph


def local_attach(ctr_mol, neighbors, prev_nodes, amap_list):
    ctr_mol = copy_edit_mol(ctr_mol)
    nei_amap = {nei.nid:{} for nei in prev_nodes + neighbors}

    for nei_id,ctr_atom,nei_atom in amap_list:
        nei_amap[nei_id][nei_atom] = ctr_atom

    ctr_mol = attach_mols(ctr_mol, neighbors, prev_nodes, nei_amap)
    return ctr_mol.GetMol()


MAX_NCAND = 200

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
   
    elif nei_mol.GetNumBonds() == 1: #neighbor is a bond
        bond = nei_mol.GetBondWithIdx(0)
        bond_val = int(bond.GetBondTypeAsDouble())
        b1,b2 = bond.GetBeginAtom(), bond.GetEndAtom()

        for atom in ctr_atoms: 
            #Optimize if atom is carbon (other atoms may change valence)
            if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() < bond_val:
                continue
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
                    if a1.GetAtomicNum() == 6 and a1.GetTotalNumHs() + a2.GetTotalNumHs() < 4:
                        continue
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
        # graph_data = []
        # for attr in [""]


        # edges_data = []
        # for n1, n2, data in self.graph.edges.data():
        #     data = json.dumps(data)
        #     edges_data.append((n1, n2, data))

        # nodes_data = []
        # for n, data in self.graph.nodes.data():
        #     data = json.dumps(data)
        #     nodes_data.append((n, data))
        # return hash((tuple(edges_data), tuple(nodes_data)))

def enum_assemble(node, neighbors, prev_nodes=[], prev_amap=[], print_out=False):
    all_attach_confs = []
    singletons = [nei_node.nid for nei_node in neighbors + prev_nodes if nei_node.graph.number_of_nodes() == 1]

    #------------------------ GET TO THIS AFTER ALL GRAPH CONV-----------------------------
    # tri_mol_count = [nei.tri_mol.GetNumAtoms() for nei in node.neighbors]
    # if node.graph.number_of_nodes() <= 2 and len(node.neighbors) >= 3 and sum(tri_mol_count) >= 3 * len(tri_mol_count):
    #     candidates = enum_assemble_singleton_tri(node)
    #     return candidates

    cands_id = []
    def search(cur_amap, depth):
        if len(all_attach_confs) > MAX_NCAND:
            return
        if depth == len(neighbors):
            all_attach_confs.append(cur_amap)
            return

        nei_node = neighbors[depth]

        cand_amap = enum_attach(node.tri_mol, nei_node, cur_amap, singletons) # mol -> tri_mol
        # print('cand_amap', cand_amap)
        cand_smiles = set()
        cand_graphs = set()
        true_cand_graphs = []
        candidates = []
        for i, amap in enumerate(cand_amap):
            cand_mol = local_attach(node.tri_mol, neighbors[:depth+1], prev_nodes, amap) # mol -> tri_mol

            if print_out:
                # draw_mol(mol_to_nx(cand_mol), 3000 + ((depth + 1) * 100) + i, ["map_num", "bond_type", "color"])
                GM = iso.GraphMatcher(node.label_G, mol_to_nx(cand_mol), node_match=node_equal_iso, edge_match=ring_edge_equal_iso)

                # duplicated_graph = len([1 for G in true_cand_graphs if nx.is_isomorphic(node.label_G, mol_to_nx(cand_mol), node_match=node_equal_iso)])
                # print(GM.is_isomorphic(), 
                #         30000 + ((depth + 1) * 1000) + i, 
                #         Chem.MolToSmiles(cand_mol), 
                #         "graph", Graph(mol_to_nx(cand_mol)) in cand_graphs,  
                #         "smiles", Chem.MolToSmiles(cand_mol) in cand_smiles,
                #         "t_graph", bool(duplicated_graph))
                print(GM.is_isomorphic(), Chem.MolToSmiles(cand_mol))

            smiles = Chem.MolToSmiles(cand_mol)
            # if smiles in cand_smiles and Graph(mol_to_nx(cand_mol)) in cand_graphs:
            #     continue
            if Graph(mol_to_nx(cand_mol)) in cand_graphs:
                continue
            # if smiles in cand_smiles:
            #     continue
            # duplicate = len([1 for G in true_cand_graphs if nx.is_isomorphic(node.label_G, mol_to_nx(cand_mol), node_match=node_equal_iso)])
            # if duplicate: continue
            

            cand_smiles.add(smiles)
            cand_graphs.add(Graph(mol_to_nx(cand_mol)))

            duplicate = len([1 for G in true_cand_graphs if nx.is_isomorphic(node.label_G, mol_to_nx(cand_mol), node_match=node_equal_iso)])
            if not duplicate: true_cand_graphs.append(mol_to_nx(cand_mol))
            cands_id.append(3000 + ((depth + 1) * 100) + i)

            candidates.append(amap)

        if len(candidates) == 0:
            return

        for new_amap in candidates:
            search(new_amap, depth + 1)

    search(prev_amap, 0)

    if print_out: print("cands_id", cands_id)

    cand_smiles = set()
    candidates = []
    candidates_G = []
    for i, amap in enumerate(all_attach_confs):
        cand_mol = local_attach(node.tri_mol, neighbors, prev_nodes, amap)
        # cand_mol = Chem.MolFromSmiles(Chem.MolToSmiles(cand_mol))
        cand_G = mol_to_nx(cand_mol)
        smiles = Chem.MolToSmiles(cand_mol)
        if smiles in cand_smiles:
            continue
        cand_smiles.add(smiles)
        Chem.Kekulize(cand_mol)
        candidates.append( (smiles,cand_mol,amap) )
        candidates_G.append( (smiles,cand_G,amap) )

    return candidates, candidates_G


            

#Only used for debugging purpose
def dfs_assemble(cur_graph, global_amap, fa_amap, cur_node, fa_node):
    fa_nid = fa_node.nid if fa_node is not None else -1
    prev_nodes = [fa_node] if fa_node is not None else []

    children = [nei for nei in cur_node.neighbors if nei.nid != fa_nid]

    neighbors = [nei for nei in children if nei.graph.number_of_nodes() > 1]
    neighbors = sorted(neighbors, key=lambda x:x.graph.number_of_nodes(), reverse=True)
    singletons = [nei for nei in children if nei.graph.number_of_nodes() == 1]
    neighbors = singletons + neighbors

    cur_amap = [(fa_nid,a2,a1) for nid,a1,a2 in fa_amap if nid == cur_node.nid] # check if there is any atommap has been occupied previously when attaching with parent
    cands, cands_G = enum_assemble(cur_node, neighbors, prev_nodes, cur_amap)

    if (not cands) and (not cands_G): return

    # cand_smiles, cand_mols, cand_amap = zip(*cands)
    # label_idx = cand_smiles.index(cur_node.label)
    # label_amap = cand_amap[label_idx]

    cand_smiles, cand_Gs, cand_amap = zip(*cands_G)
    for i, cand_G in enumerate(cand_Gs):
        GM = iso.GraphMatcher(cur_node.label_G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)
        if GM.is_isomorphic(): 
            label_idx = i
            break
    try:
        label_amap = cand_amap[label_idx]
    except:
        print("cur_node.nid", cur_node.nid)
        draw_mol(cur_node.label_G, 1999, ["map_num", "bond_type", "color"])
        enum_assemble(cur_node, neighbors, prev_nodes, cur_amap, print_out=True)

        for i, cand_G in enumerate(cand_Gs):
            draw_mol(cand_G, 2000 + i, ["map_num", "bond_type", "color"])
            print(cand_smiles[i])

    for nei_id,ctr_atom,nei_atom in label_amap:
        if nei_id == fa_nid:
            continue
        global_amap[nei_id][nei_atom] = global_amap[cur_node.nid][ctr_atom]
    
    # cur_mol = attach_mols(cur_mol, children, [], global_amap) #father is already attached
    cur_graph = attach_graphs(cur_graph, children, [], global_amap) #father is already attached

    for nei_node in children:
        if not nei_node.is_leaf:
            dfs_assemble(cur_graph, global_amap, label_amap, nei_node, cur_node)

def get_subgraph(mol, cliques):
    mol = Chem.Mol(mol)

    fragment_atoms = set([atom for cliq in cliques for atom in cliq])

    new_mol = Chem.RWMol(Chem.MolFromSmiles(''))
    for atom in mol.GetAtoms():
        new_atom = copy_atom(atom)
        new_mol.AddAtom(new_atom)

    selected_bonds = []
    for cliq in cliques:
        subset_possible_bonds = list(itertools.combinations(cliq, 2))
        for subset in subset_possible_bonds:
            if set(subset) not in selected_bonds: selected_bonds.append(set(subset))
    
    print(selected_bonds)
    for bond in selected_bonds:
        bond = list(bond)
        a1, a2 = bond[0], bond[1]
        bond_obj = mol.GetBondBetweenAtoms(a1, a2)
        if bond_obj:
            new_mol.AddBond(a1, a2, order=bond_obj.GetBondType())
            new_bond = new_mol.GetBondBetweenAtoms(a1, a2)
            new_bond.SetBoolProp(KEY, False)
        else:
            new_mol.AddBond(a1, a2, order=Chem.BondType.SINGLE)
            new_bond = new_mol.GetBondBetweenAtoms(a1, a2)
            new_bond.SetBoolProp(KEY, True)

    subgraph = mol_to_nx(new_mol, skip_unattached=True)

    return subgraph

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
        (2, 4, 5), (5, 9, 6), (9, 11, 12)
    ]
    molTreeEdges = [
        (0, 6), (1, 6), (2, 6), (3, 6), (4, 6), (5, 6), (3, 7), (5, 8), (0, 9), (2, 10)
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

def main():
    # ------------------------------------------------------------------------
    mol, cliques, molTreeEdges, triangulated_graph = gen_mol2()


    list_of_nodes = []
    for clique in cliques:
        if isinstance(clique, tuple):
            c = list(clique)
            m = MolTreeNode(get_fragments(mol ,list(clique)), c)
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


    #---------------------------------------------------------------------------
    
    # testing
    # root = list_of_nodes[-2]
    root = list_of_nodes[1]
    # root = list_of_nodes[6]
    cur_graph = root.graph.copy()
    global_amap = [{}] + [{} for node in list_of_nodes] # nid starts with 1, thus the first dict is not used, just as offset
    global_amap[2] = {node:node for node in cur_graph.nodes()}

    dfs_assemble(cur_graph, global_amap, [], root, None)

    draw_mol(cur_graph, 2000)
    draw_mol(triangulated_graph, 2001)

main()


