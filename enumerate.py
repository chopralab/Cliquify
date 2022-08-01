from re import sub
from matplotlib.pyplot import draw
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
        return mol, None

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

    def recover(self, original_graph):
        clique = []
        clique.extend(self.clique)

        if not self.is_leaf:
            for cidx in self.clique:
                original_graph.nodes[cidx]["MapNum"] = self.nid
            
            print(original_graph.nodes)
            raise
            




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
                triangulated_mol.AddBond(a1, a2, order=Chem.BondType.TRIPLE)
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

    print()
    print(nei_amap)
    print("nei_id : {nei_atom : ctr_atom }")

    for nei_node in prev_nodes + neighbors:
        nei_id,nei_mol = nei_node.nid,nei_node.mol
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
                print(a1, a2)
                if ctr_mol.GetBondBetweenAtoms(a1, a2) is None:
                    ctr_mol.AddBond(a1, a2, bond.GetBondType())
                elif nei_id in prev_nids: #father node overrides
                    ctr_mol.RemoveBond(a1, a2)
                    ctr_mol.AddBond(a1, a2, bond.GetBondType())
    return ctr_mol

def local_attach(ctr_mol, neighbors, prev_nodes, amap_list):
    ctr_mol = copy_edit_mol(ctr_mol)
    nei_amap = {nei.nid:{} for nei in prev_nodes + neighbors}

    for nei_id,ctr_atom,nei_atom in amap_list:
        nei_amap[nei_id][nei_atom] = ctr_atom

    ctr_mol = attach_mols(ctr_mol, neighbors, prev_nodes, nei_amap)
    return ctr_mol.GetMol()


MAX_NCAND = 200

#This version records idx mapping between ctr_mol and nei_mol
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


def local_attach_graph(cand_G, nei_G, amap):
    for node in nei_G.nodes():
        if node not in amap:
            node_attr = copy_node_attr(nei_G, node)
            amap[node] = len(cand_G.nodes)
            cand_G.add_node(len(cand_G.nodes), **node_attr)
    
    for node1, node2, data in nei_G.edges(data=True):
        a1 = amap[node1]
        a2 = amap[node2]
        if not cand_G.has_edge(a1, a2):
            cand_G.add_edge(a1, a2, **data)


def enum_attach_single_bond(ctr_G, nei_G):
    cands_G = []
    for b1 in ctr_G.edges():
        b1_st, b1_ed = b1[0], b1[1]
        for b2 in nei_G.edges():
            b2_st, b2_ed = b2[0], b2[1]

            if ring_edge_equal(ctr_G, nei_G, b1, b2):
                cand_G = ctr_G.copy()
                amap = {b2_st : b1_st, b2_ed: b1_ed}
                local_attach_graph(cand_G, nei_G, amap)

                # duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                # if not duplicate: cands_G.append(cand_G)
                cands_G.append(cand_G)

            elif ring_edge_equal(ctr_G, nei_G, b1, b2, reverse=True):
                cand_G = ctr_G.copy()
                amap = {b2_st : b1_ed, b2_ed: b1_st}
                local_attach_graph(cand_G, nei_G, amap)
                
                # duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                # if not duplicate: cands_G.append(cand_G)
                print("here")
                cands_G.append(cand_G)

    return cands_G
                
def enum_attach_double_bond(ctr_G, nei_G):
    ctr_pair_bonds = []
    for node in ctr_G.nodes():
        nei_n_ctr = [(nei, node) for nei in ctr_G.neighbors(node)]
        subset_possible_bonds = list(itertools.combinations(nei_n_ctr, 2))
        ctr_pair_bonds.extend(subset_possible_bonds)

    nei_pair_bonds = []
    for node in nei_G.nodes():
        nei_n_ctr = [(nei, node) for nei in nei_G.neighbors(node)]
        subset_possible_bonds = list(itertools.combinations(nei_n_ctr, 2))
        nei_pair_bonds.extend(subset_possible_bonds)

    cands_G = []
    for b1, b2 in ctr_pair_bonds:
        ctr_ctr_node = list(set(b1).intersection(set(b2)))[0]
        ctr_left_node = list(set(b1) - set(b1).intersection(set(b2)))[0]
        ctr_right_node = list(set(b2) - set(b1).intersection(set(b2)))[0]
        
        for b3, b4 in nei_pair_bonds:
            nei_ctr_node = list(set(b3).intersection(set(b4)))[0]
            nei_left_node = list(set(b3) - set(b3).intersection(set(b4)))[0]
            nei_right_node = list(set(b4) - set(b3).intersection(set(b4)))[0]
    
            if ring_edge_equal(ctr_G, nei_G, b1, b3) and ring_edge_equal(ctr_G, nei_G, b2, b4):
                cand_G = ctr_G.copy()
                amap = {nei_ctr_node : ctr_ctr_node, nei_left_node: ctr_left_node, nei_right_node: ctr_right_node}
                local_attach_graph(cand_G, nei_G, amap)

                duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                if not duplicate: cands_G.append(cand_G)
                # cands_G.append(cand_G)

            elif ring_edge_equal(ctr_G, nei_G, b1, b3, reverse=True) and ring_edge_equal(ctr_G, nei_G, b2, b4, reverse=True):
                cand_G = ctr_G.copy()
                amap = {nei_ctr_node : ctr_ctr_node, nei_left_node: ctr_right_node, nei_right_node: ctr_left_node}
                local_attach_graph(cand_G, nei_G, amap)

                duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                if not duplicate: cands_G.append(cand_G)
                # cands_G.append(cand_G)

    return cands_G
                
def enum_assemble_singleton_tri(node):

    neighbors = node.neighbors
    possible_cands_G = []
    cur_graph = neighbors[0].graph.copy()

    def search(cur_graph, depth):
        # print('current depth', depth)

        if depth == len(neighbors):
            possible_cands_G.append(cur_graph)
            return

        nei_node = neighbors[depth]
        cands_G = enum_attach_single_bond(cur_graph, nei_node.graph)

        if node.tri_mol.GetNumAtoms() == 1 and depth >= 2: # atom as singleton
            cands_G.extend(enum_attach_double_bond(cur_graph, nei_node.graph))

        if node.tri_mol.GetNumAtoms() == 2 and depth >= 3: # bond as singleton
            cands_G.extend(enum_attach_double_bond(cur_graph, nei_node.graph))

        # cands_G.extend(enum_attach_double_bond(cur_graph, nei_node.graph))

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
                
            filtered_graph.append(cand_G)
        
        if len(filtered_graph) == 0:
            return
        
        for cand_G in filtered_graph:
            search(cand_G, depth + 1)

    search(cur_graph, 1)

    print("candidates:", len(possible_cands_G))
    print("\n")

    draw_mol(node.label_G, 100000)
    for i, cand_G in enumerate(possible_cands_G):
        try: mol = nx_to_mol(cand_G)
        except: continue

        cand_smiles = Chem.MolToSmiles(mol)
        # graph_match = nx.is_isomorphic(node.label_G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)
        GM = iso.GraphMatcher(node.label_G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)
        graph_match = GM.is_isomorphic()
        smiles_match = node.label == cand_smiles
        if smiles_match and graph_match:
            print(cand_smiles)
            print("smiles_match", smiles_match)
            print("graph_match", graph_match)
            print()
            # draw_mol(cand_G, i)
    raise

    return possible_cands

def enum_assemble(node, neighbors, prev_nodes=[], prev_amap=[]):
    all_attach_confs = []
    singletons = [nei_node.nid for nei_node in neighbors + prev_nodes if nei_node.mol.GetNumAtoms() == 1]

    
    tri_mol_count = [nei.tri_mol.GetNumAtoms() for nei in node.neighbors]
    if node.mol.GetNumAtoms() <= 2 and len(node.neighbors) >= 3 and sum(tri_mol_count) >= 3 * len(tri_mol_count):
        candidates = enum_assemble_singleton_tri(node)
        # candidates = enum_assemble_singleton_tri(neighbors + prev_nodes)

        return candidates

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
        candidates = []
        for amap in cand_amap:
            cand_mol = local_attach(node.tri_mol, neighbors[:depth+1], prev_nodes, amap) # mol -> tri_mol
            # cand_test = remapping(cand_mol)
            # print(Chem.MolToSmiles(cand_test))
            # cand_mol = sanitize(cand_mol)
            # if cand_mol is None:
            #     continue
            # smiles = get_smiles(cand_mol)
            # if smiles in cand_smiles:
            #     continue
            # cand_smiles.add(smiles)
            # candidates.append(amap)

        # if len(candidates) == 0:
        #     return

        # for new_amap in candidates:
        #     search(new_amap, depth + 1)

    search(prev_amap, 0)
    # cand_smiles = set()
    # candidates = []
    # for amap in all_attach_confs:
    #     cand_mol = local_attach(node.mol, neighbors, prev_nodes, amap)
    #     cand_mol = Chem.MolFromSmiles(Chem.MolToSmiles(cand_mol))
    #     smiles = Chem.MolToSmiles(cand_mol)
    #     if smiles in cand_smiles:
    #         continue
    #     cand_smiles.add(smiles)
    #     Chem.Kekulize(cand_mol)
    #     candidates.append( (smiles,cand_mol,amap) )

    # return candidates

#Only used for debugging purpose
def dfs_assemble(cur_mol, global_amap, fa_amap, cur_node, fa_node):
    fa_nid = fa_node.nid if fa_node is not None else -1
    prev_nodes = [fa_node] if fa_node is not None else []

    children = [nei for nei in cur_node.neighbors if nei.nid != fa_nid]
    neighbors = [nei for nei in children if nei.mol.GetNumAtoms() > 1]
    neighbors = sorted(neighbors, key=lambda x:x.mol.GetNumAtoms(), reverse=True)
    singletons = [nei for nei in children if nei.mol.GetNumAtoms() == 1]
    neighbors = singletons + neighbors

    cur_amap = [(fa_nid,a2,a1) for nid,a1,a2 in fa_amap if nid == cur_node.nid] # check if there is any atommap has been occupied previously when attaching with parent
    cands = enum_assemble(cur_node, neighbors, prev_nodes, cur_amap)

    # if not cands: return

    # cand_smiles,cand_amap = zip(*cands)
    # label_idx = cand_smiles.index(cur_node.label)
    # label_amap = cand_amap[label_idx]

    # for nei_id,ctr_atom,nei_atom in label_amap:
    #     if nei_id == fa_nid:
    #         continue
    #     global_amap[nei_id][nei_atom] = global_amap[cur_node.nid][ctr_atom]
    
    # cur_mol = attach_mols(cur_mol, children, [], global_amap) #father is already attached
    # for nei_node in children:
    #     if not nei_node.is_leaf:
    #         dfs_assemble(cur_mol, global_amap, label_amap, nei_node, cur_node)

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

def main():
    # ------------------------------------------------------------------------

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
    cliques = [
        (9, 15, 5), (9, 15, 14), (9, 12, 14), (14, 2, 12),
        (15, 2, 14), (15, 5, 2), (15, 14)
    ]

    # original_graph = get_subgraph(mol, cliques)
    original_graph = mol_to_nx(mol)
    ghost_edges = [(9, 5), (9, 6), (9, 7), (9, 14), (9, 12), (9, 11), (12, 2), (12, 1), (12, 0), (2, 15), (2, 5), (2, 4)]
    triangulated_graph = add_ghost_edges(original_graph, ghost_edges)

    molTreeEdges = [
        (0, 6), (1, 6), (2, 6), (3, 6), (4, 6), (5, 6)
    ]

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
        node.is_leaf = (len(node.neighbors) == 1)

        # # label only for singleton 
        # tri_mol_count = [nei.tri_mol.GetNumAtoms() for nei in node.neighbors]
        # if node.mol.GetNumAtoms() <= 2 and len(node.neighbors) >= 3 and sum(tri_mol_count) >= 3 * len(tri_mol_count):
        #     cand_label = Chem.MolFromSmiles ("CC(C)C(C)C")
        #     for bond in cand_label.GetBonds():
        #         bond.SetBoolProp(KEY, False)

        #     cand_label = copy_edit_mol(cand_label)
        #     ghost_edges = [(1, 5), (0, 5), (4, 5), (2, 3), (2, 4), (0, 2)]
        #     # ghost_edges = [(1, 5), (0, 5), (4, 5), (2, 3), (2, 4)       ]

        #     for edge in ghost_edges:
        #         a1, a2 = edge
        #         cand_label.AddBond(a1, a2, order=Chem.BondType.SINGLE)
        #         new_bond = cand_label.GetBondBetweenAtoms(a1, a2)
        #         new_bond.SetBoolProp(KEY, True)

        #     cand_G_label = mol_to_nx(cand_label.GetMol())
        #     cand_smiles = Chem.MolToSmiles(cand_label.GetMol())
        #     node.label = cand_smiles
        #     node.label_G = cand_G_label
        #     print(cand_smiles)
            # draw_mol(node.label_G, 10000001)
            # raise

    for node in list_of_nodes:
        node.recover(triangulated_graph)

    raise


    #---------------------------------------------------------------------------
    
    # testing
    root = list_of_nodes[-2] # 
    cur_mol = copy_edit_mol(root.tri_mol)
    global_amap = [{}] + [{} for node in list_of_nodes]
    global_amap[1] = {atom.GetIdx():atom.GetIdx() for atom in cur_mol.GetAtoms()}

    dfs_assemble(cur_mol, global_amap, [], root, None)

main()
