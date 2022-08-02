import networkx as nx
from rdkit import Chem
# import matplotlib
import matplotlib.pyplot as plt
import itertools
import networkx.algorithms.isomorphism as iso



KEY = "gh"
def mol_to_nx(mol, skip_unattached=False):
    G = nx.Graph()


    for atom in mol.GetAtoms():
        if skip_unattached and not atom.GetNeighbors(): continue
        G.add_node(atom.GetIdx(),
                   symbol=atom.GetSymbol(),
                   formal_charge=atom.GetFormalCharge(),
                   chiral_tag=atom.GetChiralTag(),
                   hybridization=atom.GetHybridization(),
                   num_explicit_hs=atom.GetNumExplicitHs(),
                   is_aromatic=atom.GetIsAromatic(),
                   map_num=atom.GetAtomMapNum())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType(),
                   ghost=bond.GetBoolProp(KEY),
                   color='r' if bond.GetBoolProp(KEY) else 'b',
                   )
    return G

def nx_to_mol(G):
    mol = Chem.RWMol()
    atomic_nums = nx.get_node_attributes(G, 'symbol')
    chiral_tags = nx.get_node_attributes(G, 'chiral_tag')
    formal_charges = nx.get_node_attributes(G, 'formal_charge')
    node_is_aromatics = nx.get_node_attributes(G, 'is_aromatic')
    node_hybridizations = nx.get_node_attributes(G, 'hybridization')
    num_explicit_hss = nx.get_node_attributes(G, 'num_explicit_hs')
    map_nums = nx.get_node_attributes(G, 'map_num')
    node_to_idx = {}
    for node in G.nodes():
        a=Chem.Atom(atomic_nums[node])
        a.SetChiralTag(chiral_tags[node])
        a.SetFormalCharge(formal_charges[node])
        a.SetIsAromatic(node_is_aromatics[node])
        a.SetHybridization(node_hybridizations[node])
        a.SetNumExplicitHs(num_explicit_hss[node])
        a.SetAtomMapNum(map_nums[node])
        idx = mol.AddAtom(a)
        node_to_idx[node] = idx

    bond_types = nx.get_edge_attributes(G, 'bond_type')
    ghost = nx.get_edge_attributes(G, 'ghost')
    for edge in G.edges():
        first, second = edge
        ifirst = node_to_idx[first]
        isecond = node_to_idx[second]
        bond_type = bond_types[first, second]
        mol.AddBond(ifirst, isecond, bond_type)
        new_bond = mol.GetBondBetweenAtoms(ifirst, isecond)
        new_bond.SetBoolProp(KEY, ghost[first, second])
        
    try: Chem.SanitizeMol(mol)
    except: pass
    return mol

# ON NETWORKX GRAPH

def node_equal_iso(node1, node2):
    return node1["symbol"] == node2["symbol"] and node1["formal_charge"] == node2["formal_charge"] \
        and node1["map_num"] == node2["map_num"]


def ring_edge_equal_iso(edge1, edge2):
    return edge1["bond_type"] == edge2["bond_type"] and \
        edge1["ghost"] == edge2["ghost"]

def ring_edge_equal_iso2(edge1, edge2):
    return edge1["bond_type"] == edge2["bond_type"]


def copy_node_attr(G, idx):
    val = {
        "symbol": G.nodes[idx]["symbol"],
        "chiral_tag": G.nodes[idx]["chiral_tag"],
        "formal_charge": G.nodes[idx]["formal_charge"],
        "is_aromatic": G.nodes[idx]["is_aromatic"],
        "hybridization": G.nodes[idx]["hybridization"],
        "num_explicit_hs": G.nodes[idx]["num_explicit_hs"],
    }

    return val

def node_equal(a1, a2):
    return a1["symbol"] == a2["symbol"] and a1["formal_charge"] == a2["formal_charge"]


def ring_edge_equal(G1, G2, b1, b2, reverse=False):
    bond_prop = G1.get_edge_data(*b1) == G2.get_edge_data(*b2)
    if reverse: b2 = b2[::-1]

    return node_equal(G1.nodes[b1[0]], G2.nodes[b2[0]]) and node_equal(G1.nodes[b1[1]], G2.nodes[b2[1]]) and bond_prop

def draw_mol(cand_G, numb=0, attr=['symbol', 'bond_type', 'color']):
    plt.clf()
    symbol, bond_type, color = attr

    pos = nx.spring_layout(cand_G)
    nx.draw(cand_G, pos)
    node_labels = nx.get_node_attributes(cand_G, symbol)
    node_labels = {k : "       ({})".format(v) for k, v in node_labels.items()}
    nx.draw_networkx_labels(cand_G, pos, node_labels)
    edge_labels = nx.get_edge_attributes(cand_G, bond_type)
    nx.draw_networkx_edge_labels(cand_G, pos, edge_labels)
    colors_cand_G = nx.get_edge_attributes(cand_G, color).values()
    nx.draw(cand_G, pos, node_color="yellow", with_labels=True, edge_color=colors_cand_G)
    # plt.show()
    plt.savefig("subgraph/show{}.png".format(numb))
    plt.clf()


    return


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