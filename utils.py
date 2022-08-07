import networkx as nx
from rdkit import Chem
# import matplotlib
import matplotlib.pyplot as plt
import itertools
import networkx.algorithms.isomorphism as iso
import copy


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

def node_equal_iso2(node1, node2):
    return node1["symbol"] == node2["symbol"] and node1["formal_charge"] == node2["formal_charge"]

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
        "map_num": G.nodes[idx]["map_num"],
    }

    return val

def node_equal(a1, a2):
    return a1["symbol"] == a2["symbol"] and a1["formal_charge"] == a2["formal_charge"]


def ring_edge_equal(G1, G2, b1, b2, reverse=False):
    # bond_prop = G1.get_edge_data(*b1)["ghost"] == G2.get_edge_data(*b2)["ghost"]
    bond_prop = G1.get_edge_data(*b1) == G2.get_edge_data(*b2)
    if reverse: b2 = b2[::-1]

    return node_equal(G1.nodes[b1[0]], G2.nodes[b2[0]]) and node_equal(G1.nodes[b1[1]], G2.nodes[b2[1]]) and bond_prop

def draw_mol(cand_G, numb=0, attr=['symbol', 'bond_type', 'color'], folder="subgraph", label=""):
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
    plt.savefig("{}/show{}_{}.png".format(folder, numb, label))
    plt.clf()


    return


def attach_graphs(ctr_graph, neighbors, prev_nodes, nei_amap, print_out=False):
    prev_nids = [node.nid for node in prev_nodes]

    # print()
    # print(nei_amap)
    # print("nei_nid : {nei_atom : ctr_atom }")
    # print([nei_node.nid for nei_node in prev_nodes + neighbors])

    for nei_node in prev_nodes + neighbors:
        nei_id,nei_graph = nei_node.nid, nei_node.graph # mol -> tri_mol
        
        try: amap = nei_amap[nei_id]
        except: continue
        # amap = nei_amap[nei_id]
        for node in nei_graph.nodes():
            if node not in amap:
                node_attr = copy_node_attr(nei_graph, node)
                amap[node] = len(ctr_graph.nodes)
                # print(node_attr)
                ctr_graph.add_node(len(ctr_graph.nodes), **node_attr)

        if nei_graph.number_of_edges() == 0:
            nei_atom = nei_graph.nodes[0]
            ctr_atom = ctr_graph.nodes[amap[0]]
            ctr_atom["map_num"] = nei_atom["map_num"]
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

def label_amap_formation(ctr_G, nei_node, amap, global_amap, print_out=False):
    label_amap_G = []
    for nei_atom, ctr_atom in amap.items():
        nei_id = nei_node.nid
        ctr_id = ctr_G.nodes[ctr_atom]["map_num"] 
        if print_out: 
            print(amap)
            print((nei_id, ctr_atom, nei_atom, ctr_id))
            print(global_amap[ctr_id], ctr_id, ctr_atom)
            print(global_amap)
            print()

        ctr_atom = [fragment_idx for fragment_idx, full_graph_index in global_amap[ctr_id].items() if full_graph_index == ctr_atom ].pop()
        label_amap_G.append((nei_id, ctr_atom, nei_atom, ctr_id))
    return label_amap_G

count = 0
def enum_attach_single_bond(ctr_G, nei_node, global_amap):
    nei_G = nei_node.graph
    cands_G = []
    cands_G_amap = []
    for b1 in ctr_G.edges():
        b1_st, b1_ed = b1[0], b1[1]
        for b2 in nei_G.edges():
            b2_st, b2_ed = b2[0], b2[1]

            if ring_edge_equal(ctr_G, nei_G, b1, b2):
                cand_G = ctr_G.copy()
                cand_global_amap = copy.deepcopy(global_amap)

                amap = {b2_st : b1_st, b2_ed: b1_ed}
                try:
                    label_amap_G = label_amap_formation(ctr_G, nei_node, amap, cand_global_amap, print_out=False)
                    local_attach_graph(cand_G, nei_G, amap)

                    duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                    if not duplicate: 
                        cands_G.append(cand_G)
                        cands_G_amap.append(label_amap_G) #, print(label_amap_G)
                    # cands_G.append(cand_G)
                except: pass

            if ring_edge_equal(ctr_G, nei_G, b1, b2, reverse=True):
                cand_G = ctr_G.copy()
                cand_global_amap = copy.deepcopy(global_amap)

                amap = {b2_st : b1_ed, b2_ed: b1_st}
                try:
                    label_amap_G = label_amap_formation(ctr_G, nei_node, amap, cand_global_amap, print_out=False)
                    local_attach_graph(cand_G, nei_G, amap)
                    
                    duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                    if not duplicate: 
                        cands_G.append(cand_G)
                        cands_G_amap.append(label_amap_G) #, print(label_amap_G)
                # cands_G.append(cand_G)
                except: pass

    return cands_G, cands_G_amap
                
def enum_attach_double_bond(ctr_G, nei_node, global_amap):
    nei_G = nei_node.graph

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
    cands_G_amap = []
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
                cand_global_amap = copy.deepcopy(global_amap)

                amap = {nei_ctr_node : ctr_ctr_node, nei_left_node: ctr_left_node, nei_right_node: ctr_right_node}
                try:
                    label_amap_G = label_amap_formation(ctr_G, nei_node, amap, cand_global_amap)
                    local_attach_graph(cand_G, nei_G, amap)

                    duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                    if not duplicate: 
                        cands_G.append(cand_G)
                        cands_G_amap.append(label_amap_G) #, print(label_amap_G)
                except:
                    pass
                # cands_G.append(cand_G)

            if ring_edge_equal(ctr_G, nei_G, b1, b3, reverse=True) and ring_edge_equal(ctr_G, nei_G, b2, b4, reverse=True):
                cand_G = ctr_G.copy()
                cand_global_amap = copy.deepcopy(global_amap)

                amap = {nei_ctr_node : ctr_ctr_node, nei_left_node: ctr_right_node, nei_right_node: ctr_left_node}
                try:
                    label_amap_G = label_amap_formation(ctr_G, nei_node, amap, cand_global_amap)
                    local_attach_graph(cand_G, nei_G, amap)

                    duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                    if not duplicate: 
                        cands_G.append(cand_G)
                        cands_G_amap.append(label_amap_G) #, print(label_amap_G)
                except:
                    pass
                # cands_G.append(cand_G)

    return cands_G, cands_G_amap

def enum_attach_double_bond2(ctr_G, nei_node, global_amap):
    nei_G = nei_node.graph

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
    cands_G_amap = []
    for b1, b2 in ctr_pair_bonds:
        ctr_ctr_node = list(set(b1).intersection(set(b2)))[0]
        ctr_left_node = list(set(b1) - set(b1).intersection(set(b2)))[0]
        ctr_right_node = list(set(b2) - set(b1).intersection(set(b2)))[0]
        
        for b3, b4 in nei_pair_bonds:
            nei_ctr_node = list(set(b3).intersection(set(b4)))[0]
            nei_left_node = list(set(b3) - set(b3).intersection(set(b4)))[0]
            nei_right_node = list(set(b4) - set(b3).intersection(set(b4)))[0]
                
            if ring_edge_equal(ctr_G, nei_G, (ctr_left_node, ctr_ctr_node), (nei_left_node, nei_ctr_node)) and ring_edge_equal(ctr_G, nei_G, (ctr_ctr_node, ctr_right_node), (nei_ctr_node, nei_right_node)):
                cand_G = ctr_G.copy()
                cand_global_amap = copy.deepcopy(global_amap)
                amap = {nei_ctr_node : ctr_ctr_node, nei_left_node: ctr_left_node, nei_right_node: ctr_right_node}
                # try:
                #     label_amap_G = label_amap_formation(ctr_G, nei_node, amap, cand_global_amap)
                #     local_attach_graph(cand_G, nei_G, amap)

                #     duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                #     if not duplicate: 
                #         cands_G.append(cand_G)
                #         cands_G_amap.append(label_amap_G) #, print(label_amap_G)
                # except: pass

                label_amap_G = label_amap_formation(ctr_G, nei_node, amap, cand_global_amap, print_out=False)
                local_attach_graph(cand_G, nei_G, amap)

                duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                if not duplicate: 
                    cands_G.append(cand_G)
                    cands_G_amap.append(label_amap_G) #, print(label_amap_G)

            if ring_edge_equal(ctr_G, nei_G, (ctr_right_node, ctr_ctr_node), (nei_left_node, nei_ctr_node)) and ring_edge_equal(ctr_G, nei_G, (ctr_ctr_node, ctr_left_node), (nei_ctr_node, nei_right_node)):
            # if ring_edge_equal(ctr_G, nei_G, b1, b3, reverse=True) and ring_edge_equal(ctr_G, nei_G, b2, b4):
                cand_G = ctr_G.copy()
                cand_global_amap = copy.deepcopy(global_amap)

                amap = {nei_ctr_node : ctr_ctr_node, nei_left_node: ctr_right_node, nei_right_node: ctr_left_node}
                # try:
                #     label_amap_G = label_amap_formation(ctr_G, nei_node, amap, cand_global_amap)
                #     local_attach_graph(cand_G, nei_G, amap)

                #     duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                #     if not duplicate: 
                #         cands_G.append(cand_G)
                #         cands_G_amap.append(label_amap_G) #, print(label_amap_G)
                # except: pass
                label_amap_G = label_amap_formation(ctr_G, nei_node, amap, cand_global_amap,  print_out=False)
                local_attach_graph(cand_G, nei_G, amap)

                duplicate = len([1 for G in cands_G if nx.is_isomorphic(G, cand_G, node_match=node_equal_iso, edge_match=ring_edge_equal_iso)])
                if not duplicate: 
                    cands_G.append(cand_G)
                    cands_G_amap.append(label_amap_G) #, print(label_amap_G)
        # print("-----------------------------------------------------")

    return cands_G, cands_G_amap

