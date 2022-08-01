import networkx as nx
from rdkit import Chem
# import matplotlib
import matplotlib.pyplot as plt


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
                   is_aromatic=atom.GetIsAromatic())
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
    node_to_idx = {}
    for node in G.nodes():
        a=Chem.Atom(atomic_nums[node])
        a.SetChiralTag(chiral_tags[node])
        a.SetFormalCharge(formal_charges[node])
        a.SetIsAromatic(node_is_aromatics[node])
        a.SetHybridization(node_hybridizations[node])
        a.SetNumExplicitHs(num_explicit_hss[node])
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