import rdkit
import rdkit.Chem as Chem
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from collections import defaultdict
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
# from Cliquify.utils import *
from utils import *
import networkx as nx
import copy

MST_MAX_WEIGHT = 100 

def add_ghost_edges(G, ghost_edges):

    for edge in ghost_edges:
        G.add_edge(edge[0],
                edge[1],
                bond_type=Chem.BondType.SINGLE,
                ghost=True,
                color='r',
                )
    return G

def tree_decomp(mol):
    n_atoms = mol.GetNumAtoms()
    for bond in mol.GetBonds(): bond.SetBoolProp(KEY, False)
    graph = mol_to_nx(mol)

    # ssr = sorted([sorted(list(x)) for x in Chem.GetSymmSSSR(mol)])
    # print(ssr)
    cliques = []
    non_ring_nodes = set()
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        if not bond.IsInRing():
            cliques.append([a1,a2])
            non_ring_nodes.add(a1)
            non_ring_nodes.add(a2)

    ssr_list = nx.cycle_basis(graph)

    cliques.extend(copy.deepcopy(ssr_list))

    nei_list = [[] for i in range(n_atoms)] # store your neighboring clique index
    nei_list2 = [[] for i in range(n_atoms)]
    for i in range(len(cliques)):
        for atom in cliques[i]:
            nei_list[atom].append(i)
            nei_list2[atom].append((cliques[i]))

    bridge_rings = set()
    #Merge Rings with intersection > 2 atoms
    for i in range(len(cliques)):
        if len(cliques[i]) <= 2: continue
        for atom in cliques[i]:
            for j in nei_list[atom]:
                if i >= j or len(cliques[j]) <= 2: continue
                inter = set(cliques[i]) & set(cliques[j])
                if len(inter) > 2:
                    bridge_rings.add(tuple(cliques[i]))
                    bridge_rings.add(tuple(cliques[j]))

                    cliques[i].extend(cliques[j])
                    cliques[i] = list(set(cliques[i]))
                    # combine two rings together and get their set of atoms, this would result in the bridge compound itself
                    cliques[j] = [] # will be cleared later
  

    ssr_dict = {}
    chosen_sele = set()
    for ring1 in ssr_list:
        if tuple(ring1) in bridge_rings: continue
        intersect_list = []
        for ring2 in ssr_list:
            if ring1 == ring2: continue
            # print(ring1, tuple(ring1) in bridge_rings)
            inter = set(ring1) & set(ring2)
            if len(inter) > 0 and inter not in intersect_list:
                intersect_list.append(inter)
        
        if not intersect_list: 
            ssr_dict[tuple(ring1)] = ring1[0]
        elif len(intersect_list) == 1:

            # [OPTIONAL] increase deterministic of ortho, meta, para attachment
            non_ring_nodes = set(ring1) & non_ring_nodes
            if non_ring_nodes:
                non_ring_nei = set(graph.neighbors(non_ring_nodes.pop())) & set(ring1)
                ssr_dict[tuple(ring1)] = non_ring_nei.pop()
                continue
            
            # spiro
            inter_atoms = list(intersect_list[0])
            if len(inter_atoms) == 1:
                neis = graph.neighbors(inter_atoms.pop())
                ssr_dict[tuple(ring1)] = list(set(neis) & set(ring1)).pop() # choose neighbors in ring
            else: # fused
                ssr_dict[tuple(ring1)] = inter_atoms[0] # choose intersect atoms
        else:
            # print(ring1, intersect_list)

            inter_comb = list(itertools.combinations(intersect_list, 2))
            inter_atoms = set() # intersection between ring intersection 
            union_inter = set()
            for comb in inter_comb:
                inter_atoms |= comb[0] & comb[1]
                union_inter.update(comb[0])
                union_inter.update(comb[1])

            if inter_atoms:
                # honeycomb structure
                avail_sele = (union_inter - chosen_sele) - inter_atoms
                if not avail_sele: 
                    ssr_dict[tuple(ring1)] = list(union_inter - inter_atoms).pop()
                    continue
                
                for inter in intersect_list:
                    if avail_sele & inter and inter_atoms & inter:
                        seleNode = inter - inter_atoms
                        ssr_dict[tuple(ring1)] = seleNode.pop()
                        break
            else:
                # non honeycomb
                ssr_dict[tuple(ring1)] = union_inter.pop()
        

    original_graph = graph.copy()
    for ring, seleNode in ssr_dict.items():
        neis = graph.neighbors(seleNode)
        ghost_nodes = set(ring) - set(neis) - {seleNode}
        ghost_edges = [(seleNode, node) for node in list(ghost_nodes)]
        add_ghost_edges(graph, ghost_edges) 
    
    original_cliques = copy.deepcopy(cliques)
    cliques = [cliq for cliq in cliques if cliq and tuple(cliq) not in list(ssr_dict.keys())]
    for ring in ssr_dict:
        temp_graph = graph.subgraph(ring)
        basis = nx.cycle_basis(temp_graph)
        cliques.extend(basis)

    # draw_mol(graph, 2, folder="tree_decomp_img")


    # ---------------------------------------------------------------------
    cliques = [c for c in cliques if len(c) > 0] # triangulated cliques
    nei_list = [[] for i in range(n_atoms)]
    nei_list2 = [[] for i in range(n_atoms)]
    for i in range(len(cliques)):
        for atom in cliques[i]:
            nei_list[atom].append(i)
            nei_list2[atom].append((cliques[i]))
    
    original_cliques = [c for c in original_cliques if len(c) > 0] # before triangulation cliques
    # print('ori cliques', original_cliques)
    original_nei_list = [[] for i in range(n_atoms)]
    for i in range(len(original_cliques)):
        for atom in original_cliques[i]:
            original_nei_list[atom].append(i)

    #Build edges and add singleton cliques
    edges = defaultdict(int)
    for atom in range(n_atoms):
        if len(nei_list[atom]) <= 1: 
            continue
        original_cnei = original_nei_list[atom]
        cnei = nei_list[atom]
        # cnei2 = nei_list2[atom]
        bonds = [c for c in original_cnei if len(original_cliques[c]) == 2]
        rings = [c for c in original_cnei if len(original_cliques[c]) >= 3]
        # rings = [c for c in cnei if len(cliques[c]) > 4]
        # print('rings', rings, 'cnei', cnei2, 'atom', atom)
        if len(bonds) > 2:
        # if len(bonds) > 2 or (len(bonds) == 2 and len(original_cnei) > 2): #In general, if len(cnei) >= 3, a singleton should be added, but 1 bond + 2 ring is currently not dealt with.
            cliques.append([atom])
            c2 = len(cliques) - 1
            for c1 in cnei:
                edges[(c1,c2)] = 1
        elif len(rings) > 2: #Multiple (n>2) complex rings
            # print('in here atom2', atom)
            if len(graph[atom]) > 3: continue # 3 triangles but not all connected

            second_node = None
            neis = original_graph.neighbors(atom)
            for nei in neis:
                temp_cnei = original_nei_list[nei]
                rings = [c for c in temp_cnei if len(cliques[c]) >= 3]
                if len(rings) > 2:
                    second_node = nei
                    break

            if second_node:
                if not sorted([atom, second_node]) in cliques:
                    second_cnei = nei_list[second_node]
                    cliques.append(sorted([atom, second_node]))
                    c2 = len(cliques) - 1
                    for c1 in cnei + second_cnei:
                        # print(cliques[c1], cliques[c2])
                        edges[(c1,c2)] = MST_MAX_WEIGHT - 1
            else:
                cliques.append([atom])
                c2 = len(cliques) - 1
                for c1 in cnei:
                    edges[(c1,c2)] = MST_MAX_WEIGHT - 1
        else:
            # print('in here atom3', atom)
            for i in range(len(cnei)):
                for j in range(i + 1, len(cnei)):
                    c1,c2 = cnei[i],cnei[j]
                    inter = set(cliques[c1]) & set(cliques[c2])
                    if edges[(c1,c2)] < len(inter):
                        edges[(c1,c2)] = len(inter) #cnei[i] < cnei[j] by construction
    # print()
    # # print(cliques)
    cliques = [tuple(cliq) for cliq in cliques]
    # for k, v in edges.items():
    #     # print(k, v)
    #     print(cliques[k[0]], cliques[k[1]], v)

    # edges = [u + (MST_MAX_WEIGHT-v,) for u,v in edges.iteritems()]
    edges = [u + (MST_MAX_WEIGHT-v,) for u,v in edges.items()]
    # print(edges)
    if len(edges) == 0:
        return cliques, edges

    #Compute Maximum Spanning Tree
    row,col,data = zip(*edges)
    n_clique = len(cliques)
    clique_graph = csr_matrix( (data,(row,col)), shape=(n_clique,n_clique) )
    junc_tree = minimum_spanning_tree(clique_graph)
    row,col = junc_tree.nonzero()
    edges = [(row[i],col[i]) for i in range(len(row))]
    return (cliques, edges, graph)    

if __name__ == "__main__":
    # mol = Chem.MolFromSmiles("ON=C1C[C@H]2CC3(C[C@@H](C1)c1ccccc12)OCCO3")
    # mol = Chem.MolFromSmiles("CC(C)(C)c1ccc(C(=O)N[C@H]2CCN3CCCc4cccc2c43)cc1")
    # mol = Chem.MolFromSmiles("O=c1c2ccc3c(=O)n(-c4nccs4)c(=O)c4ccc(c(=O)n1-c1nccs1)c2c34")
    # mol = Chem.MolFromSmiles("C1CCC2(CC1)CCCC1CCCCC21")
    # mol = Chem.MolFromSmiles("C12C3C1C23")
    # mol = Chem.MolFromSmiles("C1CC23CCC1(CC2)CC3")
    # mol = Chem.MolFromSmiles("C1CC2CCC1C1CCC2CC1")
    tree_decomp(mol)


