from rdkit import Chem
import networkx as nx
import matplotlib.pyplot as plt
# matplotlib.use("Agg")
import json




# Get the atoms of the symmetric smallest set of rings
def get_SSSR_from_mol(mol):
	# Call SSR method from rdkit
	ssr = Chem.GetSymmSSSR(mol)

	ssr_conv = []

	# For each SSR make it into a list
	for ring in ssr:
		l = list(ring)
		ssr_conv.append(l)

	# Return the SSR
	return ssr_conv

def conv_mol_to_nx(mol, SMILES, ssr):

    molGraph = nx.Graph()

    atomIR = set().union(*ssr)
    atomNNIR = []

    # add atoms into networkx graph
    for atom in mol.GetAtoms():
        index_atom = "{}_{}".format(atom.GetIdx(), atom.GetSymbol())

        if atom.GetIdx() not in atomIR:
            atomNNIR.append(index_atom)
        
        molGraph.add_nodes_from([index_atom])
        
    # bond that is not involved ring
    bondNNIR = []
    bondIR = []
    for ring in ssr:
        for i in range(len(ring)):
            bondIR.append({ring[i], ring[(i+1) % len(ring)]})

    # add bonds into networkx graph
    for bond in mol.GetBonds():

        startAtom = "{}_{}".format(
            bond.GetBeginAtomIdx(), 
            bond.GetBeginAtom().GetSymbol()
        )
        endAtom = "{}_{}".format(
            bond.GetEndAtomIdx(), 
            bond.GetEndAtom().GetSymbol()
        )
        bondType = bond.GetBondTypeAsDouble()

        if set([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]) not in bondIR:
            # bondNNIR.append([startAtom, endAtom, bondType])
            bondNNIR.append([startAtom, endAtom])

        # print(startAtom, endAtom, bondType)
        molGraph.add_edge(startAtom, endAtom, bondType=bondType)

    return molGraph, bondNNIR, atomNNIR


def rankSeleNode(molGraph, ring, intersection_atom_list):

    ring_sele = []
    for atom in ring:
        if atom in intersection_atom_list:
            ring_sele.append([atom, 0])	
        # atom neighbors not involved in intersection atom list
        elif len(set(molGraph.neighbors(atom)).intersection(set(intersection_atom_list)))==0:
            ring_sele.append([atom,3])

        # atom neighbors involved in intersection atom list
        elif len(set(molGraph.neighbors(atom)).intersection(set(intersection_atom_list)))==1:
            ring_sele.append([atom,2])

        # atom involved in intersection list
        else:
            ring_sele.append([atom,1])

    # seleNode for baseRing
    ring_sele.sort(key = lambda x: x[1])
    base_seleNode = ring_sele[-1][0]

    return base_seleNode

# add ghost edges to molecules to aid in decomposition
# 1. rank the rings connectivity
# 2. then add seleNode for each rings to start forming ghost edge
def ring_labeling(molGraph, SMILES, ssr):

    mol = Chem.MolFromSmiles(SMILES)

    labeled_ssr = []
    for ring in ssr:
        temp_ring = []
        for index in ring:
            atom_index = "{}_{}".format(
                index, mol.GetAtomWithIdx(index).GetSymbol()
            )
            temp_ring.append(atom_index)
        labeled_ssr.append(temp_ring)


    ring_pair_dictionary = dict()

    # keep track of bridge_ring pair
    bridge_ring_pair = []

    for ring in labeled_ssr:

        bond_list_dict = {
            "B": [],
            "F": [],
            "S": [],
            "N": []
        }

        for ring2 in labeled_ssr:
            if ring == ring2:
                continue

            # assumption: find intersection nodes between two rings
            # if they are 2 intersecting nodes, it cannot be two rings attanched to one ring
            #  due to the constraint of just comparing two rings
            intersection_atoms = set(ring).intersection(set(ring2))

            # bridge connection
            if len(intersection_atoms) >= 3:
                bond_list_dict["B"].append(intersection_atoms)
                bridge_ring_pair.append([set(ring), set(ring2)])

            # fused connection
            elif len(intersection_atoms) == 2:
                bond_list_dict["F"].append(intersection_atoms)

            # fused connection
            elif len(intersection_atoms) == 1:
                bond_list_dict["S"].append(intersection_atoms)

        atoms_involved_in_3_bonds = set()
        for atom_index in ring:
            # add starting atom which has 3 or more bonds
            val = [n for n in molGraph.neighbors(atom_index)]
            if len(val) > 2:
                atoms_involved_in_3_bonds.add(atom_index)

        # combining all sets
        temp_bond_list = sum(bond_list_dict.values(), [])
        current_intersect_atoms = set().union(*temp_bond_list)
        
        # minus off atoms associated with fused, bridged, spiro
        # left with atoms associated with bonds outside of ring 
        atoms_nnir = atoms_involved_in_3_bonds - current_intersect_atoms
        atoms_nnir = [[atom] for atom in atoms_nnir]

        if atoms_nnir:
            bond_list_dict["N"] = atoms_nnir
                
        ring_pair_dictionary[tuple(ring)] = bond_list_dict

    print("ring_pair_dictionary", ring_pair_dictionary)

    return ring_pair_dictionary, bridge_ring_pair

def formTriClique(molGraph, SMILES, ssr):

    ring_pair_dictionary, bridge_ring_pair = ring_labeling(molGraph, SMILES, ssr)

    # start adding ghost edge
    ghost_edge_list = []

    for ring, intersection_dict in ring_pair_dictionary.items():

        # seleNode, node chosen to start drawing ghost edges
        bridge_seleNode_list = []
        if intersection_dict.get("B"):
            
            # print("R", set(ring))

            # find seleNode for bridge:
            temp_ring_pair = []
            for ring_pair in bridge_ring_pair:
                temp_ring_pair = ring_pair if set(ring) in ring_pair else []
                break
            
            # another_ring = ring_pair[1] if set(ring).intersection(ring_pair[0]) == ring else ring_pair[0]

            intersection_atom_list = intersection_dict["B"][0]
            sum_of_bridge_atom = set().union(*temp_ring_pair)

            # choosing selenode for clique forming in the bridge
            for atom in intersection_atom_list:
                neighbors_list = [n for n in molGraph.neighbors(atom)]
                if len(sum_of_bridge_atom.intersection(set(neighbors_list))) == 3:
                    bridge_seleNode_list.append(atom)

            # print("bridge_seleNode_list" , bridge_seleNode_list)

            intersection_dict["seleNode"] = [bridge_seleNode_list[0]]
            
            # adding appropriate ghost edge
            for atom in intersection_atom_list:
                if set([bridge_seleNode_list[0], atom]) not in ghost_edge_list and atom != bridge_seleNode_list[0] and atom not in molGraph.neighbors(bridge_seleNode_list[0]):
                    ghost_edge_list.append(set([bridge_seleNode_list[0], atom]))

            non_bridge_atoms = set(ring) - set(intersection_atom_list)

            base_seleNode = rankSeleNode(molGraph, ring, intersection_atom_list)

            # print("BS", base_seleNode)

            # skip adding ghost edge to non bridge_selenode
            val = set.union(non_bridge_atoms, set(bridge_seleNode_list))
            val = val - set(base_seleNode)
            # print("val", val)

            intersection_dict["seleNode"].append(base_seleNode)

            for atom in val:
                if set([base_seleNode, atom]) not in ghost_edge_list and atom != base_seleNode and atom not in molGraph.neighbors(base_seleNode):
                    ghost_edge_list.append(set([base_seleNode, atom]))

        elif intersection_dict.get("F"):

            # find connection part of the ring
            conn_type_list = []
            for conn_type, intersection_atom_list in intersection_dict.items():
                # intersection_atom_list
                if intersection_atom_list:
                    conn_type_list.append(conn_type)

            # # if ring spiro and fused nodes connected
            if "S" in conn_type_list:
                comb_fused_spiro_atoms = []
                comb_fused_spiro_atoms.extend(intersection_dict.get("F"))
                comb_fused_spiro_atoms.extend(intersection_dict.get("S"))

                comb_fused_spiro_atoms = set().union(*comb_fused_spiro_atoms)
                seleNode = rankSeleNode(molGraph, ring, list(comb_fused_spiro_atoms))
                print("F_seleNode0", seleNode)

            else:
                # only fused is present
                intersection_atom_list = intersection_dict['F']
                # if the ring is connected by two or more rings (fused)
                if len(intersection_atom_list) > 1:
                    nodes_part_of_conn = set().union(*intersection_atom_list)

                    seleNode = rankSeleNode(molGraph, ring , list(nodes_part_of_conn))
                    print("F_seleNode1", seleNode)
                else:
                    intersection_atom_list = set().union(*intersection_atom_list)
                    seleNode = rankSeleNode(molGraph, ring, list(intersection_atom_list))

                    print("F_seleNode2", seleNode)

            intersection_dict["seleNode"] = [seleNode]

            for atom in ring:
                if set([seleNode, atom]) not in ghost_edge_list and atom != seleNode and atom not in molGraph.neighbors(seleNode):
                    ghost_edge_list.append(set([seleNode, atom]))

        elif intersection_dict.get("S"):
            intersection_atom_list = intersection_dict['S']
            if len(intersection_atom_list) > 1:
                nodes_part_of_conn = set().union(*intersection_atom_list)

                seleNode = rankSeleNode(molGraph, ring , list(nodes_part_of_conn))
                print("S_seleNode0", seleNode)
            else:
                intersection_atom_list = set().union(*intersection_atom_list)
                seleNode = list(intersection_atom_list)[0]
                print("S_seleNode1", seleNode)
            
            intersection_dict["seleNode"] = [seleNode]

            for atom in ring:
                if set([seleNode, atom]) not in ghost_edge_list and atom != seleNode and atom not in molGraph.neighbors(seleNode):
                    ghost_edge_list.append(set([seleNode, atom]))
        
        elif intersection_dict.get("N"):
            intersection_atom_list = intersection_dict['N']

            # if ring has more than one non ring connection
            if len(intersection_atom_list) > 1:
                nodes_part_of_conn = set().union(*intersection_atom_list)

                seleNode = rankSeleNode(molGraph, ring , list(nodes_part_of_conn))
                print("N_seleNode0", seleNode)
            else:
            # else, choose one atom beside the connection
                intersection_atom_list = set().union(*intersection_atom_list)
                intersectingNode = list(intersection_atom_list)[0]
                seleList = list(molGraph.neighbors(intersectingNode))
                for atom in seleList:
                    if atom in ring:
                        seleNode = atom
                        break
                
                print("N_seleNode1", seleNode)

            intersection_dict["seleNode"] = [seleNode]

            for atom in ring:
                if set([seleNode, atom]) not in ghost_edge_list and atom != seleNode and atom not in molGraph.neighbors(seleNode):
                    ghost_edge_list.append(set([seleNode, atom]))
            
    print("ghost_edge_list", ghost_edge_list)

    for edge in ghost_edge_list:
        molGraph.add_edge(*tuple(edge), bondType=0.0)

    return molGraph, ring_pair_dictionary


def checkOverlap(cliqList):
    cliqPairList = []

    for i in range(len(cliqList)):
        for j in range(i, len(cliqList)):
            cliq = cliqList[i]
            cliq2 = cliqList[j]

            if cliq == cliq2:
                continue
            elif len(set(cliq).intersection(set(cliq2))) == 2:
                cliqPairList.append([[cliq, cliq2] , 2])
            elif len(set(cliq).intersection(set(cliq2))) == 1:
                cliqPairList.append([[cliq, cliq2] , 1])
            else:
                cliqPairList.append([[cliq, cliq2] , 0])

    cliqPairList.sort(key=lambda x:x[1], reverse=True)
    # print("cliqPairList" ,cliqPairList)
    print()  

    return cliqPairList

def bond_encoding(treeNode, bond_dict):

    tempBondList = []
    for i in range(len(treeNode)):
        tempBondList.append([treeNode[i], treeNode[(i+1) % len(treeNode)]])

    # print("bond_dict", bond_dict)

    bondEncode = ""
    for bond in tempBondList:

        if bond_dict.get(tuple(bond)) is not None:
            bondEncode += str(bond_dict.get(tuple(bond)))
        else:
            bondEncode += str(bond_dict.get(tuple(bond[::-1])))
        
        bondEncode += "|"

    treeNodeEncode = "|".join(treeNode)
    # print("treeNodeEncode", treeNodeEncode)
    
    return [treeNodeEncode, bondEncode]

def conv_graph_to_tree(molGraph, ring_pair_dictionary, bondNNIR, atomNNIR):
    masterTree = nx.Graph()

    print()
    bond_dict = {}
    for line in nx.generate_edgelist(molGraph, data=['bondType']):
        startAtom, endAtom, bondType = line.split(" ")
        bond_dict[tuple([startAtom, endAtom])] = float(bondType)

    # print("bond_dict", bond_dict)

    edgeList = []

    for ring, intersection_dict in ring_pair_dictionary.items():
        print(ring)
        print(intersection_dict)

        seleList = intersection_dict.get("seleNode")
        cliqsList = []
        for seleNode in seleList:
            cliqList = nx.cliques_containing_node(molGraph, nodes=seleNode)
            cliqPairList = checkOverlap(cliqList)

            print("cliqList", cliqList)
            for cliqPair in cliqPairList:
                
                treeNodeStart = bond_encoding(cliqPair[0][0], bond_dict)
                treeNodeEnd = bond_encoding(cliqPair[0][1], bond_dict)

                edgeList.append([treeNodeStart, treeNodeEnd])

                if not cliqList:
                    break
                if cliqPair[0][0] in cliqList:
                    cliqList.remove(cliqPair[0][0])
                if cliqPair[0][1] in cliqList:
                    cliqList.remove(cliqPair[0][1])

            
        if intersection_dict.get("N"):
            nonRing_list = intersection_dict.get("N")
            for atom_nonRing in nonRing_list:
                print(atom_nonRing)
                cliqList_nonRing = nx.cliques_containing_node(molGraph, nodes=atom_nonRing[0])
                print("cliqList_N", cliqList_nonRing)

            raise
            print("N atom", intersection_dict.get("N"))
    
    print()
    print(bondNNIR)
    print(atomNNIR)
    print()

    print(edgeList)

    #         for cliqPair in cliqPairList:
                
    #         print("ring", ring)
    #         print("cliqs", cliqList)
    #         print("seleNode", seleNode)
    #         print()

    # for cliq in nx.enumerate_all_cliques(molGraph):
    #     print(cliq)




def draw_mol_graph(molGraph, figName):
	plt.clf()
	#cmap2 = matplotlib.colors.ListedColormap(matplotlib.cm.get_cmap("Set1").colors[:12])
	#cmap2 = plt.cm.Set1
	# nx.draw(molGraph, node_color=color_list, with_labels=True, cmap=color_list, font_weight='bold')
	nx.draw(molGraph, node_color="yellow", with_labels=True, font_weight='bold')    
	plt.savefig(figName)
	plt.clf()




def test():

    figname = "test_molecule"

    # SMILES = "CCC1CCC(=O)C2(C)CCCCC12" # basic molecule
    # SMILES = "C1CCC2CC(CCC2C1)C1CC2CCCCC2C2CCCCC12" # multiple ring test
    SMILES = "C1CCC2(CC1)CCC1C(CC(C3CCC4CCCCC4C3)C3CCCCC13)C2" # multiple ring test
    # SMILES = "CC1CCC2CCC3(CCCCC3)CC2C1" # multiple ring test
    SMILES = "C1CCC2CCCCC(C1)C2" # multiple ring test
    SMILES = "CC1C2CCCCC1CCCC2" # multiple ring test
    SMILES = "C1CCC2(CC1)CCCC1CCCCC21" # multiple ring test
    # SMILES = "CC1CCC2C(C1)C(C)CCC21CCCC(C)C1" # multiple ring test
    # SMILES = "C1CC2CCC1C2"
    # SMILES = "C1CCC2CC3CCCCC3CC2C1"
    # SMILES = "CC1CCCC2CC3C(C)CCC(C)C3CC12"
    # SMILES = "C1CCC2(CC1)CCCC1CCCCC21"
    # SMILES = "C1CCC2(CC1)CCC1CCCCC1C2"
    # SMILES = "C1CCC2(CC1)CCC1(CCCCC1)CC2"
    # SMILES = "C(C1CCCCC1)C1CCCC(C1)C1CCCCC1"

    # SMILES = "C1CCC(CC1)C1CCC2CCCCC(C2)C1" # basic bridge + NNIR
    SMILES = "C(C1CCCCC1)C1CCC2CCCCC(C1)C2" # basic bridge + NNIR


    # SMILES = "C1CCC2(CC1)CCCCC2" # basic spiro
    # SMILES = "C1CCC2CCCCC2C1" # basic fused
    # SMILES = "C1CCC2CCCCC(C1)C2" # basic bridged
    # SMILES = "C1CCC(CC1)C1CCCCC1" # basic one bond in between
    # SMILES = "C(C1CCCCC1)C1CCCCC1" # basic two bonds in between
    # SMILES = "C1CCC2CC3CCCCC3CC2C1" # basic ring in between


    mol = Chem.MolFromSmiles(SMILES)
    ssr = get_SSSR_from_mol(mol)
    # print(ssr)
    molGraph, bondNNIR, atomNNIR = conv_mol_to_nx(mol, SMILES, ssr)
    draw_mol_graph(molGraph=molGraph, figName=figname+".png")
    cliqGraph = molGraph.copy()
    cliqGraph, ring_pair_dictionary = formTriClique(cliqGraph, SMILES, ssr)
    tree = draw_mol_graph(molGraph=cliqGraph, figName=figname+"_out.png")

    tree = cliqGraph.copy()
    conv_graph_to_tree(tree, ring_pair_dictionary, bondNNIR, atomNNIR)
    
    


    # draw_mol_graph(molGraph=molGraph, figName=figname+"_out"+".png")





if __name__ == "__main__":
    test()