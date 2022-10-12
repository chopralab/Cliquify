      
from Cliquify.utils import *

def debug_singleton_3nei(cur_graph, cur_node, children, prev_nodes, fa_nid, global_amap):

    print([(child.clique, child.nid) for child in children])
    print()
    # temp_global_amap = [amap if i == cur_node.nid or i == fa_nid else {} for i, amap in enumerate(global_amap) ]
    temp_global_amap = [amap for i, amap in enumerate(global_amap)]
    nid_graph = {child.nid: child.graph.copy() for child in children}

    cur_clique = list(global_amap[cur_node.nid].values())
    fa_clique = list(global_amap[fa_nid].values())
    label_starting_clique = set(cur_clique + fa_clique)
    non_label_node = list(set(cur_graph.nodes()) - label_starting_clique)

    starting_graph = cur_graph.subgraph(label_starting_clique).copy()
    
    total_label_amap = []

    print('a', temp_global_amap)
    
    # --------------------------------------------------------------------------------
    print("\First iteration")
    # cur_graph = validate_graph.copy()
    # print('a', temp_global_amap)
    cands_G, cands_G_amap = enum_attach_single_bond(starting_graph, children[0], temp_global_amap)
    candidate_graphs = [starting_graph.copy() for _ in cands_G_amap]
    temp_global_amaps = [copy.deepcopy(temp_global_amap) for _ in cands_G_amap]
    
    chosen_amap = None
    chosen_graph = None
    draw_mol(starting_graph, 990, ["map_num", "bond_type", "color"], folder="extension", label="starting_graph")
    draw_mol(children[0].graph, 991, ["map_num", "bond_type", "color"], folder="extension", label="children[0]")
    draw_mol(cur_node.label_G, 992, ["map_num", "bond_type", "color"], folder="extension", label="label_G")

    for i, label_amap in enumerate(cands_G_amap):
        cand_graph = candidate_graphs[i]
        temp_global_amap = temp_global_amaps[i]

        for nei_id, ctr_atom,nei_atom, ctr_id in label_amap: # nei nid, atom idx, atom idx
            if nei_id == fa_nid: continue

            try: temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]
            except:
                # add artificial nodes
                temp_global_amap[ctr_id][ctr_atom] = len(cand_graph.nodes)
                temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]

                node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
                cand_graph.add_node(len(cand_graph.nodes), **node_attr)
        
        cand_graph = attach_graphs(cand_graph, [children[0]], prev_nodes, temp_global_amap, print_out=True) #father is already attached
        
        # validate_graph = cand_graph.copy()
        # validate_graph.remove_nodes_from(non_label_node)

        GM = iso.GraphMatcher(cur_node.label_G, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso2)
        GM2 = iso.GraphMatcher(cur_node.label_G, cands_G[i], node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        if GM.subgraph_is_isomorphic(): 
        # if (GM.subgraph_is_isomorphic() and GM2.subgraph_is_isomorphic()) or GM2.subgraph_is_isomorphic(): 
            # draw_mol(cands_G[i], 1000 + i, ["map_num", "bond_type", "color"], folder="extension")
            # draw_mol(cand_graph, 1000 + i, ["map_num", "bond_type", "color"], folder="extension")
            print(1000 + i, True, label_amap)
            print(temp_global_amap)
        else:
            print(1000 + i, False)
            draw_mol(cand_graph, 1000 + i, ["map_num", "bond_type", "color"], folder="extension")
            if i == 0:
                print(1000 + i, True)
                total_label_amap.append(label_amap)
                chosen_amap = copy.deepcopy(temp_global_amap)
                chosen_graph = cand_graph.copy()

    # --------------------------------------------------------------------------------
    temp_global_amap = chosen_amap
    starting_graph = chosen_graph
    print("\nSecond iteration")
    cands_G, cands_G_amap = enum_attach_double_bond2(starting_graph, children[1], temp_global_amap)
    candidate_graphs = [starting_graph.copy() for _ in cands_G_amap]
    temp_global_amaps = [copy.deepcopy(temp_global_amap) for _ in cands_G_amap]
    draw_mol(starting_graph, 1990, ["map_num", "bond_type", "color"], folder="extension", label="starting_graph")
    draw_mol(children[1].graph, 1991, ["map_num", "bond_type", "color"], folder="extension", label="children[1]")
    draw_mol(cur_node.label_G, 1992, ["map_num", "bond_type", "color"], folder="extension", label="label_G")
    
    chosen_amap = None
    chosen_graph = None
    for i, label_amap in enumerate(cands_G_amap):
        cand_graph = candidate_graphs[i]
        temp_global_amap = temp_global_amaps[i]

        for nei_id, ctr_atom,nei_atom, ctr_id in label_amap: # nei nid, atom idx, atom idx
            if nei_id == fa_nid: continue

            try: temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]
            except:
                # add artificial nodes
                temp_global_amap[ctr_id][ctr_atom] = len(cand_graph.nodes)
                temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]

                node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
                cand_graph.add_node(len(cand_graph.nodes), **node_attr)
        
        cand_graph = attach_graphs(cand_graph, [children[1]], prev_nodes, temp_global_amap, print_out=True) #father is already attached
        
        # validate_graph = cand_graph.copy()
        # validate_graph.remove_nodes_from(non_label_node)

        GM = iso.GraphMatcher(cur_node.label_G, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        GM2 = iso.GraphMatcher(cur_node.label_G, cands_G[i], node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        if GM.subgraph_is_isomorphic(): 
        # if (GM.subgraph_is_isomorphic() and GM2.subgraph_is_isomorphic()) or GM2.subgraph_is_isomorphic(): 
            # draw_mol(cands_G[i], 2000 + i, ["map_num", "bond_type", "color"], folder="extension")
            print(2000 + i, True, label_amap)
            draw_mol(cand_graph, 2000 + i, ["map_num", "bond_type", "color"], folder="extension")
            print(temp_global_amap)
        else:
            print(2000 + i, False, label_amap)
            draw_mol(cand_graph, 2000 + i, ["map_num", "bond_type", "color"], folder="extension")
            # if i == 3: 
            #     total_label_amap.append(label_amap)
            #     chosen_amap = copy.deepcopy(temp_global_amap)
            #     chosen_graph = cand_graph.copy()

    return

def debug_singleton_6nei(cur_graph, cur_node, children, prev_nodes, fa_nid, global_amap):

    print([(child.clique, child.nid) for child in children])
    print()
    # temp_global_amap = [amap if i == cur_node.nid or i == fa_nid else {} for i, amap in enumerate(global_amap) ]
    temp_global_amap = [amap for i, amap in enumerate(global_amap)]
    nid_graph = {child.nid: child.graph.copy() for child in children}

    cur_clique = list(global_amap[cur_node.nid].values())
    fa_clique = list(global_amap[fa_nid].values())
    label_starting_clique = set(cur_clique + fa_clique)
    non_label_node = list(set(cur_graph.nodes()) - label_starting_clique)

    starting_graph = cur_graph.subgraph(label_starting_clique).copy()
    
    total_label_amap = []

    print('a', temp_global_amap)
    
    # -----------------------------------------------------------------------------------
    # cands_G, cands_G_amap = enum_attach_double_bond2(starting_graph, children[0], temp_global_amap)

    # label_amap = cands_G_amap[1]
    # # print('bef', temp_global_amap)
    # for nei_id,ctr_atom,nei_atom, ctr_id in label_amap: # nei nid, atom idx, atom idx
    #     if nei_id == fa_nid:
    #         continue
    #     try:
    #         temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]
    #     except:
    #         # add artificial nodes
    #         temp_global_amap[ctr_id][ctr_atom] = len(starting_graph.nodes)
    #         temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]

    #         node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
    #         starting_graph.add_node(len(starting_graph.nodes), **node_attr)
    
    # # print('aft', temp_global_amap)
    # starting_graph = attach_graphs(starting_graph, [children[0]], prev_nodes, temp_global_amap) #father is already attached

    # # draw_mol(starting_graph, 1003, ["map_num", "bond_type", "color"], folder="extension")
    # validate_graph = starting_graph.copy()
    # # validate_graph.remove_nodes_from(non_label_node)
    # GM = iso.GraphMatcher(cur_node.label_G, validate_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
    # if GM.subgraph_is_isomorphic(): print(True, i)

    # total_label_amap.append(label_amap)
    # print('b', temp_global_amap)

    # --------------------------------------------------------------------------------
    print("\First iteration")
    # cur_graph = validate_graph.copy()
    # print('a', temp_global_amap)
    cands_G, cands_G_amap = enum_attach_double_bond2(starting_graph, children[0], temp_global_amap)
    candidate_graphs = [starting_graph.copy() for _ in cands_G_amap]
    temp_global_amaps = [copy.deepcopy(temp_global_amap) for _ in cands_G_amap]
    
    chosen_amap = None
    chosen_graph = None
    draw_mol(starting_graph, 990, ["map_num", "bond_type", "color"], folder="extension", label="starting_graph")
    draw_mol(children[0].graph, 991, ["map_num", "bond_type", "color"], folder="extension", label="children[0]")
    draw_mol(cur_node.label_G, 992, ["map_num", "bond_type", "color"], folder="extension", label="label_G")

    for i, label_amap in enumerate(cands_G_amap):
        cand_graph = candidate_graphs[i]
        temp_global_amap = temp_global_amaps[i]

        for nei_id, ctr_atom,nei_atom, ctr_id in label_amap: # nei nid, atom idx, atom idx
            if nei_id == fa_nid: continue

            try: temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]
            except:
                # add artificial nodes
                temp_global_amap[ctr_id][ctr_atom] = len(cand_graph.nodes)
                temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]

                node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
                cand_graph.add_node(len(cand_graph.nodes), **node_attr)
        
        cand_graph = attach_graphs(cand_graph, [children[0]], prev_nodes, temp_global_amap, print_out=True) #father is already attached
        
        # validate_graph = cand_graph.copy()
        # validate_graph.remove_nodes_from(non_label_node)

        GM = iso.GraphMatcher(cur_node.label_G, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        GM2 = iso.GraphMatcher(cur_node.label_G, cands_G[i], node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        if GM.subgraph_is_isomorphic(): 
        # if (GM.subgraph_is_isomorphic() and GM2.subgraph_is_isomorphic()) or GM2.subgraph_is_isomorphic(): 
            # draw_mol(cands_G[i], 1000 + i, ["map_num", "bond_type", "color"], folder="extension")
            # draw_mol(cand_graph, 1000 + i, ["map_num", "bond_type", "color"], folder="extension")
            # print(1000 + i, True, label_amap)
            print(temp_global_amap)
            if i == 3: 
                total_label_amap.append(label_amap)
                chosen_amap = copy.deepcopy(temp_global_amap)
                chosen_graph = cand_graph.copy()
        else:
            print(1000 + i, False)
            draw_mol(cand_graph, 1000 + i, ["map_num", "bond_type", "color"], folder="extension")

    
    # --------------------------------------------------------------------------------
    print("\nSecond iteration")
    temp_global_amap = chosen_amap
    starting_graph = chosen_graph
    cands_G, cands_G_amap = enum_attach_single_bond(starting_graph, children[1], temp_global_amap)
    candidate_graphs = [starting_graph.copy() for _ in cands_G_amap]
    temp_global_amaps = [copy.deepcopy(temp_global_amap) for _ in cands_G_amap]
    
    chosen_amap = None
    chosen_graph = None
    for i, label_amap in enumerate(cands_G_amap):
        cand_graph = candidate_graphs[i]
        temp_global_amap = temp_global_amaps[i]

        for nei_id, ctr_atom,nei_atom, ctr_id in label_amap: # nei nid, atom idx, atom idx
            if nei_id == fa_nid: continue

            try: temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]
            except:
                # add artificial nodes
                temp_global_amap[ctr_id][ctr_atom] = len(cand_graph.nodes)
                temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]

                node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
                cand_graph.add_node(len(cand_graph.nodes), **node_attr)
        
        cand_graph = attach_graphs(cand_graph, [children[1]], prev_nodes, temp_global_amap, print_out=True) #father is already attached
        
        # validate_graph = cand_graph.copy()
        # validate_graph.remove_nodes_from(non_label_node)

        GM = iso.GraphMatcher(cur_node.label_G, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        GM2 = iso.GraphMatcher(cur_node.label_G, cands_G[i], node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        if GM.subgraph_is_isomorphic(): 
        # if (GM.subgraph_is_isomorphic() and GM2.subgraph_is_isomorphic()) or GM2.subgraph_is_isomorphic(): 
            # draw_mol(cand_graph, 2000 + i, ["map_num", "bond_type", "color"], folder="extension")
            # draw_mol(cands_G[i], 2000 + i, ["map_num", "bond_type", "color"], folder="extension")
            print(2000 + i, True, label_amap)
            print(temp_global_amap)
            if i == 3: 
                total_label_amap.append(label_amap)
                chosen_amap = copy.deepcopy(temp_global_amap)
                chosen_graph = cand_graph.copy()

    #-------------------------------------------------------------------------------------------
    print("\Third iteration")
    temp_global_amap = chosen_amap
    starting_graph = chosen_graph
    print('c', temp_global_amap)

    cands_G, cands_G_amap = enum_attach_single_bond(starting_graph, children[2], temp_global_amap)
    candidate_graphs = [starting_graph.copy() for _ in cands_G_amap]
    temp_global_amaps = [copy.deepcopy(temp_global_amap) for _ in cands_G_amap]
    
    chosen_amap = None
    chosen_graph = None
    for i, label_amap in enumerate(cands_G_amap):
        cand_graph = candidate_graphs[i]
        temp_global_amap = temp_global_amaps[i]

        for nei_id,ctr_atom,nei_atom, ctr_id in label_amap: # nei nid, atom idx, atom idx
            if nei_id == fa_nid: continue

            try: temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]
            except:
                # add artificial nodes
                temp_global_amap[ctr_id][ctr_atom] = len(cand_graph.nodes)
                temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]

                node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
                cand_graph.add_node(len(cand_graph.nodes), **node_attr)
        
        cand_graph = attach_graphs(cand_graph, [children[2]], prev_nodes, temp_global_amap, print_out=True) #father is already attached
        # print(temp_global_amap)

        # GM = iso.GraphMatcher(cur_node.label_G, cand_graph, edge_match=ring_edge_equal_iso)
        # GM2 = iso.GraphMatcher(cur_node.label_G, cands_G[i], node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        # if GM.subgraph_is_isomorphic(): 
        #     draw_mol(cand_graph, 3000 + i, ["map_num", "bond_type", "color"], folder="extension")
        #     # draw_mol(cands_G[i], 3000 + i, ["map_num", "bond_type", "color"], folder="extension")
        #     print(3000 + i, True, label_amap)
        #     # if i == 3: 
        #     #     total_label_amap.append(label_amap)
        #     # chosen_amap = copy.deepcopy(temp_global_amap)
        #     #     starting_graph = cand_graph.copy()

        if i == 10:
            total_label_amap.append(label_amap)
            chosen_amap = copy.deepcopy(temp_global_amap)
            chosen_graph = cand_graph.copy()

    # ----------------------------------------------------------------------------------------------------------------
    print("\Fourth iteration")
    temp_global_amap = chosen_amap
    starting_graph = chosen_graph
    print('c', temp_global_amap)

    cands_G, cands_G_amap = enum_attach_double_bond2(starting_graph, children[3], temp_global_amap)
    candidate_graphs = [starting_graph.copy() for _ in cands_G_amap]
    temp_global_amaps = [copy.deepcopy(temp_global_amap) for _ in cands_G_amap]
    
    chosen_amap = None
    for i, label_amap in enumerate(cands_G_amap):
        cand_graph = candidate_graphs[i]
        temp_global_amap = temp_global_amaps[i]

        for nei_id,ctr_atom,nei_atom, ctr_id in label_amap: # nei nid, atom idx, atom idx
            if nei_id == fa_nid: continue

            try: temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]
            except:
                # add artificial nodes
                temp_global_amap[ctr_id][ctr_atom] = len(cand_graph.nodes)
                temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]

                node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
                cand_graph.add_node(len(cand_graph.nodes), **node_attr)
        
        cand_graph = attach_graphs(cand_graph, [children[3]], prev_nodes, temp_global_amap, print_out=True) #father is already attached
        # draw_mol(cand_graph, 4000 + i, ["map_num", "bond_type", "color"], folder="extension")

        # GM = iso.GraphMatcher(cur_node.label_G, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        # GM2 = iso.GraphMatcher(cur_node.label_G, cands_G[i], node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        # if GM.subgraph_is_isomorphic(): 
        #     draw_mol(cand_graph, 4000 + i, ["map_num", "bond_type", "color"], folder="extension")
        #     # draw_mol(cands_G[i], 3000 + i, ["map_num", "bond_type", "color"], folder="extension")
        #     print(4000 + i, True, label_amap)
        #     # if i == 3: 
        #     #     total_label_amap.append(label_amap)
        #     # chosen_amap = copy.deepcopy(temp_global_amap)
        #     #     starting_graph = cand_graph.copy()

        if i == 5:
            total_label_amap.append(label_amap)
            chosen_amap = copy.deepcopy(temp_global_amap)
            chosen_graph = cand_graph.copy()

    # raise
    # ----------------------------------------------------------------------------------------------------------------
    print("\Fifth iteration")
    temp_global_amap = chosen_amap
    starting_graph = chosen_graph
    print('c', temp_global_amap)

    # draw_mol(starting_graph, 4999, ["map_num", "bond_type", "color"], folder="extension")

    cands_G, cands_G_amap = enum_attach_double_bond2(starting_graph, children[4], temp_global_amap)
    candidate_graphs = [starting_graph.copy() for _ in cands_G_amap]
    temp_global_amaps = [copy.deepcopy(temp_global_amap) for _ in cands_G_amap]
    
    chosen_amap = None
    for i, label_amap in enumerate(cands_G_amap):
        cand_graph = candidate_graphs[i]
        temp_global_amap = temp_global_amaps[i]

        for nei_id,ctr_atom,nei_atom, ctr_id in label_amap: # nei nid, atom idx, atom idx
            if nei_id == fa_nid: continue

            try: temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]
            except:
                # add artificial nodes
                temp_global_amap[ctr_id][ctr_atom] = len(cand_graph.nodes)
                temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]

                node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
                cand_graph.add_node(len(cand_graph.nodes), **node_attr)
        
        cand_graph = attach_graphs(cand_graph, [children[4]], prev_nodes, temp_global_amap, print_out=True) #father is already attached
        # draw_mol(cand_graph, 5000 + i, ["map_num", "bond_type", "color"], folder="extension")
        # print(temp_global_amap[children[4].nid])

        GM = iso.GraphMatcher(cur_node.label_G, cand_graph, node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        GM2 = iso.GraphMatcher(cur_node.label_G, cands_G[i], node_match=node_equal_iso2, edge_match=ring_edge_equal_iso)
        if GM2.is_isomorphic(): 
            # draw_mol(cand_graph, 5000 + i, ["map_num", "bond_type", "color"], folder="extension")
            print(5000 + i, True, label_amap)
            if i == 5: 
                total_label_amap.append(label_amap)
                chosen_amap = copy.deepcopy(temp_global_amap)
                starting_graph = cand_graph.copy()

    # draw_mol(cur_node.label_G, 5300, ["map_num", "bond_type", "color"], folder="extension")


    # ---------------------------FINAL---------------------------------------------------------


            # temp_global_amap = copy.deepcopy(global_amap)
            # val_graph = cur_graph.copy()


            # temp_global_amap[neighbors[0].nid] = {node:node for node in neighbors[0].graph.nodes()} # allow global amap for loop to feed value in from ctr_node/cur_node
            # for nei_id,ctr_atom,nei_atom, ctr_id in concat_label_amap: # nei nid, atom idx, atom idx
            #     if nei_id == fa_nid: continue

            #     try: temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]
            #     except:
            #         # add artificial nodes
            #         temp_global_amap[ctr_id][ctr_atom] = len(val_graph.nodes)
            #         temp_global_amap[nei_id][nei_atom] = temp_global_amap[ctr_id][ctr_atom]

            #         node_attr = copy_node_attr(nid_graph[nei_id], nei_atom)
            #         val_graph.add_node(len(val_graph.nodes), **node_attr)
            
            # val_graph = attach_graphs(val_graph, neighbors, prev_nodes, temp_global_amap)
            # draw_mol(val_graph, 6000 + i, ["map_num", "bond_type", "color"], folder="extension")
            # concat_label_amap = [label_amap for label_amap_list in total_label_amap for label_amap in label_amap_list]