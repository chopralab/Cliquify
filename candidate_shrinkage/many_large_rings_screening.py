
import sys
sys.path.append('../')
import networkx as nx

from rdkit import Chem

with open("../zinc/all.txt") as f:
    smiles_list = f.readlines()

def bridge_checking(cliques):
    for i in range(len(cliques)):
        for j in range(i + 1, len(cliques)):
            if cliques[i] == cliques[j]: continue
            inter = set(cliques[i]) & set(cliques[j])
            if len(inter) > 2:
                return True
    return False

def many_large_rings():
    with open("many_large_rings_increment(idx).txt", "a") as myfile:
        for num_rings in range(2, 10):
            count = 0
            for smiles in smiles_list:
                mol = Chem.MolFromSmiles(smiles)

                if count == 500: break

                ssr = [list(x) for x in Chem.GetSymmSSSR(mol)]
                if bridge_checking(ssr): continue

                large_ssr = [list(x) for x in Chem.GetSymmSSSR(mol) if len(list(x)) >= 5]


                # if len(ssr) > 5:
                if len(large_ssr) == num_rings:
                    # print(smiles[:-1])
                    # myfile.writelines("{}".format(smiles))
                    myfile.writelines("{},{}\n".format(smiles.strip(), num_rings))
                    count += 1
many_large_rings()

def bridge_rings():            
    with open("bridge_rings.txt", "a") as myfile:
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)

            # if count == 1000: break

            ssr = [list(x) for x in Chem.GetSymmSSSR(mol)]
            if bridge_checking(ssr):
                myfile.writelines("{}".format(smiles))
                count += 1

    print(count)