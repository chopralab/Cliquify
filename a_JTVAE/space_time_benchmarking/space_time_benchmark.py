import sys
sys.path.append("../")
import rdkit.Chem as Chem
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

from jtnn.mol_tree import *
from jtnn.chemutils import *

import timeit
from multiprocessing import Pool

if __name__ == "__main__":

    filename = "fail_mol_list(multiprocessing).txt"

    with open("../../zdata/zinc/all.txt") as f:
        smiles_list = f.readlines()

    def decode_test():
        
        wrong = 0
        desc_count = 0
        for tot,s in enumerate(smiles_list[:100_000]):


            tree = MolTree(s)
            tree.recover()

            label_idx = 0
            try: tree.nodes[label_idx]
            except: continue
            cur_mol = copy_edit_mol(tree.nodes[label_idx].mol)
            global_amap = [{}] + [{} for node in tree.nodes]
            global_amap[label_idx + 1] = {atom.GetIdx():atom.GetIdx() for atom in cur_mol.GetAtoms()}

            dfs_assemble(cur_mol, global_amap, [], tree.nodes[label_idx], None)

            cur_mol = cur_mol.GetMol()
            cur_mol = Chem.MolFromSmiles(Chem.MolToSmiles(cur_mol))
            if not cur_mol:
                # print(gold_smiles)
                # print(None)
                # print(wrong, tot + 1, "idx:", tot)
                wrong += 1
                print(tot, wrong, gold_smiles, "decompose")
                # with open("fail_mol_list.txt", "a") as myfile:
                #     # print(gold_smiles, tot)
                #     myfile.writelines("{},{},decompose\n".format(gold_smiles, tot))
                
                continue

            set_atommap(cur_mol)
            dec_smiles = Chem.MolToSmiles(cur_mol)

            # gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(s))
            gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=False)
            
            # TEST FOR CHIRALITY DURING DECODING
            # if s != dec_smiles:
            #     print(tot, wrong, gold_smiles, "chirality")

            if gold_smiles != dec_smiles:
                # print(gold_smiles)
                # print(dec_smiles)
                wrong += 1
                print(tot, wrong, gold_smiles, "reconstruct")
                # with open("fail_mol_list.txt", "a") as myfile:
                #     # print(gold_smiles, tot)
                #     myfile.writelines("{},{},reconstruct\n".format(gold_smiles, tot))

                # if len(gold_smiles) == len(dec_smiles):
                #     desc_count += 1
                #     # print(tot, gold_smiles, wrong, "desc prob?:", len(gold_smiles) == len(dec_smiles), desc_count)
                    
                #     with open("descendant_orient_awareness.txt", "a") as myfile:
                #         # print(gold_smiles, tot)
                #         myfile.writelines("{},{},{}\n".format(gold_smiles, dec_smiles,  tot))
    
    def decode_test2(idx):

        s = smiles_list[idx]

        try:
            tree = MolTree(s)
        except:
            return "{},{},Decompose".format(s, idx)
        tree.recover()

        label_idx = 0
        try: tree.nodes[label_idx]
        except: return
        cur_mol = copy_edit_mol(tree.nodes[label_idx].mol)
        global_amap = [{}] + [{} for node in tree.nodes]
        global_amap[label_idx + 1] = {atom.GetIdx():atom.GetIdx() for atom in cur_mol.GetAtoms()}

        dfs_assemble(cur_mol, global_amap, [], tree.nodes[label_idx], None)

        cur_mol = cur_mol.GetMol()
        cur_mol = Chem.MolFromSmiles(Chem.MolToSmiles(cur_mol))
        if not cur_mol:
            return "{},{},Reconstruct".format(s, idx)
            
        set_atommap(cur_mol)
        dec_smiles = Chem.MolToSmiles(cur_mol)

        gold_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=False)
        
        if gold_smiles != dec_smiles:
            return "{},{},Reconstruct".format(gold_smiles, idx)
            # print(idx, gold_smiles, "reconstruct")
        return

    def multiprocess_decode():
        pool = Pool()
        fail_list = pool.map(decode_test2, range(len(smiles_list)))
        return fail_list

    def write_file(smiles_list):
        smiles_list = [smiles for smiles in smiles_list if smiles]
        with open(filename, "a") as myfile:
            myfile.writelines(smiles_list)

    # print(timeit.Timer(decode_test).timeit(number=1))
    print(timeit.Timer(multiprocess_decode).timeit(number=1))

