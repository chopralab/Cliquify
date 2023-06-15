from rdkit import Chem

# Specify the path to the SMILES file
file_path = "COCONUT_DB.smi"

# Read the SMILES file
suppl = Chem.SmilesMolSupplier(file_path)

# Iterate over the molecules in the file
count = 0
with open("all.txt", "w") as f:
    for molecule in suppl:
        # Perform desired operations with each molecule
        # For example, you can access properties or perform calculations
        try:
            smiles = "{}\n".format(Chem.MolToSmiles(molecule))
            f.write(smiles)
        except: count += 1

print("Fail Mol", count)