with open("all.txt", "r") as f:
    smiles = f.readline()
    while smiles:
        if "." in smiles: 
            smiles = f.readline()
            continue
        with open("all2.txt", "a") as af:
            af.write(smiles)
        smiles = f.readline()