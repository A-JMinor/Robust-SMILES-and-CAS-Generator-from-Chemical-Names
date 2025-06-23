# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 01:20:09 2025

@author: Ann-Joelle
"""

from tqdm import tqdm
from smiles_helper import fetch_best_smiles



# Example/test
if __name__ == "__main__":
    test_chemicals = [
        "adamantane@carboxylic_acid",             # underscore, run-together
        "diethylpropanedioicacid",               # run-together, dicarboxylic
        "n-butylbenzenehydroperoxide",           # prefix, peroxide
        "cis-9-octadecenoicacid",                # cis, dashes, run-together
        "ethylbenzoate",                         # easy case
        "1,1'-thiobisethene",                    # prime, special char
        "alpha-aminoisobutanoicacid",            # greek prefix, run-together
        "2,4,6-trinitrotoluene",                 # classic explosive (TNT)
        "p-methoxyphenylaceticacid",             # para, run-together
        "tert-butylformate",                     # tert prefix
        "benzoicacidanhydride",                  # anhydride
        "dimethyl-2,4,5-trichlorophenylthiophosphate",  # long, complex
        "trans-crotonic-acid",                   # trans, hyphens, acid
        "l-pyroglutamic_acid",                   # stereochemistry, underscore
        "isoamylchloride",                       # iso prefix
        "n'-ethylmethanediamine",                # n' with prime
        "o-methioaniline",                       # ortho prefix, spelling variant
        "cis-1,2,3,6-tetrahydrophthalicanhydride", # complex ring system
        "d-glucuronic_acid",                     # d- prefix, underscore
        "hexanedioicacid-bis(2-ethylhexyl)ester" # bis group, dicarboxylic ester
    ]
    for chem in tqdm(test_chemicals, desc="SMILES Lookup"):
        smiles, status, matched_name = fetch_best_smiles(chem)
        print(f"\nChemical: {chem}")
        print(f"  SMILES: {smiles}")
        print(f"  Status: {status}")
        print(f"  Name/Synonym/Variant Used: {matched_name}")

