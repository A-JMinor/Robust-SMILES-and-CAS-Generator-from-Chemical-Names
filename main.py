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
        "adamantane@carboxylic_acid",            
        "diethylpropanedioicacid",                            
        "ethylbenzoate",                         
        "1,1'-thiobisethene",                             
        "2,4,6-trinitrotoluene",                 
        "p-methoxyphenylaceticacid",                             
        "benzoicacidanhydride",                  
        "dimethyl-2,4,5-trichlorophenylthiophosphate",                                                                                      
        "cis-1,2,3,6-tetrahydrophthalicanhydride",                   
        "hexanedioicacid-bis(2-ethylhexyl)ester", 
        "1,1-ethanediamine",                     
        "2-pyridinamine,1-oxide",             
        "1,1'-sulfonylbis(4-methyl-)benzene",   
        "1-piperazinecarboxylicacidethylester"
        # 
    ]
    for chem in tqdm(test_chemicals, desc="SMILES Lookup"):
        smiles, status, matched_name = fetch_best_smiles(chem)
        print(f"\nChemical: {chem}")
        print(f"  SMILES: {smiles}")
        print(f"  Status: {status}")
        print(f"  Name/Synonym/Variant Used: {matched_name}")

