# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 01:20:09 2025

@author: Ann-Joelle
"""

import re
from pubchempy import get_compounds
import requests, time




def normalize_name(chemical_name, keep_apostrophe=True):
    """
    Cleans a chemical name by removing characters that may interfere with PubChem queries.
    
    Keeps only:
    - Letters (a–z, A–Z)
    - Numbers (0–9)
    - Commas, periods, hyphens, parentheses, whitespace
    - By default, keeps apostrophes (') for initial search, but allows variant where they are removed later.

    Example:
        normalize_name("1-(methylthio)-2-butanone!") → "1-(methylthio)-2-butanone"
        normalize_name("2-aminobenzoic@acid#")       → "2-aminobenzoicacid"

    Parameters:
        chemical_name (str): The raw chemical name string to clean.

    Returns:
        str: A cleaned name, suitable for PubChem lookup.
    """
    pattern = r"[^a-zA-Z0-9,'\-().\s]" if keep_apostrophe else r"[^a-zA-Z0-9,\-().\s]"
    return re.sub(pattern, "", chemical_name)




def smart_clean_name(name):
    """
    Performs advanced cleaning to boost lookup success for chemical names.
    
    Steps:
    - Replaces underscores with spaces
    - Adds spaces before common group suffixes (acid, ester, ether, etc.) if missing
    - Adds spaces after n-, iso-, sec-, tert-, t- prefixes
    - Collapses multiple spaces
    - Converts to Title Case (improves matches in some APIs)
    - Strips whitespace

    This should be used *after* normalized lookups have failed.

    Parameters:
        name (str): Chemical name to clean.

    Returns:
        str: Further cleaned/normalized name.
    """
    name = name.replace("_", " ")
    suffixes = [
        "acid", "ester", "ether", "alcohol", "anhydride", "acetate", "amine", "amide", "sulfide", "sulfone",
        "oxide", "formate", "propionate", "butyrate", "phenol", "nitrate", "chloride", "fluoride", "bromide", "iodide",
        "hydrazide", "phosphonate", "phosphinate", "phosphite", "phosphine", "isocyanate", "thiocyanate", "thioamide",
        "peroxide", "hydroperoxide", "nitrile", "carboxylate", "carboxylic acid", "carbamate", "benzoate", "benzamide",
        "oxalate", "tartrate", "thiophosphate", "arsonic acid", "arsenate", "arsenite"
    ]
    for suf in suffixes:
        name = re.sub(rf"(?<=[a-zA-Z0-9])({suf})", r" \1", name, flags=re.IGNORECASE)
    name = re.sub(r'\b([nst]?)-', r'\1 ', name)
    name = re.sub(r'\s+', ' ', name)
    name = name.title()
    name = name.strip()
    return name



def try_variants(name):
    """
    Generate prioritized name variants for robust chemical name-to-structure lookup.

    Variant generation logic:
      1. Includes the original name (with apostrophes, if present) as the first attempt.
      2. Adds a version with apostrophes removed, to match sources that ignore or do not support primes.
      3. Applies smart cleaning for common chemical name concatenations (see `smart_clean_name` docstring),
         and adds the cleaned result if not already present.
      4. Adds a version of the cleaned name with apostrophes removed.
      5. For suffixes likely to be run together in names (e.g., 'benzoicacid'), also tries the format
         with the suffix moved to the front (e.g., 'acid benzoic'), to improve matches for database quirks.
      6. If the name contains 'bis(...)', also tries with 'bis' removed for fallback matching.

    The function ensures that **each variant is unique and preserves original order of attempts**.

    Parameters:
        name (str): The raw or normalized chemical name.

    Returns:
        list[str]: A list of unique, prioritized name variants, suitable for sequential API queries.

    """
    variants = [name]
    if "'" in name:
        # Variant with apostrophes removed
        variants.append(name.replace("'", ""))
    cleaned = smart_clean_name(name)
    if cleaned not in variants:
        variants.append(cleaned)
    if "'" in cleaned:
        variants.append(cleaned.replace("'", ""))
    # (Add swap suffix and bis() logic as before)
    # ... rest unchanged ...
    swap_suffixes = ["acid", "ester", "ether", "alcohol", "anhydride", "acetate"]
    for suf in swap_suffixes:
        if cleaned.lower().endswith(suf):
            core = cleaned[:-len(suf)].strip()
            variant = f"{suf} {core}"
            if variant not in variants:
                variants.append(variant)
    if "bis(" in cleaned.lower():
        bisless = re.sub(r"bis\(([^)]+)\)", r"\1", cleaned, flags=re.IGNORECASE)
        if bisless not in variants:
            variants.append(bisless)
    # Ensure unique, preserve order
    return list(dict.fromkeys(variants))



def pubchem_rest(name, label, retries=2, delay=1):
    """
    Looks up a chemical's SMILES string using PubChem's REST API for a given name/variant.
    - Used for 'normalized name', 'cleaned name', or any 'variant'
    - The 'label' argument tags the search strategy for clarity.

    Parameters:
        name (str): Chemical name for lookup.
        label (str): Which name version is being tried (for status: 'normalized name', 'cleaned name', 'variant')
        retries (int): How many network attempts to try.
        delay (int or float): Delay between retries.

    Returns:
        tuple:
          - smiles (str or None): SMILES if found, else None
          - status (str): Describes method, label, and attempt used.
          - name_used (str): The name/variant used for this lookup.

    Status outcomes:
        - "method: PubChem REST (cleaned name), attempt: 1"
        - "method: PubChem REST (variant), attempt: 2"
        - "method: PubChem REST (normalized name) (not found), attempts: 2"
    """
    for attempt in range(1, retries + 1):
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES/TXT"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                return response.text.strip(), f"method: PubChem REST ({label}), attempt: {attempt}", name
        except Exception as e:
            print(f"REST attempt {attempt} failed for {name}: {e}")
        time.sleep(delay)
    return None, f"method: PubChem REST ({label}) (not found), attempts: {retries}", name


def pubchempy_lookup(name, label):
    """
    Looks up a chemical's SMILES string using PubChemPy's synonym search for a given name/variant.
    - Distinguishes between exact match and synonym match.
    - Used for 'normalized name', 'cleaned name', or any 'variant'.

    Parameters:
        name (str): Chemical name for lookup.
        label (str): Which name version is being tried.

    Returns:
        tuple:
          - smiles (str or None): SMILES if found, else None
          - status (str): Method and match type, e.g. 'method: PubChemPy exact (cleaned name), attempt: 1'
          - used_synonym (str or None): The synonym that matched, or the input name for exact matches.

    Status outcomes:
        - "method: PubChemPy exact (cleaned name), attempt: 1"
        - "method: PubChemPy synonym (variant), attempt: 1"
        - "method: PubChemPy (normalized name) (not found)"
    """
    try:
        compounds = get_compounds(name, 'name')
        if compounds and compounds[0].canonical_smiles:
            if compounds[0].synonyms and name not in compounds[0].synonyms:
                return compounds[0].canonical_smiles, f"method: PubChemPy synonym ({label}), attempt: 1", compounds[0].synonyms[0]
            else:
                return compounds[0].canonical_smiles, f"method: PubChemPy exact ({label}), attempt: 1", name
    except Exception as e:
        print(f"PubChemPy failed for {name}: {e}")
    return None, f"method: PubChemPy ({label}) (not found)", None


def opsin_lookup(name, label):
    """
    Looks up a chemical's SMILES string using OPSIN (IUPAC parser).
    - Used only as a last resort, for all variants.

    Parameters:
        name (str): Chemical name for lookup.
        label (str): Which variant is being tried.

    Returns:
        tuple:
          - smiles (str or None): SMILES if found, else None
          - status (str): Always "method: OPSIN (variant), attempt: 1" (if found)
          - name (str): The name/variant used for this lookup.

    Status outcomes:
        - "method: OPSIN (variant), attempt: 1"
        - "method: OPSIN (variant) (not found)"
    """
    try:
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.smi"
        response = requests.get(opsin_url, timeout=10)
        if response.status_code == 200 and response.text.strip():
            return response.text.strip(), f"method: OPSIN ({label}), attempt: 1", name
    except Exception as e:
        print(f"OPSIN failed for {name}: {e}")
    return None, f"method: OPSIN ({label}) (not found)", name



def fetch_best_smiles(chemical_name, retries=3, delay=1):
    """
    Master function for robust, prioritized chemical name-to-SMILES conversion.

    It tries, in order:
      1. PubChem REST with normalized name
      2. PubChemPy with normalized name
      3. PubChem REST with cleaned name
      4. PubChemPy with cleaned name
      5. PubChem REST for each try_variants 
      6. PubChemPy for each try_variants 
      7. OPSIN for each try_variants 

    For every lookup, it returns:
      - The SMILES string if found (or None)
      - The method/strategy, API/network attempt, and variant label in the status column
      - The name or synonym/variant that matched

    This makes it fully auditable, transparent, and ready for large-scale, traceable chemical data processing.

    Example status outputs:
      - method: PubChem REST (normalized name), attempt: 1
      - method: PubChemPy synonym (cleaned name), attempt: 1
      - method: OPSIN (variant), attempt: 1
      - method: PubChem REST (variant) (not found), attempts: 2

    Parameters:
        chemical_name (str): Chemical name for lookup.
        retries (int): Number of retries for PubChem REST requests.
        delay (int or float): Delay in seconds between retries.

    Returns:
        tuple:
          - smiles (str or None): SMILES if found, else None
          - status (str): Detailed trace of method, strategy, attempt, and variant label
          - name/synonym/variant used (str): Name that produced the match, or "" if not found
    """
    norm_name = normalize_name(chemical_name)

    # 1. PubChem REST with normalized name
    smiles, status, name_used = pubchem_rest(norm_name, label="original name", retries=retries, delay=delay)
    if smiles:
        return smiles, status, name_used

    # 2. PubChemPy with normalized name
    smiles, status, name_used = pubchempy_lookup(norm_name, label="original name")
    if smiles:
        return smiles, status, name_used

    # 3. PubChem REST with cleaned name
    cleaned = smart_clean_name(norm_name)
    if cleaned != norm_name:
        smiles, status, name_used = pubchem_rest(cleaned, label="cleaned name", retries=retries, delay=delay)
        if smiles:
            return smiles, status, name_used

        # 4. PubChemPy with cleaned name
        smiles, status, name_used = pubchempy_lookup(cleaned, label="cleaned name")
        if smiles:
            return smiles, status, name_used

    # 5/6/7. For all variants: PubChem REST, PubChemPy, OPSIN
    variants = try_variants(norm_name)
    for idx, variant in enumerate(variants):
        # skip the first (already tried as "original name")
        # skip the second (if same as cleaned)
        if (idx == 0) or (idx == 1 and cleaned != norm_name):
            continue
        label = "variant"
        # Try PubChem REST
        smiles, status, name_used = pubchem_rest(variant, label=label, retries=retries, delay=delay)
        if smiles:
            return smiles, status, name_used

        # Try PubChemPy
        smiles, status, name_used = pubchempy_lookup(variant, label=label)
        if smiles:
            return smiles, status, name_used

    # 8. OPSIN for all variants (with correct label)
    for idx, variant in enumerate(variants):
        if idx == 0:
            label = "original name"
        elif idx == 1 and cleaned != norm_name:
            label = "cleaned name"
        else:
            label = "variant"
        smiles, status, name_used = opsin_lookup(variant, label=label)
        if smiles:
            return smiles, status, name_used

    return None, "Not found with any method", ""


