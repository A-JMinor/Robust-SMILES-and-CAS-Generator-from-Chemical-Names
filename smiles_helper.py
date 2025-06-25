# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 01:20:09 2025

@author: Ann-Joelle
"""

import re
import requests, time
import pubchempy as pcp




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



def comma_space_variants(name):
    """
    Returns unique variants of a chemical name with and without spaces after commas.

    Useful for database lookups where "A,B" and "A, B" may be treated differently.

    Parameters:
        name (str): Chemical name.

    Returns:
        list[str]: [name_with_spaces, name_without_spaces] (duplicates removed).

    Example:
        >>> comma_space_variants("2-pyridinamine,1-oxide")
        ['2-pyridinamine, 1-oxide', '2-pyridinamine,1-oxide']
    """
    with_space = re.sub(r",\s*", ", ", name)
    without_space = re.sub(r",\s*", ",", name)
    variants = [with_space]
    if without_space != with_space:
        variants.append(without_space)
    return list(dict.fromkeys(variants))



def smart_clean_name(name, title_case=False):
    """
    Minimal cleaning for chemical names without altering their structure.

    - Replaces underscores with spaces.
    - Collapses multiple spaces and trims whitespace.
    - (Optional) Converts to Title Case if title_case=True.
    - Does not split glued suffixes or do aggressive normalization.

    Parameters:
        name (str): The raw chemical name.
        title_case (bool): If True, applies .title() to the result.

    Returns:
        str: The minimally cleaned name.
    
    Example:
        >>> smart_clean_name("ethyl_acetate")
        'ethyl acetate'
        >>> smart_clean_name("ethyl_acetate", title_case=True)
        'Ethyl Acetate'
    """
        
    name = name.replace("_", " ")
    name = re.sub(r'\s+', ' ', name)
    name = name.strip()
    if title_case:
        name = name.title()
    return name



def split_suffix_variants(name):
    """
    Generates all variants of a chemical name with recursively split suffixes (no extra cleaning).
    Example:
        split_suffix_variants("1-piperazinecarboxylicacidethylester")
        -> [
            '1-piperazinecarboxylicacidethylester',
            '1-piperazinecarboxylicacidethyl Ester',
            '1-piperazinecarboxylic Acidethyl Ester',
            '1-piperazine Carboxylic Acidethyl Ester',
            '1-piperazine Carboxylic Acid Ethyl Ester'
        ]
    Parameters:
        name (str): The chemical name to split.
    Returns:
        list[str]: All split variants, from least to most split.
    """
    import re

    suffixes = [
        "carboxylic acid", "benzenesulfonic acid", "arsonic acid", "sulfonic acid", "phosphonic acid", "phosphinic acid",
        "acid", "ester", "ether", "alcohol", "anhydride", "acetate", "amine", "amide", "sulfide", "sulfone",
        "oxide", "formate", "propionate", "butyrate", "phenol", "nitrate", "chloride", "fluoride", "bromide", "iodide",
        "hydrazide", "phosphonate", "phosphinate", "phosphite", "phosphine", "isocyanate", "thiocyanate", "thioamide",
        "peroxide", "hydroperoxide", "nitrile", "carboxylate", "carbamate", "benzoate", "benzamide",
        "oxalate", "tartrate", "thiophosphate", "arsenate", "arsenite", "ethyl", "methyl", "propyl", "butyl", "phenyl"
    ]
    suffixes = sorted(suffixes, key=len, reverse=True)

    def recursive_splits(s, start=0, seen=None):
        if seen is None:
            seen = set()
        result = set([s])
        for suf in suffixes:
            for m in re.finditer(rf'(?i)([a-zA-Z0-9\-])({suf})', s):
                idx = m.start(2)
                if s[idx-1] not in (' ', '-'):
                    before = s[:idx]
                    after = s[idx:]
                    new = before + " " + after
                    new = re.sub(r'\s+', ' ', new).strip()
                    if new not in seen:
                        seen.add(new)
                        result |= recursive_splits(new, idx+1, seen)
        return result

    variants = list(recursive_splits(name))
    variants = sorted(variants, key=lambda x: (x.count(' '), len(x)))
    return variants



def try_variants(name):
    
    """
    Generate prioritized, unique chemical name variants for robust structure lookup (e.g., in PubChem).

    This function produces likely name variants to maximize SMILES retrieval rates from structure databases,
    compensating for differences in formatting, punctuation, and group-core ordering.

    Variant generation strategy:
        1. Original input name.
        2. Version with all apostrophes (') removed.
        3. Minimally cleaned name (spaces, underscores, no case change).
        4. Cleaned name with apostrophes removed.
        5. Title-case cleaned name and its de-apostrophized version.
        6. For select group suffixes (see swap_cores, e.g. "benzene", "oxide", etc.) at the end:
            a. Suffix moved to the front with a comma.
            b. Suffix moved to the front with a dash.
        7. For names with "bis(...)", tries a dash-connected variant.
        8. All outputs are unique and preserve the order generated.

    Parameters:
        name (str): The raw or normalized chemical name.

    Returns:
        list[str]: Ordered, unique name variants for sequential lookup.

    Examples:
        >>> try_variants("1,1'-sulfonylbis(4-methyl-)benzene")
        [
          "1,1'-sulfonylbis(4-methyl-)benzene",
          "1,1-sulfonylbis(4-methyl-)benzene",
          "1,1'-Sulfonylbis(4-Methyl-)Benzene",
          "1,1-Sulfonylbis(4-Methyl-)Benzene",
          "Benzene, 1,1'-Sulfonylbis(4-Methyl-)",
          "Benzene-1,1'-Sulfonylbis(4-Methyl-)",
          "1,1'-Sulfonyl-4-Methyl-Benzene"
        ]

        >>> try_variants("2-pyridinamine, 1-oxide")
        [
          "2-pyridinamine, 1-oxide",
          "2-Pyridinamine, 1-Oxide",
          "1-oxide, 2-pyridinamine",
          "1-oxide-2-pyridinamine"
        ]

    Notes:
        - The function preserves dashes and parentheses, as these often matter for correct structure assignment.
        - The swap-core logic covers common inorganic/organic groups frequently swapped in PubChem synonym lists.
        - The bis(...) logic handles "bis(substituent)" names for linked/bridged compounds.
        - Numeric locants (e.g., "1-", "2-") are preserved except as moved by swap-core logic.

    Usage:
        Use the returned list as a sequence of fallback queries for chemical structure lookup.
    """
    
    
    variants = [name]
    if "'" in name:
        variants.append(name.replace("'", ""))
    cleaned = smart_clean_name(name)
    if cleaned not in variants:
        variants.append(cleaned)
    if "'" in cleaned:
        variants.append(cleaned.replace("'", ""))

    # Swapping core at the end, keeping dashes!
    swap_cores = [
        "oxide", "chloride", "bromide", "fluoride", "iodide", 
        "hydrochloride", "hydrobromide",
        "sulfate", "phosphate", "nitrate", "acetate", "formate", "benzoate",
        "sulfonate",
        "carboxylic acid", "anhydride",
        "amide", "ester", "ether", "alcohol",
        "benzene", "phenol", "pyridine", "pyrimidine", "pyrrole", "quinoline", "aniline"
    ]


    for core in swap_cores:
        pattern = rf"^(.*?)[, ]*((?:\d+-|N-)?{core})(-?)$"
        m = re.match(pattern, cleaned, re.IGNORECASE)
        if m:
            rest = m.group(1).strip(" ,-")
            core_part = m.group(2)
            # dash = m.group(3)   # <--- NICHT mehr verwenden!
            swapped1 = f"{core_part.capitalize()}, {rest}"
            # Füge mit Bindestrich (ohne Komma) nur, wenn beides nicht leer ist
            if rest:
                swapped2 = f"{core_part.capitalize()}-{rest}"
                if swapped2 not in variants and swapped2.strip():
                    variants.append(swapped2)
            if swapped1 not in variants and swapped1.strip():
                variants.append(swapped1)
    

    # Bis-Logik wie gehabt
    if "bis(" in cleaned.lower():
        # Füge ein Minus nach "Sulfonyl" hinzu, falls nach der Klammer (bei bis()) keins steht
        bisless = re.sub(
            r"bis\(([^)]+)\)",
            lambda m: ("-" if not m.group(1).startswith("-") else "") + m.group(1).rstrip("-") + "-",
            cleaned,
            flags=re.IGNORECASE
        )

        if bisless not in variants:
            variants.append(bisless)

    cleaned = smart_clean_name(name)
    cleaned_title = smart_clean_name(name, title_case=True)
    if cleaned not in variants:
        variants.append(cleaned)
    if cleaned_title not in variants:
        variants.append(cleaned_title)
    if "'" in cleaned:
        variants.append(cleaned.replace("'", ""))
    if "'" in cleaned_title:
        variants.append(cleaned_title.replace("'", ""))

            
    return list(dict.fromkeys(variants))



def pubchem_rest(name, label, retries=3, delay=1):
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
          - status (str): If multi-structure, prefix with 'Probably wrong;'
          - name (str): The name/variant used for this lookup.

    Status outcomes:
        - "method: OPSIN (variant), attempt: 1"
        - "Probably wrong; method: OPSIN (variant), attempt: 1" (if SMILES contains '.')
        - "method: OPSIN (variant) (not found)"
    """
    try:
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.smi"
        response = requests.get(opsin_url, timeout=10)
        if response.status_code == 200 and response.text.strip():
            smiles = response.text.strip()
            status = f"method: OPSIN ({label}), attempt: 1"
            # Check for multi-structure ('.')
            if '.' in smiles:
                status = "Probably wrong; " + status
            return smiles, status, name
    except Exception as e:
        print(f"OPSIN failed for {name}: {e}")
    return None, f"method: OPSIN ({label}) (not found)", name




def fetch_best_smiles(chemical_name, retries=3, delay=0.5):
    """
    Robust, prioritized chemical name-to-SMILES resolver using PubChem and OPSIN, designed for challenging or inconsistent chemical name data.

    The function systematically queries both PubChem and OPSIN using a wide range of name cleaning and variant-generation strategies,
    attempting to maximize the chance of successful resolution, even for "messy" or unconventional chemical names.

    **Lookup and variant order:**
      1. **PubChem** (REST API and PubChemPy) with the normalized name (removes problematic symbols, standardizes allowed characters).
      2. **PubChem** (REST and Py) for all variants with and without spaces after commas (some databases store both styles).
      3. **PubChem** (REST and Py) with the minimally cleaned name (removes underscores, extra spaces).
      4. **PubChem** (REST and Py) for all suffix-split variants (aggressively inserts spaces before known chemical suffixes; applied to each comma-space variant).
      5. **PubChem** (REST and Py) for all name variants generated by `try_variants` (includes swapped group cores, de-apostrophized, bis logic, and title-case forms).
      6. **OPSIN** with the normalized name.
      7. **OPSIN** for all comma-space variants.
      8. **OPSIN** with the minimally cleaned name.
      9. **OPSIN** for all suffix-split variants (for each comma-space variant).
     10. **OPSIN** for all `try_variants` of the normalized name.

    For each query, the function returns as soon as a result is found, providing:
        - The first matching SMILES (as a string), or None if not found.
        - A detailed status string (including method, variant type, attempt, and a warning if the result is likely incorrect, e.g., if OPSIN returns a SMILES with a ".").
        - The exact name/variant used for the successful query.

    **Parameters:**
        chemical_name (str): The chemical name to resolve.
        retries (int): Number of retries for PubChem REST queries.
        delay (float): Delay (seconds) between PubChem REST retries.

    **Returns:**
        tuple:
          - smiles (str or None): Canonical SMILES if found, else None.
          - status (str): Trace of method, strategy, attempt, and variant label. If the returned SMILES contains a period (".") in an OPSIN result, "Probably wrong;" is prepended to the status.
          - name_used (str): The name/variant used for the successful query, or "" if not found.

    **Notes:**
    - The function tries both PubChem REST (direct HTTP) and PubChemPy (Python library), as each can resolve different synonym styles.
    - Comma-space variants address databases where comma and space are inconsistently handled in names.
    - Suffix splitting increases robustness for names with glued suffixes, e.g., "benzoicacid" → "benzoic acid".
    - The `try_variants` function tries various fallback name manipulations, including swapping group cores, bis logic, and title-casing.
    - OPSIN is always used last, and a warning is added if the returned SMILES is a disconnected structure (contains ".").


    **Example usage:**
        smiles, status, name_used = fetch_best_smiles("n-butylbenzenehydroperoxide")
        print(smiles, status, name_used)
    """

    # 1. PubChem with normalized name (REST and Py)
    norm_name = normalize_name(chemical_name)
    for pubchem_func, label in [(pubchem_rest, "normalized name"), (pubchempy_lookup, "normalized name")]:
        smiles, status, name_used = pubchem_func(norm_name, label=label, retries=retries, delay=delay) if pubchem_func == pubchem_rest else pubchem_func(norm_name, label=label)
        if smiles:
            return smiles, status, name_used

    # 2. PubChem with comma_space_variants
    for variant in comma_space_variants(norm_name):
        for pubchem_func in [pubchem_rest, pubchempy_lookup]:
            smiles, status, name_used = pubchem_func(variant, label="comma_space_variant", retries=retries, delay=delay) if pubchem_func == pubchem_rest else pubchem_func(variant, label="comma_space_variant")
            if smiles:
                return smiles, status, name_used

    # 3. PubChem with smart_clean_name
    cleaned = smart_clean_name(norm_name)
    for pubchem_func in [pubchem_rest, pubchempy_lookup]:
        smiles, status, name_used = pubchem_func(cleaned, label="cleaned", retries=retries, delay=delay) if pubchem_func == pubchem_rest else pubchem_func(cleaned, label="cleaned")
        if smiles:
            return smiles, status, name_used

    # 4. PubChem with all split_suffix_variants (for each comma-space variant)
    for cvariant in comma_space_variants(norm_name):
        for split in split_suffix_variants(cvariant):
            for pubchem_func in [pubchem_rest, pubchempy_lookup]:
                smiles, status, name_used = pubchem_func(split, label="split_suffix", retries=retries, delay=delay) if pubchem_func == pubchem_rest else pubchem_func(split, label="split_suffix")
                if smiles:
                    return smiles, status, name_used

    # 5. PubChem with try_variants
    for variant in try_variants(norm_name):
        for pubchem_func in [pubchem_rest, pubchempy_lookup]:
            smiles, status, name_used = pubchem_func(variant, label="try_variant", retries=retries, delay=delay) if pubchem_func == pubchem_rest else pubchem_func(variant, label="try_variant")
            if smiles:
                return smiles, status, name_used

    # 6. OPSIN with normalized name
    smiles, status, name_used = opsin_lookup(norm_name, label="normalized name")
    if smiles:
        if '.' in smiles:
            status = "Probably wrong; " + status
        return smiles, status, name_used

    # 7. OPSIN with comma_space_variants
    for variant in comma_space_variants(norm_name):
        smiles, status, name_used = opsin_lookup(variant, label="comma_space_variant")
        if smiles:
            if '.' in smiles:
                status = "Probably wrong; " + status
            return smiles, status, name_used

    # 8. OPSIN with smart_clean_name
    smiles, status, name_used = opsin_lookup(cleaned, label="cleaned")
    if smiles:
        if '.' in smiles:
            status = "Probably wrong; " + status
        return smiles, status, name_used

    # 9. OPSIN with all split_suffix_variants (for each comma-space variant)
    for cvariant in comma_space_variants(norm_name):
        for split in split_suffix_variants(cvariant):
            smiles, status, name_used = opsin_lookup(split, label="split_suffix")
            if smiles:
                if '.' in smiles:
                    status = "Probably wrong; " + status
                return smiles, status, name_used

    # 10. OPSIN with try_variants
    for variant in try_variants(norm_name):
        smiles, status, name_used = opsin_lookup(variant, label="try_variant")
        if smiles:
            if '.' in smiles:
                status = "Probably wrong; " + status
            return smiles, status, name_used

    return None, "Not found with any method", ""


# === 🔎 Get CAS numbers via PubChem ===

def smiles_to_cas(smiles):
    """
    Given a SMILES string, query PubChem to extract CAS numbers from synonyms.

    Parameters:
        smiles (str): A chemical's SMILES string.

    Returns:
        list of str: A list of CAS Registry Numbers found in the synonyms,
                     or [''] if no match is found or if an error occurs.
    """
    cas_pattern = re.compile(r"(\d{2,7})-(\d{2})-(\d)")
    try:
        compound = pcp.get_compounds(smiles, 'smiles')[0]
        synonyms = compound.synonyms
        matches = [m for s in synonyms for m in cas_pattern.findall(s)]
        return ['-'.join(match) for match in matches]
    except Exception as e:
        print(f"❌ Error for {smiles}: {e}")
        return ['']
