# Robust-SMILES-Generator-from-Chemical-Names

In practice, chemical names from laboratory inventory lists, COSMO-RS solvent databases, procurement catalogs, or old experimental records are often inconsistent, concatenated, have typos, ambiguous punctuation, non-IUPAC spellings, or legacy formats entered by different people over years. Most tools or APIs expect well-formed, IUPAC-compliant names or common synonyms. This package is built for robust, automated generation of canonical SMILES from chemical names, even with inconsistent formatting, typos, and non-standard nomenclature. The workflow systematically cleans names, generates intelligent name variants, and performs multi-tiered database lookups (PubChem REST, PubChemPy synonym, and OPSIN) to maximize your chances of getting the correct SMILES. Results are fully traceable, so you always know which variant matched and how.

---

## Features

- **Automated structure assignment** from raw or inconsistent chemical names.
- **Intelligent name cleaning and normalization:** Handles underscores, concatenated names, apostrophes, misplaced/extra spaces, ambiguous commas, and other data quirks.
- **Combinatorial variant generation:** Systematically increases match rate by considering multiple spellings, punctuation, core-group swaps, and suffix splitting for each name.
- **Multi-tiered lookup:** Queries multiple authoritative sources (PubChem REST, PubChemPy, and OPSIN) in order of reliability.
- **Comprehensive traceability:** Records the exact method, variant, and lookup status for each assignment.

---

## How It Works

### 1. Name Extraction and Preparation

- Cleans each name, preserving original punctuation and spacing for initial queries.
- Generates prioritized lists of **name variants**:
    - Removing/adding spaces after commas (`comma_space_variants`)
    - Minimal “smart” cleaning (removes underscores, extra spaces, optional Title Case)
    - Progressive suffix splitting (recursively separates glued functional groups)
    - Core-group swaps (moves “oxide”, “benzene”, etc. to the front if present)
    - Handles apostrophes and case normalization

### 2. Multi-Tiered Structure Assignment

- **Tier 1: PubChem REST API and PubChemPy Synonym Search**  
  Direct query for canonical SMILES using each variant and any registered synonym (alternate name, IUPAC, trade name) matches in PubChem’s database. If found, the pipeline records the result and stops further lookups for that name.
- **Tier 2: OPSIN (Open Parser for Systematic IUPAC Nomenclature)**  
  For unresolved names, attempts systematic IUPAC parsing via OPSIN, a specialized open-source tool for structure-from-name conversion.

### 3. Metadata and Traceability

- For each chemical, the following are recorded:
  - **Original Name**
  - **Assigned SMILES**
  - **SMILES_Status** (describes the method, variant, and number of attempts)
  - **SMILES_NameUsed** (the exact name or synonym that yielded the match)
- If all attempts fail, the status is labeled as “not found with any method.”

---

