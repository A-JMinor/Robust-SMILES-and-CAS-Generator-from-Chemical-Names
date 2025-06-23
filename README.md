# Robust-SMILES-Generator-from-Chemical-Names
Robust, automated generation of canonical SMILES from chemical names, even with inconsistent formatting, typos, and non-standard nomenclature. The workflow systematically cleans names, generates intelligent variants, and performs multi-tiered database lookups (PubChem REST, PubChemPy synonym, and OPSIN).

---

## Features
- **Automated structure assignment** from raw or inconsistent chemical names.
- **Intelligent name cleaning and normalization:** Handles underscores, concatenated names, apostrophes, and other data quirks.
- **Systematic variant generation:** Increases match rate by considering multiple spellings and arrangements of each name.
- **Multi-tiered lookup:** Queries multiple authoritative sources (PubChem REST, PubChemPy, and OPSIN) in order of reliability.
- **Comprehensive traceability:** Records the exact method, variant, and lookup status for each assignment.

---

## How It Works

### 1. Name Extraction and Preparation
- Cleans each name, preserving original punctuation for the initial attempt.
- Generates a prioritized list of name variants by removing/adding spaces, handling apostrophes, and considering common functional group positions and naming conventions.

### 2. Multi-Tiered Structure Assignment
- **Tier 1: PubChem REST API**  
  Direct query for canonical SMILES using each variant. If found, the pipeline records the result and stops further lookups for that name.
- **Tier 2: PubChemPy Synonym Search**  
  Checks if any registered synonym (alternate name, IUPAC, trade name) matches in PubChem’s database.
- **Tier 3: OPSIN (Open Parser for Systematic IUPAC Nomenclature)**  
  For unresolved names, attempts systematic IUPAC parsing via OPSIN, a specialized open-source tool for structure-from-name conversion.

### 3. Metadata and Traceability
- For each chemical, the following are recorded:
  - **Original Name**
  - **Assigned SMILES**
  - **SMILES_Status** (describes the method, variant, and number of attempts)
  - **SMILES_NameUsed** (the exact name or synonym that yielded the match)
- If all attempts fail, the status is labeled as “not found with any method.”

---

## Requirements

- Python 3.7+
- [`requests`](https://pypi.org/project/requests/)
- [`pubchempy`](https://pypi.org/project/PubChemPy/)
- [`tqdm`](https://pypi.org/project/tqdm/)
- [`pandas`](https://pypi.org/project/pandas/)
