# TurboID


### **Step 1 (processing of raw census-out files)**

Files were processed with the Vinogradova Lab pipeline and were treated as “whole proteome” files – i.e., we quantify protein expression changes and corresponding peptide values were summed to receive the respective raw signal intensity per protein

--> Processing details
- Files were processed twice: 
  - with the restriction that at least 1 peptide must be present per protein
  - with the restriction that at least 2 peptides must be present per protein 
- Additionally following filters were used:
  - Non unique peptides were removed
  - Peptides that are not fully tryptic were removed
  - Peptides with more than 1 missed cleavage site were removed


### **Step 2 (mouse_fasta_to_FP_TP_tables.py)**
- create main table (based off of our mouse fasta file used for IP2 search)
  - read in fasta file to get all uniprot IDs
  - use uniprot.org to get all annotations for these proteins
  - label them as TP or FP based on following criteria in subcellular location column: 
    - **TP**: "Endoplasmic Reticulum", "Secreted", "secreted", "Endoplasmic reticulum", "Rough endoplasmic reticulum", "endoplasmic reticulum"
    - **FP**: "Mitochondrion", "mitochondrion"
    - additionally make sure that if there is a SignalP annotation (ie if SignalP column is not empty) for a protein marked as FP then remove the FP annotation 
  - crossreference main list with spleen, adipose tissue etc. (files provided by Ken)
  - save table

### Step 3 ()