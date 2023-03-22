# TurboID

### Step 1 (mouse_fasta_to_FP_TP_tables.py)
- create main table (based off of our mouse fasta file used for IP2 search)
  - read in fasta file to get all uniprot IDs
  - use uniprot.org to get all annotations for these proteins
  - label them as TP or FP based on following criteria in subcellular location column: 
    - **TP**: "Endoplasmic Reticulum", "Secreted", "secreted", "Endoplasmic reticulum", "Rough endoplasmic reticulum", "endoplasmic reticulum"
    - **FP**: "Mitochondrion", "mitochondrion"
    - additionally make sure that if there is a SignalP annotation (ie if SignalP column is not empty) for a protein marked as FP then remove the FP annotation 
  - crossreference main list with spleen, adipose tissue etc. (files provided by Ken)
  - save table

### Step 2 ()