import pandas as pd
import numpy as np

def remove_duplicates_from_aggregate(row): 
    list_of_cols = row.index.tolist()
    for col in list_of_cols: 
        row[col] = ', '.join(map(str, set(row[col])))
    return row


def annotate_TP(subcelullar_location):
    tp_term_list = ["Endoplasmic Reticulum", "Secreted", "secreted", "Endoplasmic reticulum", "Rough endoplasmic reticulum", "endoplasmic reticulum"]
    # "Golgi", "Extracellular space", "extracellular space", 
    subcelullar_location = str(subcelullar_location)
    if "Note:" in subcelullar_location:
        subcelullar_location = subcelullar_location.split("Note:")[0]
    for item in tp_term_list: 
        if item in subcelullar_location:
            return "TP" 
   

def annotate_FP_subcellular_loc(subcelullar_location):
    fp_term_list = ["Mitochondrion", "mitochondrion"] # removed: 'cytoskeleton' "Cytoskeleton", "Nucleus", "nucleus", 
    subcelullar_location = str(subcelullar_location)
    if "Note:" in subcelullar_location:
        subcelullar_location = subcelullar_location.split("Note:")[0]
    for item in fp_term_list: 
        if item in subcelullar_location:
            return "FP" 

def annotate_FP_mitomatrix(gene_names, mitomatrix_list):
    gene_names = str(gene_names).upper()
    gene_name_list = gene_names.split(',')
    for gene in gene_name_list:
        if gene_names in mitomatrix_list:
            return "FP"

def conclude_annotation(df):   
    if df["TP"] == "TP":
        return "TP" 
    # if signal peptide is present then even though annotation says FP, change annotation to nan
    elif df["FP"] == "FP" and pd.notna(df["Signal peptide"]):
        return np.nan
    elif df["FP"] == "FP" and pd.isnull(df["Signal peptide"]):
        return "FP"