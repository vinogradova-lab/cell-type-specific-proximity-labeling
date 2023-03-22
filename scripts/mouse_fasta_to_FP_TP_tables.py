# %%
from Bio import SeqIO
import pandas as pd
import numpy as np
from pathlib import Path
import openpyxl

#%%
#import functions
def remove_duplicates_from_aggregate(row): 
    list_of_cols = row.index.tolist()
    for col in list_of_cols: 
        row[col] = ', '.join(map(str, set(row[col])))
    return row

def annotate_TP(subcelullar_location):
    tp_term_list = ["Endoplasmic Reticulum", "Secreted", "secreted", "Endoplasmic reticulum", "Rough endoplasmic reticulum", "endoplasmic reticulum"]
    #"Golgi", "Extracellular space", "extracellular space", 
    #
    subcelullar_location = str(subcelullar_location)
    if "Note:" in subcelullar_location:
        subcelullar_location = subcelullar_location.split("Note:")[0]
    for item in tp_term_list: 
        if item in subcelullar_location:
            return("TP")
    
def annotate_FP(subcelullar_location):
    fp_term_list = ["Mitochondrion", "mitochondrion"] #'cytoskeleton' "Cytoskeleton", "Nucleus", "nucleus", 
    subcelullar_location = str(subcelullar_location)
    if "Note:" in subcelullar_location:
        subcelullar_location = subcelullar_location.split("Note:")[0]
    for item in fp_term_list: 
        if item in subcelullar_location:
            return("FP")

def annotation(df):   
    if df["TP"] == "TP":
        return("TP")
    #if signal peptide is present then even though annotation says FP, change annotation to nan
    elif df["FP"] == "FP" and pd.notna(df["Signal peptide"]):
        return(np.nan)
    elif df["FP"] == "FP" and pd.isnull(df["Signal peptide"]):
        return("FP")

# %%
#read in all paths 
mouse_fasta_path = Path("../data/fasta_files/radus_database__cravattlab_mouse_nonredundant_07-11-2017_reversed.fasta")
tissue_localisation_path = Path("/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Reference datasets/Mouse/Tissue localisation")


# %%
#read in fasta file as dict 
fasta_dict = {record.id : record.description for record in SeqIO.parse(mouse_fasta_path, "fasta")}

# %%
#create main fasta table from fasta entries, remove reverse
main_fasta_table = pd.DataFrame.from_dict(fasta_dict, orient='index')
main_fasta_table = main_fasta_table.reset_index()
main_fasta_table.columns = ["uniprot_id", "description"]
main_fasta_table["id"] = main_fasta_table["uniprot_id"]
main_fasta_table[["db","uniprot_id","entry_name"]] = main_fasta_table["uniprot_id"].str.split("|", expand=True)
#filter out reverse entires
main_fasta_table = main_fasta_table[~main_fasta_table["id"].str.contains("Reverse_")]
main_fasta_table = main_fasta_table.set_index("uniprot_id")
print(main_fasta_table.shape)
main_fasta_table.head()

# %%
#run uniprot ID mapping with this file - (uniprot.org ID mapping)
main_fasta_table.to_csv("../data/main_fasta_table/01_before_uniprot_search/from_fasta_to_mouse_main_fasta_table_before_uniprot_search.csv") 
initial_number_of_entires = len(main_fasta_table)

# %%
#import uniprot id mapping results 
uniprot_id_mapping_path = Path("../data/main_fasta_table/02_after_uniprot_mapping/mouse_mapping_uniprot-download_true_fields_accession_2Creviewed_2Cid_2Cprotein_nam-2022.08.31-15.17.09.71.csv")

# %%
#read in uniprot id mapping results, annotate keratins
uniprot_mapping = pd.read_csv(uniprot_id_mapping_path)
uniprot_mapping["keratin"] = uniprot_mapping["Protein names"].str.contains("keratin", case=False)
uniprot_mapping.rename({"From":"uniprot_id"}, axis=1, inplace=True)
uniprot_mapping = uniprot_mapping.set_index("uniprot_id")
print(uniprot_mapping.shape)
uniprot_mapping.head()

# %%
#merge with main fasta table
merged_main_fasta_table = main_fasta_table.join(uniprot_mapping)
#change dtypes so that we can concat later 
merged_main_fasta_table = merged_main_fasta_table.astype({'Length': 'object', 'Mass': 'object'})
merged_main_fasta_table.shape

# %%
#aggregate duplicates
#main_fasta_table has 21651 rows 
#uniprot_mapping has 21407 rows - because not everything could be found 
#when merged we see 21659 rows - because there are durplicates (6 uniprot ids, some of them triplicates)
#so we remove duplicates, aggregate them into a single row and put them back into our main table
duplicates_list = merged_main_fasta_table[merged_main_fasta_table.index.duplicated()].index.tolist()
duplicates_df = merged_main_fasta_table[merged_main_fasta_table.index.isin(duplicates_list)]
deduplicated_entries = duplicates_df.groupby("uniprot_id").agg(list) #.reset_index()
deduplicated_entries = deduplicated_entries.apply(remove_duplicates_from_aggregate, axis=1)

merged_main_fasta_table = merged_main_fasta_table.drop(index=duplicates_list)
merged_main_fasta_table = pd.concat([merged_main_fasta_table, deduplicated_entries])
merged_main_fasta_table.shape

#%% 
#check if initial number of proteins is the same
assert initial_number_of_entires == len(merged_main_fasta_table)

# %%
#annotate TP and FP
merged_main_fasta_table["TP"] = merged_main_fasta_table["Subcellular location [CC]"].apply(annotate_TP) 
merged_main_fasta_table["FP"] = merged_main_fasta_table["Subcellular location [CC]"].apply(annotate_FP) 
print(merged_main_fasta_table["TP"].value_counts())
print(merged_main_fasta_table["FP"].value_counts())

# %%
#clean up annotation by making sure that if there is signal peptide then FP annotation is removed 
merged_main_fasta_table["annotation"] = merged_main_fasta_table.apply(annotation, axis=1) 
merged_main_fasta_table["annotation"].value_counts()

# %%
#double check that there are no FPs with Signal peptide annotation
assert merged_main_fasta_table[merged_main_fasta_table["annotation"]=="FP"]["Signal peptide"].notna().sum() == 0 

# %%
#number of TP with Signal peptide
merged_main_fasta_table[merged_main_fasta_table["annotation"]=="TP"]["Signal peptide"].value_counts().sum()

# %%
#check spleen, adipose tissue, etc file names
files_in_folder = list(tissue_localisation_path.iterdir())
[file_path.stem for file_path in files_in_folder]

# %%
#crossreference main list with spleen, adipose tissue etc.
for file_path in files_in_folder:
    column_name = file_path.stem.split("_", 1)[0]
    df = pd.read_excel(file_path, engine='openpyxl')
    entry_list = df.Entry.tolist()
    merged_main_fasta_table[column_name] = merged_main_fasta_table.index.isin(entry_list)

# %%
#check size of final table
assert initial_number_of_entires == len(merged_main_fasta_table)
print(merged_main_fasta_table.shape)
merged_main_fasta_table.head()

# %%
#save final table 
merged_main_fasta_table.to_csv("../data/main_fasta_table/03_table_for_analysis/main_fasta_table.csv")
