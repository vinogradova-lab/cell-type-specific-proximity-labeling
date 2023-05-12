# %%
from Bio import SeqIO
import pandas as pd
import numpy as np
from pathlib import Path
import openpyxl
import json
from functools import reduce 
from operator import concat
from src.fasta_table_funcs import *
%load_ext autoreload
%autoreload 2

# %% 
with open('paths.json') as paths_file:
    file_contents = paths_file.read()

paths_json = json.loads(file_contents)
paths_dict = paths_json["mouse_fasta_to_FP_TP_tables_paths"]

# %%
# read in all paths 
mouse_fasta_path = Path(paths_dict['mouse_fasta_path']) 
tissue_localisation_path = Path(paths_dict['tissue_localisation_path']) 
mitomatrix_protein_path = Path(paths_dict['mitomatrix_protein_path']) 
#mitomatrix_mouse_ortholog_path = Path(paths_dict['mitomatrix_mouse_ortholog_path']) 
mart_export_path = Path(paths_dict['mart_export_path']) 
output_folder_path = Path(paths_dict['output_folder_path']) 
secretion_prediction_resource_path = Path(paths_dict['secretion_prediction_resource_path'])

# %% 
spr_df = pd.read_csv(secretion_prediction_resource_path)
spr_df = spr_df.rename(columns={"Accession":"uniprot_sp_key", "Organism":"organism"})

# %%
# read in fasta file as dict 
fasta_dict = {record.id : record.description for record in SeqIO.parse(mouse_fasta_path, "fasta")}

# %%
# create main fasta table from fasta entries, remove reverse
main_fasta_table = pd.DataFrame.from_dict(fasta_dict, orient='index')
main_fasta_table = main_fasta_table.reset_index()
main_fasta_table.columns = ["uniprot_id", "description"]
main_fasta_table["id"] = main_fasta_table["uniprot_id"]
main_fasta_table[["db","uniprot_id","entry_name"]] = main_fasta_table["uniprot_id"].str.split("|", expand=True)
# filter out reverse entires
main_fasta_table = main_fasta_table[~main_fasta_table["id"].str.contains("Reverse_")]
main_fasta_table = main_fasta_table.set_index("uniprot_id")
print(main_fasta_table.shape)
main_fasta_table.head()

# %%
# run uniprot ID mapping with this file - (uniprot.org ID mapping)
main_fasta_table.to_csv(output_folder_path / "01_before_uniprot_search" / "from_fasta_to_mouse_main_fasta_table_before_uniprot_search.csv") 
initial_number_of_entires = len(main_fasta_table)

# %%
# import uniprot id mapping results 
uniprot_id_mapping_path = Path(paths_dict['uniprot_id_mapping_path'])

# %%
# read in uniprot id mapping results, annotate keratins
uniprot_mapping = pd.read_csv(uniprot_id_mapping_path)
uniprot_mapping["keratin"] = uniprot_mapping["Protein names"].str.contains("keratin", case=False)
uniprot_mapping.rename({"From":"uniprot_id"}, axis=1, inplace=True)
uniprot_mapping = uniprot_mapping.set_index("uniprot_id")
print(uniprot_mapping.shape)
uniprot_mapping.head()

# %%
# merge with main fasta table
merged_main_fasta_table = main_fasta_table.join(uniprot_mapping)
# change dtypes so that we can concat later 
merged_main_fasta_table = merged_main_fasta_table.astype({'Length': 'object', 'Mass': 'object'})
merged_main_fasta_table.shape

# %%
# aggregate duplicates
# main_fasta_table has 21651 rows 
# uniprot_mapping has 21407 rows - because not everything could be found 
# when merged we see 21659 rows - because there are duplicates (6 uniprot ids, some of them triplicates)
# so we remove duplicates, aggregate them into a single row and put them back into our main table
duplicates_list = merged_main_fasta_table[merged_main_fasta_table.index.duplicated()].index.tolist()
duplicates_df = merged_main_fasta_table[merged_main_fasta_table.index.isin(duplicates_list)]
deduplicated_entries = duplicates_df.groupby("uniprot_id").agg(list) #.reset_index()
deduplicated_entries = deduplicated_entries.apply(remove_duplicates_from_aggregate, axis=1)
deduplicated_entries["Gene Names"] = deduplicated_entries["Gene Names"].str.replace(",", "")

merged_main_fasta_table = merged_main_fasta_table.drop(index=duplicates_list)
merged_main_fasta_table = pd.concat([merged_main_fasta_table, deduplicated_entries])
merged_main_fasta_table.shape

# %% 
# check if initial number of proteins is the same
assert initial_number_of_entires == len(merged_main_fasta_table)


#crossreference with secretion prediction resource shared by Corey 
sp_list = spr_df.reset_index().uniprot_sp_key.tolist()
merged_main_fasta_table = merged_main_fasta_table.reset_index()
merged_main_fasta_table['Entry'] = merged_main_fasta_table['Entry'].astype(str)
merged_main_fasta_table['uniprot_id'] = merged_main_fasta_table['uniprot_id'].astype(str)
merged_main_fasta_table["uniprot_sp_key"] = merged_main_fasta_table.apply(check_uniprot, args=(sp_list,), axis=1)
merged_main_fasta_table = merged_main_fasta_table.merge(spr_df, on="uniprot_sp_key", how="left")
merged_main_fasta_table = merged_main_fasta_table.set_index('uniprot_id')

# %% 
# read in mitomatrix gene names - mitomatrix table is human, we work with mouse so we crossreference by gene name  
mitomatrix_df = pd.read_csv(mitomatrix_protein_path)
mitomatrix_df['Gene Names'] = mitomatrix_df['Gene Names'].apply(lambda x: x.split(";"))
mitomatrix_df = mitomatrix_df.explode("Gene Names")
mitomatrix_df = mitomatrix_df.rename(columns={"Gene Names":"Gene name"})
mitomatrix_df.head()

# we need to be careful when crossreferencing by gene names
# https://www.biostars.org/p/149115/
# https://pypi.org/project/mousipy/
# https://bioinformatics.stackexchange.com/questions/17486/converting-mouse-genes-to-human-genes


# %%
# follow https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/ turotial and read in results
mart_human_mouse_orthologs = pd.read_csv(mart_export_path)

# %% 
# merge mart results to get mouse gene list 
mitomatrix_mart_merge_df = mitomatrix_df.merge(mart_human_mouse_orthologs, how="left", on="Gene name")
mouse_mitomatrix_genes_list = mitomatrix_mart_merge_df["Mouse gene name"].tolist()

# %%
# annotate TP and FP
# merged_main_fasta_table = merged_main_fasta_table.reset_index()
merged_main_fasta_table["TP"] = merged_main_fasta_table.apply(annotate_TP_signalp, axis=1) 

# %% 
# merged_main_fasta_table["FP"] = merged_main_fasta_table["Subcellular location [CC]"].apply(annotate_FP) 
list_of_found_genes = []
merged_main_fasta_table["FP"] = merged_main_fasta_table["Gene Names"].apply(annotate_FP_mitomatrix, args=(mouse_mitomatrix_genes_list, )) 
print(merged_main_fasta_table["TP"].value_counts())
print(merged_main_fasta_table["FP"].value_counts())

# %%
# clean up annotation 
merged_main_fasta_table["annotation"] = merged_main_fasta_table.apply(conclude_annotation_signalp, axis=1) 
merged_main_fasta_table["annotation"].value_counts()

# %%
# this is not necessary 
# double check that there are no FPs with Signal peptide annotation
#assert merged_main_fasta_table[merged_main_fasta_table["annotation"]=="FP"]["Signal peptide"].notna().sum() == 0 
# number of TP with Signal peptide
#merged_main_fasta_table[merged_main_fasta_table["annotation"]=="TP"]["Signal peptide"].value_counts().sum()

# %%
# check spleen, adipose tissue, etc file names
files_in_folder = list(tissue_localisation_path.iterdir())
[file_path.stem for file_path in files_in_folder]

# %%
# crossreference main list with spleen, adipose tissue etc.
for file_path in files_in_folder:
    column_name = file_path.stem.split("_", 1)[0]
    df = pd.read_excel(file_path, engine='openpyxl')
    entry_list = df.Entry.tolist()
    merged_main_fasta_table[column_name] = merged_main_fasta_table.index.isin(entry_list)

# %%
# check size of final table
assert initial_number_of_entires == len(merged_main_fasta_table)
print(merged_main_fasta_table.shape)
merged_main_fasta_table.head()

# %%
# save final table 
merged_main_fasta_table.to_csv(output_folder_path / "03_table_for_analysis" / "main_fasta_table_with_signal_p.csv")