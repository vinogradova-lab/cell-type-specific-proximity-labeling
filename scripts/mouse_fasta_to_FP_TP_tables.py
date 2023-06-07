# %%

import pandas as pd
import numpy as np
from pathlib import Path
import openpyxl
import json
from functools import reduce 
from operator import concat
from src.fasta_table_funcs import *
import math
import warnings 
warnings.simplefilter("ignore")
%load_ext autoreload
%autoreload 2

# %% 
# read in json file as dict
with open('paths.json') as paths_file:
    file_contents = paths_file.read()

paths_json = json.loads(file_contents)
paths_dict = paths_json["mouse_fasta_to_FP_TP_tables_paths"]
with_signalp = True

# %%
# read in all paths 
mouse_fasta_path = Path(paths_dict['mouse_fasta_path']) 
tissue_localisation_path = Path(paths_dict['tissue_localisation_path']) 
uniprot_id_mapping_path = Path(paths_dict['uniprot_id_mapping_path'])
mitomatrix_protein_path = Path(paths_dict['mitomatrix_protein_path']) 
mart_export_path = Path(paths_dict['mart_export_path']) 
output_folder_path = Path(paths_dict['output_folder_path']) 
secretion_prediction_resource_path = Path(paths_dict['secretion_prediction_resource_path'])
human_mouse_homology_mapping_path = Path(paths_dict['human_mouse_homology_mapping_path'])
gprofiler_mouse_geneid_mapping_path = Path(paths_dict['gprofiler_mouse_geneid_mapping_path'])
leftover_uniprot_geneid_mapping_path = Path(paths_dict['leftover_uniprot_geneid_mapping_path'])

# %%
# run uniprot ID mapping with this file - (uniprot.org ID mapping)
main_fasta_table = create_df_from_fasta(mouse_fasta_path)
print(f"Size of main fasta table: {main_fasta_table.shape}.")
main_fasta_table.to_csv(output_folder_path / "01_before_uniprot_search" / "from_fasta_to_mouse_main_fasta_table_before_uniprot_search.csv") 
initial_number_of_entires = len(main_fasta_table)

# import uniprot id mapping results, merge with fasta table and remove duplicates 
uniprot_mapping = import_uniprot_annot(uniprot_id_mapping_path)
# print(f"Size of uniprot mapping table: {uniprot_mapping.shape}.")

# merge with main fasta table
merged_main_fasta_table = main_fasta_table.join(uniprot_mapping)
# print(f"Size of main merged table: {merged_main_fasta_table.shape}.")

# aggregate duplicates
# main_fasta_table has 21651 rows 
# uniprot_mapping has 21407 rows - because not everything could be found 
# when merged we see 21659 rows - because there are duplicates (6 uniprot ids, some of them triplicates)
# so we remove duplicates, aggregate them into a single row and put them back into our main table
merged_main_fasta_table = clean_up_fasta_table_merge(merged_main_fasta_table)
# print(f"Size of cleaned merged fasta table: {merged_main_fasta_table.shape}.")

# check if initial number of proteins is the same
assert initial_number_of_entires == len(merged_main_fasta_table)
print("Fasta table imported and annotated!")

# %% 
# crossreference table with secretion prediction resource shared by Corey 
merged_main_fasta_table = add_secretion_prediction_to_fasta_table(merged_main_fasta_table, secretion_prediction_resource_path)

# read in mitomatrix gene names - mitomatrix table is human, we work with mouse so we need to get corresponding mapping 
mouse_mitomatrix_genes_list = get_human_mitomatrix_mouse_ortholog(mitomatrix_protein_path, mart_export_path)

# annotate TP and FP
merged_main_fasta_table = annotate_TP_FP(with_signalp, merged_main_fasta_table, mouse_mitomatrix_genes_list)

if with_signalp: 
    assert len(merged_main_fasta_table.loc[merged_main_fasta_table["annotation"] == "TP"]) == 1602
    assert len(merged_main_fasta_table.loc[merged_main_fasta_table["annotation"] == "FP"]) == 438
else: 
    assert len(merged_main_fasta_table.loc[merged_main_fasta_table["annotation"] == "TP"]) == 2806
    assert len(merged_main_fasta_table.loc[merged_main_fasta_table["annotation"] == "FP"]) == 435
    # double check that there are no FPs with Signal peptide annotation
    assert merged_main_fasta_table[merged_main_fasta_table["annotation"]=="FP"]["Signal peptide"].notna().sum() == 0 
    # number of TP with Signal peptide
    assert merged_main_fasta_table[merged_main_fasta_table["annotation"]=="TP"]["Signal peptide"].value_counts().sum() == 1735

print("FP and TP annotated")

# %%
# crossreference main list with spleen, adipose tissue etc.
merged_main_fasta_table = add_tissue_annotation_from_ken(tissue_localisation_path, merged_main_fasta_table)
assert initial_number_of_entires == len(merged_main_fasta_table)
print("Added spleen, adipose tissue etc (files provided by Ken)")

# %% 
# read in gprofiler uniprot id to gene id mapping 
merged_main_fasta_table = add_gene_id_mapping(gprofiler_mouse_geneid_mapping_path, merged_main_fasta_table)
# add gene name column
merged_main_fasta_table['gene_name'] = merged_main_fasta_table["Gene Names"].str.split(" ").str[0] 
assert merged_main_fasta_table['gene_id'].count() == 19520
assert merged_main_fasta_table['gene_name'].count() == 21250
assert initial_number_of_entires == len(merged_main_fasta_table)
print("Gene id and gene name column added!")

# %% 
# add additional gene_ids from uniprot mapping of not mapped uniprots with gprofiler 
#gene_id_na = merged_main_fasta_table.loc[merged_main_fasta_table.gene_id.isna()]
leftover_mapping_df = pd.read_csv(leftover_uniprot_geneid_mapping_path)
leftover_mapping_df.columns = ["uniprot_id", "gene_id_leftover"]
leftover_mapping_df = leftover_mapping_df.groupby("uniprot_id").agg(list).reset_index()
leftover_mapping_df["gene_id_leftover"] = leftover_mapping_df['gene_id_leftover'].str[0]
leftover_mapping_df = leftover_mapping_df.set_index("uniprot_id")
merged_main_fasta_table = merged_main_fasta_table.join(leftover_mapping_df)
merged_main_fasta_table['gene_id'] = merged_main_fasta_table['gene_id'].fillna(merged_main_fasta_table['gene_id_leftover']) 
merged_main_fasta_table = merged_main_fasta_table.drop("gene_id_leftover", axis=1)
assert merged_main_fasta_table['gene_id'].count() == 20203

# %% 
# add human mouse homology uniprot ids 
gene_name_mapping_dict, uniprot_mapping_dict, human_mouse_homology_mapping_rearranged_df = add_human_uniprot_mapping(human_mouse_homology_mapping_path)
merged_main_fasta_table  = merged_main_fasta_table.reset_index()
merged_main_fasta_table["DB Class Key_jackson_homology_db"] = merged_main_fasta_table.apply(add_db_key, args=(uniprot_mapping_dict, gene_name_mapping_dict), axis=1)
merged_main_fasta_table = merged_main_fasta_table.merge(human_mouse_homology_mapping_rearranged_df, on="DB Class Key_jackson_homology_db", how="left")
assert initial_number_of_entires == len(merged_main_fasta_table)
print("Human uniprot ids from jacksons db added!")

# %% 
assert merged_main_fasta_table["DB Class Key_jackson_homology_db"].count() == 19272

# %% 
# add additional gene id mapping 
human_mouse_homology_mapping_df = pd.read_csv(human_mouse_homology_mapping_path, delimiter="\t")
mouse_df = human_mouse_homology_mapping_df.loc[human_mouse_homology_mapping_df["Common Organism Name"] == 'mouse, laboratory'].set_index("DB Class Key").add_suffix("_mouse")
mouse_gene_id_df = mouse_df["EntrezGene ID_mouse"].groupby("DB Class Key").agg(list).str[0].reset_index()
mouse_gene_id_df.columns = ["DB Class Key_jackson_homology_db", "EntrezGene ID_mouse"]
merged_main_fasta_table = merged_main_fasta_table.merge(mouse_gene_id_df, on="DB Class Key_jackson_homology_db", how="left")
merged_main_fasta_table['gene_id'] = merged_main_fasta_table['gene_id'].fillna(merged_main_fasta_table['EntrezGene ID_mouse']).dropna()
merged_main_fasta_table = merged_main_fasta_table.drop("EntrezGene ID_mouse", axis=1)
assert initial_number_of_entires == len(merged_main_fasta_table)
assert merged_main_fasta_table['gene_id'].count() == 20723
# %%
# save final table 
if with_signalp:
    merged_main_fasta_table.to_csv(output_folder_path / "03_table_for_analysis" / "main_fasta_table_with_signal_p.csv")
else:
    merged_main_fasta_table.to_csv(output_folder_path / "03_table_for_analysis" / "main_fasta_table_without_signal_p.csv")

print("File saved, with signalp = ", with_signalp)
# %%
