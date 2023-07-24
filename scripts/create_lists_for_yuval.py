# %%
import pandas as pd 
import subprocess 
import json 
from pathlib import Path
import logging
import os

# %% 
# read in paths from json
with open('paths.json') as paths_file:
    file_contents = paths_file.read()

paths_json = json.loads(file_contents)
paths_dict = paths_json["turbo_id_analysis_pipeline_paths"]

# %%
# paths % configs 
commit_sha = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip() 
input_folder_path = Path(paths_dict['input_folder_path'])
output_folder_path = Path(paths_dict['output_folder_path'])
fasta_table_path = Path(paths_dict['fasta_table_path'])
logging.basicConfig(filename=output_folder_path / 'lists_for_yuval.log', 
                    filemode='w', 
                    encoding='utf-8', 
                    level=logging.INFO, 
                    format='%(asctime)s %(message)s', 
                    datefmt='%m/%d/%Y %I:%M:%S %p')
logging.info('Program starting')
logging.info('Gitcommit sha: %s', commit_sha)

# %% 
# get contaminant info 
fasta_table = pd.read_csv(fasta_table_path)
fasta_table = fasta_table.set_index("uniprot_id")
fasta_file_data = fasta_table[['description', 
                              'entry_name', 
                              'DB Class Key_jackson_homology_db',
                              'SWISS_PROT IDs_mouse_jackson_homology_db',
                              'Symbol_mouse_jackson_homology_db',
                              'SWISS_PROT IDs_human_jackson_homology_db',
                              'Symbol_human_jackson_homology_db',
                              'contaminant']]

# %%
# collect all final lists 
# 05_results / 01_tissue and 02_serum / for each file folder get final protein table 

dict_of_df_paths = {}

tissue_path = output_folder_path / "05_results" / "01_tissue"
serum_path = output_folder_path / "05_results" / "02_serum"

for paths in [tissue_path, serum_path]:
    for file_folder in list(paths.iterdir()):
        if file_folder.stem != '.DS_Store':
            file_name = file_folder.stem
            final_protein_table = [x for x in list(file_folder.iterdir()) if 'final_protein_table' in x.stem]
            if len(final_protein_table) == 1: 
                dict_of_df_paths[file_name] = final_protein_table[0]
            elif len(final_protein_table) == 2: 
                dict_of_df_paths[file_name] = final_protein_table[1]

# %%
list_of_values = ["Significant Down", "Significant Up"]
filtered_dfs_list = []
for file, path in dict_of_df_paths.items(): 
    logging.info(file)
    df = pd.read_csv(path)
    cols = df.filter(regex="Regulation").columns.tolist()
    df = df.set_index("uniprot_id")
    df = df[cols]
    df = df[df.isin(list_of_values).any(axis=1)]
    filtered_dfs_list.append(df)

# %%
assert len(filtered_dfs_list) == len(dict_of_df_paths)
# %%
final_list = pd.concat(filtered_dfs_list, axis=1)
cols = final_list.columns.tolist()
final_list = final_list.reset_index()

# %%
final_list = final_list.reset_index().merge(fasta_file_data.reset_index(), how="left", on="uniprot_id")
final_list = final_list.drop("index", axis= 1)

# %% 
folder_path = output_folder_path / "06_lists_for_gof_lof_variants_and_phewas"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)
# %%
final_list.to_csv(folder_path / 'list_for_gof_lof_variants_and_phewas.csv')
# %%
with_contaminans_path = "/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/230608_results_without_signalp/min_1_pep_per_protein/06_lists_for_gof_lof_variants_and_phewas/list_for_gof_lof_variants_and_phewas.csv"
without_contaminants_path = "/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/230720_results_without_signalp_crampome_removed/min_1_pep_per_protein/06_lists_for_gof_lof_variants_and_phewas/list_for_gof_lof_variants_and_phewas.csv"

with_df = pd.read_csv(with_contaminans_path)
with_df = with_df.set_index("uniprot_id")
with_df["with_contaminants"] = 1


without_df = pd.read_csv(without_contaminants_path)
without_df = without_df.set_index("uniprot_id")
without_df["without_contaminants"] = 1

# %%
comparison_df = pd.concat([with_df["with_contaminants"], without_df["without_contaminants"], fasta_table[["entry_name", "contaminant"]]], axis=1).to_csv("/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/comparison_with_without_contaminants.csv")
# %%