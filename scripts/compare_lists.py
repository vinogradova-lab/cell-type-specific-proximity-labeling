# %% 
import pandas as pd
# %%
input_path_old = "/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/050423_results/min_1_pep_per_protein/05_results/01_tissue/07302022_EV1-100A_10pl_M/final_protein_table_07302022_EV1-100A_10pl_M.csv"
input_path_new = "/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/051023_test_results_dont_use/min_1_pep_per_protein/05_results/01_tissue/07302022_EV1-100A_10pl_M/final_protein_table_07302022_EV1-100A_10pl_M.csv"

#output_path - ""

# %% 

old_df = pd.read_csv(input_path_old)
old_df = old_df.set_index("uniprot_id")
old_df = old_df.add_prefix("OLD_LISTS_")
# %%
new_df = pd.read_csv(input_path_new)
new_df = new_df.set_index("uniprot_id")
new_df = new_df.add_prefix("NEW_LISTS_")
# %%
joined_df = old_df.join(new_df)
joined_df.to_csv("100A_old_vs_new_list.csv")

# %%
# %%
input_path_old = "/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/050423_results/min_1_pep_per_protein/05_results/01_tissue/08072022_EV101_16pl_M/final_protein_table_08072022_EV101_16pl_M.csv"
input_path_new = "/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/051023_test_results_dont_use/min_1_pep_per_protein/05_results/01_tissue/08072022_EV101_16pl_M/final_protein_table_08072022_EV101_16pl_M.csv"

# %% 

old_df = pd.read_csv(input_path_old)
old_df = old_df.set_index("uniprot_id")
old_df = old_df.add_prefix("OLD_LISTS_")
# %%
new_df = pd.read_csv(input_path_new)
new_df = new_df.set_index("uniprot_id")
new_df = new_df.add_prefix("NEW_LISTS_")
# %%
joined_df = old_df.join(new_df)
joined_df.to_csv("101_old_vs_new_list.csv")

# %%
