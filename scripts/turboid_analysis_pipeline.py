# %%
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from pathlib import Path
from sklearn.decomposition import PCA
import plotly.express as px
import plotly.graph_objects as go
import plotly.express as px
import json
from src.cutoff_funcs import *
from src.filter_funcs import *
from src.norm_funcs import *
from src.readin_funcs import *
import logging
import subprocess 
import warnings 
warnings.filterwarnings('ignore')
%load_ext autoreload
%autoreload 2

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
logging.basicConfig(filename=output_folder_path / 'turboid_analysis.log', filemode='w', encoding='utf-8', level=logging.INFO)
logging.info(commit_sha)

# %%
# read in fasta table containing TP and FP annotation 
fasta_table = pd.read_csv(fasta_table_path, index_col=[0])
TP_list = fasta_table[fasta_table["annotation"] == "TP"]
FP_list = fasta_table[fasta_table["annotation"] == "FP"]

logging.info("Number of Uniport IDs in TP list: %s", len(TP_list))
logging.info("Number of Uniport IDs in FP list: %s", len(FP_list))

assert len(TP_list) == 2806
assert len(FP_list) == 435

# %%
# Get list of files, file_channel_dict and cond_dict of processed census-out files
list_of_file_paths = list(input_folder_path.iterdir())
list_of_file_paths = [x for x in list_of_file_paths if 'census-out' in x.stem]

list_of_file_names = [file_path.stem for file_path in list_of_file_paths]
file_channel_dict = get_channel_name_dict(input_folder_path, list_of_file_names)
file_condition_dict = get_cond_info(input_folder_path, list_of_file_names)

logging.info("Number of files: %s", len(list_of_file_names))
logging.info(list_of_file_names)

assert len(file_channel_dict) == len(list_of_file_names) == len(file_condition_dict)

# %%
# Assign new column names accoriding to metadata_col.csv and remove keratins in each file and annotate TP and FP
keratins = fasta_table[fasta_table.keratin == True]
keratins_list = keratins.index.tolist()

folder_path = output_folder_path / "01_raw_files_with_correct_channel_names_keratins_removed_TP_FP_annotated"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

logging.info('\n')
logging.info("Assign new column names accoriding to metadata_col.csv and remove keratins in each file and annotate TP and FP")

dfs_dict = {}
for file in list_of_file_paths:
    file_name = file.stem
    logging.info("Processing: %s", file_name)
    df = pd.read_csv(file)
    
    # remove keratins 
    before_keratin_removal = len(df)
    df = df[~df.uniprot.isin(keratins_list)]
    after_keratin_removal = len(df)
    logging.info("Removed %s keratins", before_keratin_removal-after_keratin_removal)

    # rename columns
    df = df.set_index(["uniprot", 'description', 'pep_num'])
    df.rename(columns=file_channel_dict[file_name], inplace=True)
    
    # annotate TP FP
    df = df.reset_index()
    df = df.rename(columns={"uniprot": "Entry"})
    merged = df.merge(fasta_table[["annotation", "Entry"]], on="Entry", how="left")
    merged = merged.drop_duplicates()
    merged = merged.rename(columns={"Entry": "uniprot_id"})

    assert after_keratin_removal == len(merged)

    merged_df = merged.set_index(["uniprot_id", "description", "pep_num", "annotation"])
    dfs_dict[file_name] = merged_df
    merged_df.to_csv(folder_path / (file_name.split("processed_census-out_")[1] + ".csv"))

logging.info("FOLDER COMPLETE: 01_raw_files_with_correct_channel_names_keratins_removed_TP_FP_annotated")

# %%
# Filter based on condition 
folder_path = output_folder_path / "02_filtering_per_cond_per_file" 
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

logging.info('\n')
logging.info("Filter based on cre+ channels per condition per file")

filtered_dict = {}
for file_name in list_of_file_names:
    file_folder_path = folder_path / file_name.split("processed_census-out_")[1]
    if not os.path.exists(file_folder_path):
        os.mkdir(file_folder_path)

    logging.info("Processing: %s", file_name.split("processed_census-out_")[1])
    
    df = dfs_dict[file_name]
    before_filtering = len(df)
    conditions_list = file_condition_dict[file_name]["conditions"]
    control_labelling = file_condition_dict[file_name]["control_labelling"]
    treatment_labelling = file_condition_dict[file_name]["treatment_labelling"]

    #filter each condition in file 
    filtered_sub_dfs = []
    for condition in conditions_list:
        sub_df = get_condition_df(df, condition)
        filtered_sub_df, cond_columns, ctrl_columns = filter_condition_df(sub_df, treatment_labelling, control_labelling)
        filtered_sub_df = filtered_sub_df.reset_index()
        filtered_sub_dfs.append(filtered_sub_df)
    
    merge_on_cols = ['uniprot_id','description', 'pep_num', "annotation"]
    filtered_final_table = reduce(lambda df1,df2: pd.merge(df1,df2,on=merge_on_cols, how="outer"), filtered_sub_dfs)

    raw_df = df.reset_index()
    uniprot_ids = filtered_final_table["uniprot_id"].tolist()

    sub_raw_df = raw_df[raw_df['uniprot_id'].isin(uniprot_ids)]
    sub_raw_df = sub_raw_df.set_index(["uniprot_id", 'description', 'pep_num', 'annotation'])

    logging.info("Proteins were filtered out, because they did not pass at least one condition: %s", before_filtering-len(sub_raw_df))
    logging.info("Filtering for %s done", file_name)

    filtered_dict[file_name] = sub_raw_df
    sub_raw_df.to_csv(file_folder_path / ("passed_filter_" + file_name + ".csv"))

    # save filtered out proteins
    filtered_out_df = raw_df[~raw_df['uniprot_id'].isin(uniprot_ids)]
    filtered_out_df = filtered_out_df.set_index(["uniprot_id", 'description', 'pep_num', 'annotation'])
    filtered_out_df.to_csv(file_folder_path / ("removed_" + file_name + ".csv"))


# %%
# (1) generate the output put with raw SI and annotation TP/FP; 
# (2) look at TP proteins with 2+ peptides; 
# (3) calculate median SI for each cre+ channel and median of sums; 
# (4) calculate normalization ratios by dividing (median of median)/(median cre+ channel), and 
# (5) use those normalization values for each channel;
folder_path = output_folder_path / "03_normalisation_plots"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

logging.info('\n')
logging.info("Creating normalized tables and normalisation plots")

normalized_dict = {}
for file_name in list_of_file_names:
    file_folder_path = folder_path / file_name.split("processed_census-out_")[1]
    if not os.path.exists(file_folder_path):
        os.mkdir(file_folder_path)

    logging.info("Processing %s", file_name)
    
    plot_list = []
    fig = plt.figure(figsize=(25, 13)) 
    fig.suptitle(file_name)
    ax1 = plt.subplot(121) 
    ax2 = plt.subplot(122) 
    
    df = filtered_dict[file_name]
    conditions_list = file_condition_dict[file_name]["conditions"]
    control_labelling = file_condition_dict[file_name]["control_labelling"]
    treatment_labelling = file_condition_dict[file_name]["treatment_labelling"]
    
    # RAW
    pca_fig_1 = get_pca_plot(df, "Raw")
    plot_list.append(pca_fig_1)
    
    df_log = np.log2(df)
    df_log.plot(kind='box', 
                rot=90, 
                fontsize=10, 
                title="Raw", 
                ax=ax1, 
                xlabel="Channels", 
                ylabel="log2(SI)", 
                color = "green")
    
    
    # MEDIAN OF MEDIANS TP only, 2+ peptides
    norm_factor_series = []
    for condition in conditions_list:
        sub_df = get_condition_df(df, condition)
        cond_columns = sub_df.filter(like=treatment_labelling).columns.tolist()
        ctrl_columns = sub_df.filter(like=control_labelling).columns.tolist()
        sub_df["median_"+condition+"_"+treatment_labelling] = sub_df[cond_columns].median(axis=1)
        sub_df["median_"+condition+"_"+control_labelling] = sub_df[ctrl_columns].median(axis=1)
        sub_df["median_ratio"] = sub_df["median_"+condition+"_"+treatment_labelling] / sub_df["median_"+condition+"_"+control_labelling]
        sub_df = sub_df.reset_index()
        sub_df = sub_df[(sub_df["annotation"] == "TP") & (sub_df["pep_num"] >= 2)]
        sub_df = sub_df.sort_values(by=['median_ratio'], ascending=False)
        sub_df = sub_df.drop("median_"+condition+"_"+treatment_labelling, axis=1)
        sub_df = sub_df.drop("median_"+condition+"_"+control_labelling, axis=1)
        sub_df = sub_df.drop("median_ratio", axis=1)
        sub_df = sub_df.set_index(["uniprot_id", "description", "pep_num", "annotation"])
        sub_df = sub_df[cond_columns]
        
        norm_factors = sub_df.median().median() / sub_df.median()
        norm_factor_series.append(norm_factors)

    norm_factor = pd.concat(norm_factor_series, axis=1).sum(1)
    logging.info("Normalization factors")
    logging.info(norm_factor)

    normalized_df = df[norm_factor.index] * norm_factor
    list_of_ctrl_cols = [column for column in df.columns if column not in norm_factor.index]
    untouched_ctrl_columns = df[list_of_ctrl_cols]
    
    norm_df = normalized_df.reset_index().merge(untouched_ctrl_columns.reset_index(), on=["uniprot_id", "description", "pep_num", "annotation"])
    norm_df = norm_df.set_index(["uniprot_id", "description", "pep_num", "annotation"])
    norm_df = norm_df[df.columns] # get the same order of columns 
    
    assert len(norm_df) == len(df)
    assert norm_df.columns.tolist() == df.columns.tolist()

    pca_fig_2 = get_pca_plot(norm_df, "Median of medians TP only, 2+ peptides, Cre+")
    plot_list.append(pca_fig_2)
    
    norm_df_log = np.log2(norm_df)
    norm_df_log.plot(kind='box', 
                     rot=90, 
                     fontsize=10, 
                     title="Median of medians TP only, 2+ peptides, Cre+", 
                     ax=ax2, 
                     xlabel="Channels", 
                     ylabel="log2(SI)", 
                     color = "darkblue")
    
    # this is what gets passed into cutoff analysis
    normalized_dict[file_name] = norm_df
    
    fig.savefig(file_folder_path / ("box_plots_" + file_name + '.png'))
    with open(file_folder_path / ("pca_plots_" + file_name + '.html'), 'w') as f:
        f.write(file_name)
        for fig in plot_list:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))

# %%
# save normalized df in folder
folder_path = output_folder_path / "04_normalized"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

for file_name in list_of_file_names:
    df = normalized_dict[file_name]
    df.to_csv(folder_path / (file_name.split("processed_census-out_")[1] + ".csv"))


# %% 
# split tissue and serum 
folder_path = output_folder_path / "05_results"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

tissue_folder_path = folder_path / "01_tissue"
if not os.path.exists(tissue_folder_path):
    os.mkdir(tissue_folder_path)

serum_folder_path = folder_path / "02_serum"
if not os.path.exists(serum_folder_path):
    os.mkdir(serum_folder_path)

for file_name in list_of_file_names:
    # get data
    df = normalized_dict[file_name]
    conditions_list = file_condition_dict[file_name]["conditions"]
    control_labelling = file_condition_dict[file_name]["control_labelling"]
    treatment_labelling = file_condition_dict[file_name]["treatment_labelling"]
    file_type = file_condition_dict[file_name]["file_type"]
    cutoff_dict = {}
    
    if file_type == "serum": 
        serum_file_folder_path = serum_folder_path / file_name.split("processed_census-out_")[1]
        if not os.path.exists(serum_file_folder_path):
            os.mkdir(serum_file_folder_path)
        df = add_enrichment_ratio_serum_samples(df, conditions_list, control_labelling, treatment_labelling)
        volcano_df = get_volcano_plot(conditions_list, control_labelling, treatment_labelling, df, file_name, serum_file_folder_path)
        annotated_protein_df = get_detailed_protein_annotation(df, volcano_df, fasta_table)
        annotated_protein_df.to_csv(serum_file_folder_path / ("final_protein_table_" + file_name.split("processed_census-out_")[1] + ".csv"))

    elif file_type == "tissue":
        tissue_file_folder_path = tissue_folder_path / file_name.split("processed_census-out_")[1]
        if not os.path.exists(tissue_file_folder_path):
           os.mkdir(tissue_file_folder_path)

        cutoff_roc_path = tissue_file_folder_path / "01_cutoff_ROC_filter"
        if not os.path.exists(cutoff_roc_path):
            os.mkdir(cutoff_roc_path)
        decision_table, cutoff_plots_table = get_ratios_and_cutoffs(df, conditions_list, control_labelling, treatment_labelling, cutoff_roc_path, cutoff_dict, file_name, fasta_table)
        get_tp_fp_cutoff_plots(cutoff_plots_table, cutoff_roc_path, file_name, cutoff_dict)

        pass_cutoff_true_df = get_before_after_cutoff_barplots(decision_table, tissue_file_folder_path, file_name)
        
        pass_cutoff_df_norm_data = df.reset_index()[df.reset_index()["uniprot_id"].isin(pass_cutoff_true_df.uniprot_id.tolist())].set_index(["uniprot_id", "description", "pep_num", "annotation"])
        pass_cutoff_true_df = pass_cutoff_true_df.drop("index", axis=1).set_index("uniprot_id")

        if file_name == "processed_census-out_04172023_CRW_A-5_16pl_M":
            volcano_df = get_volcano_plot_treatment_vs_control(conditions_list, control_labelling, treatment_labelling, pass_cutoff_df_norm_data, file_name, tissue_file_folder_path)
            pass_cutoff_true_df = pass_cutoff_true_df.join(volcano_df)
            pass_cutoff_true_df.to_csv(tissue_file_folder_path / ("BAT_final_protein_table" + file_name.split("processed_census-out")[1] +'.csv'))

            pass_cutoff_true_df = decision_table[decision_table[['GFAP_brain_Cre(+)_Ctrl_pass_cutoff_sum', 'GFAP_brain_Cre(+)_Fast_pass_cutoff_sum', 'GFAP_brain_Cre(+)_LPS_pass_cutoff_sum']].sum(axis=1) >= 1]
            pass_cutoff_cols = pass_cutoff_true_df.filter(like="pass_cutoff").columns.tolist()
            pass_cutoff_true_df = pass_cutoff_true_df.drop(pass_cutoff_cols, axis = 1)
            pass_cutoff_df_norm_data = df.reset_index()[df.reset_index()["uniprot_id"].isin(pass_cutoff_true_df.uniprot_id.tolist())].set_index(["uniprot_id", "description", "pep_num", "annotation"])

            #get scatterplots comaprisons: W7 vs W6, W8 vs W7, W9 vs W7, and W8 vs W9
            pass_cutoff_true_df = pass_cutoff_true_df.set_index("entry_name")
            df_scat_1 = scatterplot_plot(pass_cutoff_true_df[['ratio_GFAP_brain_Cre(+)_Ctrl_W7/GFAP_brain_Cre(-)_Ctrl_W6', 'ratio_GFAP_brain_Cre(+)_Fast_W8/GFAP_brain_Cre(-)_Ctrl_W6']], tissue_file_folder_path, "W7_vs_W8")
            df_scat_2 = scatterplot_plot(pass_cutoff_true_df[['ratio_GFAP_brain_Cre(+)_Ctrl_W7/GFAP_brain_Cre(-)_Ctrl_W6', 'ratio_GFAP_brain_Cre(+)_LPS_W9/GFAP_brain_Cre(-)_Ctrl_W6']], tissue_file_folder_path, "W7_vs_W9")
            df_scat_3 = scatterplot_plot(pass_cutoff_true_df[['ratio_GFAP_brain_Cre(+)_Fast_W8/GFAP_brain_Cre(-)_Ctrl_W6', 'ratio_GFAP_brain_Cre(+)_LPS_W9/GFAP_brain_Cre(-)_Ctrl_W6']], tissue_file_folder_path, "W8_vs_W9")

            scat_list = [df_scat_1, df_scat_2, df_scat_3]
            scat_df = pd.concat(scat_list, axis=1)
            
            pass_cutoff_true_df = pass_cutoff_true_df.join(scat_df)
            #pass_cutoff_true_df = pass_cutoff_true_df.drop(["pep_num", "annotation"], axis=1)
            pass_cutoff_true_df = pass_cutoff_true_df.reset_index().set_index(["uniprot_id", "description", "pep_num", "annotation"])
            pass_cutoff_true_df.to_csv(tissue_file_folder_path / ("GFAP_final_protein_table" + file_name.split("processed_census-out")[1] +'.csv'))
        else: 
            volcano_df = get_volcano_plot(conditions_list, control_labelling, treatment_labelling, pass_cutoff_df_norm_data, file_name, tissue_file_folder_path)
            pass_cutoff_true_df = pass_cutoff_true_df.join(volcano_df)
            pass_cutoff_true_df.to_csv(tissue_file_folder_path / ("final_protein_table" + file_name.split("processed_census-out")[1] +'.csv'))

# %%
