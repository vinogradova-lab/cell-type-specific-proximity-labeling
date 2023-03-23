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
import scipy.stats as stat
import plotly.express as px
from src.cutoff_funcs import *
from src.filter_funcs import *
from src.norm_funcs import *
from src.readin_funcs import * 

# %%
#provide folder paths
input_folder_path = Path("/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/01_census_out_processing/11_files_turbo_id_1pepperprotein_03062023/01_processed_files/")
output_folder_path  = Path('/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/032223_results/min_1_pep_per_protein')

#get master list
fasta_table = pd.read_csv("/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/02_mouse_TP_FP_table/main_fasta_table.csv", index_col=[0])
TP_list = fasta_table[fasta_table["annotation"] == "TP"]
#TP_list = TP_list["uniprot_id"].tolist()
FP_list = fasta_table[fasta_table["annotation"] == "FP"]
#FP_list = FP_list["uniprot_id"].tolist()

print("Number of Uniport IDs in TP list:", len(TP_list)) #2806
print("Number of Uniport IDs in FP list:", len(FP_list)) #1077

# %%
#Get list of files of processed census-out files
list_of_file_paths = list(input_folder_path.iterdir())
list_of_file_paths = [x for x in list_of_file_paths if 'census-out' in x.stem]

list_of_file_names = [file_path.stem for file_path in list_of_file_paths]
print("Number of files: ", len(list_of_file_names))
list_of_file_names

# %%
#get file_channel dict
file_channel_dict = get_channel_name_dict(input_folder_path, list_of_file_names)
len(file_channel_dict)

#%%
#get condition dict
file_condition_dict = get_cond_info(input_folder_path, list_of_file_names)
len(file_condition_dict)

# %%
#get list of keratins uniprot ids 
keratins = fasta_table[fasta_table.keratin == True]
keratins_list = keratins.index.tolist()
len(keratins_list)

# %%
#Assign new column names accoriding to metadata_col.csv and remove keratins in each file and annotate
folder_path = output_folder_path / "01_raw_files_with_correct_channel_names_keratins_removed_TP_FP_annotated"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

dfs_dict = {}
for file in list_of_file_paths:
    file_name = file.stem
    df = pd.read_csv(file)
    #remove keratins 
    print(len(df))
    df = df[~df.uniprot.isin(keratins_list)]
    df = df.set_index(["uniprot", 'description', 'pep_num'])
    #rename columns
    df.rename(columns=file_channel_dict[file_name], inplace=True)
    #annotate TP FP
    print(len(df))
    df = df.reset_index()
    df = df.rename(columns={"uniprot": "Entry"})
    merged = df.merge(fasta_table[["annotation", "Entry"]], on="Entry", how="left")
    merged = merged.drop_duplicates()
    merged = merged.rename(columns={"Entry": "uniprot_id"})
    #merged["annotation"].value_counts()
    print(len(merged))
    print("")
    merged_df = merged.set_index(["uniprot_id", "description", "pep_num", "annotation"])
    dfs_dict[file_name] = merged_df
    merged_df.to_csv(folder_path / (file.stem + ".csv"))


# %%
#Filter based on condition 
folder_path = output_folder_path / "02_filtering_per_cond_per_file"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

filtered_dict = {}
for file_name in list_of_file_names:
    #print(file_name)
    
    df = dfs_dict[file_name]
    conditions_list = file_condition_dict[file_name]["conditions"]
    control_labelling = file_condition_dict[file_name]["control_labelling"]
    treatment_labelling = file_condition_dict[file_name]["treatment_labelling"]

    filtered_sub_dfs = []
    for condition in conditions_list:
        sub_df = get_condition_df(df, condition)
        filtered_sub_df, cond_columns, ctrl_columns = filter_condition_df(sub_df, treatment_labelling, control_labelling)
        filtered_sub_df = filtered_sub_df.reset_index()
        filtered_sub_dfs.append(filtered_sub_df)
    print("Filtering for %s done" % (file_name))
    
    merge_on_cols = ['uniprot_id','description', 'pep_num', "annotation"]
    filtered_final_table = reduce(lambda df1,df2: pd.merge(df1,df2,on=merge_on_cols, how="outer"), filtered_sub_dfs)
    print(len(filtered_final_table))
    
    raw_df = df.reset_index()
    uniprot_ids = filtered_final_table["uniprot_id"].tolist()
    sub_raw_df = raw_df[raw_df['uniprot_id'].isin(uniprot_ids)]
    print(len(sub_raw_df))
    
    sub_raw_df = sub_raw_df.set_index(["uniprot_id", 'description', 'pep_num', 'annotation'])

    filtered_dict[file_name] = sub_raw_df
    sub_raw_df.to_csv(folder_path / (file_name + ".csv"))

# %%
# (1) generate the output put with raw SI and annotation TP/FP; (2) look at TP proteins with 2+ peptides; (3) calculate median SI for each cre+ channel and median of sums; (4) calculate normalization ratios by dividing (median of median)/(median cre+ channel), and (5) use those normalization values for each channel;
normalized_dict = {}

folder_path = output_folder_path / "03_normalisation_plots"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

for file_name in list_of_file_names:
    print(file_name)
    plot_list = []
    
    fig = plt.figure(figsize=(25, 13)) 
    fig.suptitle(file_name)

    ax1 = plt.subplot(121) 
    ax2 = plt.subplot(122) 
    
    df = filtered_dict[file_name]
    conditions_list = file_condition_dict[file_name]["conditions"]
    control_labelling = file_condition_dict[file_name]["control_labelling"]
    treatment_labelling = file_condition_dict[file_name]["treatment_labelling"]
    #turbo_id_value_per_channel = pd.Series(file_turbo_id_dict[file_name])
    
    #RAW
    #plots for raw df
    #message = pearson_corr_channels(df, file_name, folder_path, "Raw")
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
    
    
    #MEDIAN OF MEDIANS TP only, 2+ peptides
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
        #sub_df = sub_df.head(100)
        sub_df = sub_df.drop("median_"+condition+"_"+treatment_labelling, axis=1)
        sub_df = sub_df.drop("median_"+condition+"_"+control_labelling, axis=1)
        sub_df = sub_df.drop("median_ratio", axis=1)
        sub_df = sub_df.set_index(["uniprot_id", "description", "pep_num", "annotation"])
        sub_df = sub_df[cond_columns]
        norm_factors = sub_df.median().median() / sub_df.median()
        norm_factor_series.append(norm_factors)
    norm_factor = pd.concat(norm_factor_series, axis=1).sum(1)
    print(norm_factor)
    normalized_df = df[norm_factor.index] * norm_factor
    list_of_ctrl_cols = [column for column in df.columns if column not in norm_factor.index]
    untouched_ctrl_columns = df[list_of_ctrl_cols]
    norm_df = normalized_df.reset_index().merge(untouched_ctrl_columns.reset_index(), on=["uniprot_id", "description", "pep_num", "annotation"])
    norm_df = norm_df.set_index(["uniprot_id", "description", "pep_num", "annotation"])
    norm_df = norm_df[df.columns] #get the same order of columns 
    #print(len(norm_df) == len(df))
    #print(norm_df.columns == df.columns)

    #message = pearson_corr_channels(norm_df, file_name, folder_path, "Normalized")
    pca_fig_6 = get_pca_plot(norm_df, "Median of medians TP only, 2+ peptides, Cre+")
    plot_list.append(pca_fig_6)
    
    norm_df_log = np.log2(norm_df)
    norm_df_log.plot(kind='box', 
                     rot=90, 
                     fontsize=10, 
                     title="Median of medians TP only, 2+ peptides, Cre+", 
                     ax=ax2, 
                     xlabel="Channels", 
                     ylabel="log2(SI)", 
                     color = "darkblue")
    
    #this is what gets passed into cutoff analysis
    normalized_dict[file_name] = norm_df
    
    fig.savefig(folder_path / (file_name+'box_plots.png'))
    with open(folder_path / (file_name+'pca_plots.html') , 'w') as f:
        f.write(file_name)
        for fig in plot_list:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))

# %%
#save normalized df in folder
folder_path = output_folder_path / "04_normalised"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

for file_name in list_of_file_names:
    df = normalized_dict[file_name]
    print(len(df))
    print("")
    df.to_csv(folder_path / (file_name + ".csv"))


# %%
#Cutoff analysis of all files
folder_path = output_folder_path / "05_cutoff_analysis"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)
    
cutoff_dict = {}
for file_name in list_of_file_names:
    print(file_name)
    
    #create output folder
    folder_name = file_name.replace("processed_census-out_", "")
    folder_path = output_folder_path / "05_cutoff_analysis" / folder_name
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
        
    #get data
    df = normalized_dict[file_name]
    #df = filtered_dict[file_name]
    conditions_list = file_condition_dict[file_name]["conditions"]
    control_labelling = file_condition_dict[file_name]["control_labelling"]
    treatment_labelling = file_condition_dict[file_name]["treatment_labelling"]
    
    #filtering and get ratio
    ratio_dfs_dict = {}

    for condition in conditions_list:
        sub_df = get_condition_df(df, condition)
        sub_df = sub_df.dropna() #this step important because we want to only take the ones in account which pass cond filtering
        cond_columns = sub_df.filter(like=treatment_labelling).columns.tolist()
        ctrl_columns = sub_df.filter(like=control_labelling).columns.tolist()
        ratio_df = get_ratio_condition_df(sub_df, cond_columns, ctrl_columns)
        ratio_dfs_dict[condition] = ratio_df
        
    ratio_tables = list(ratio_dfs_dict.values())
    ratio_tables = [table.reset_index() for table in ratio_tables]
    col_list = ['uniprot_id', 'description', 'pep_num', 'annotation']
    ratio_tables_merged = reduce(lambda df1,df2: pd.merge(df1,df2,on=col_list, how="outer"), ratio_tables)
    #merge with normalized SI
    ratio_and_signal_intensity = ratio_tables_merged.merge(df, on=['uniprot_id', 'description', 'pep_num'], how='left')
    
    all_cond_full_tables = {}
    all_cond_cutoff_tables = {}
    
    sub_cutoff_dict = {}
    all_cond_cutoff_tables = get_cutoffs(ratio_dfs_dict, folder_path, all_cond_cutoff_tables, sub_cutoff_dict)
    #break
    cutoff_dict[file_name] = sub_cutoff_dict
    
    #For each replicate get number of proteins that pass cutoff
    for condition, table_list in all_cond_cutoff_tables.items():
        for table in table_list:
            col_name = table.filter(like='log2_norm_ratio_').columns.tolist()
            col_name = col_name[0]
            col_name = col_name.replace("log2_norm_ratio_", "")
            #total_TP = len(table[table["annotation"] == "TP"])
            #total_FP = len(table[table["annotation"] == "FP"])
            table.drop('TP', axis=1, inplace=True)
            table.drop('FP', axis=1, inplace=True)
            table.drop('TPR', axis=1, inplace=True)
            table.drop('FPR', axis=1, inplace=True)
            table.rename(columns={"pass_cutoff": "pass_cutoff_"+col_name,
                                  "TPR-FPR": "TPR-FPR_"+col_name}, inplace=True)
            
    #Merge replicates of the same condition into one table and get everything in one table (outer merge)
    merged_replicate_tables = []
    for condition, table_list in all_cond_cutoff_tables.items():
        df = reduce(lambda df1,df2: pd.merge(df1,df2,on=['uniprot_id', 'annotation'], how="outer"), table_list)
        merged_replicate_tables.append(df)
        #total_TP = len(df[df["annotation"] == "TP"])
        #total_FP = len(df[df["annotation"] == "FP"])
            
    
    #Merge cutoff ratio from all conditions in one file
    df_all = reduce(lambda df1,df2: pd.merge(df1,df2,on=['uniprot_id', 'annotation'], how="outer"), merged_replicate_tables)
    #df_all.to_csv(folder_path / ("log2_FPR_TPR.csv"), index=False)
    
    #merge ratio and signal intensity with cutoff values
    ratio_and_signal_intensity.set_index(['uniprot_id', 'description', 'pep_num', 'annotation'], inplace=True)
    
    for condition in conditions_list:
        ratio_condition_cols = ratio_and_signal_intensity.filter(like="ratio_"+condition).columns.tolist()
        ratio_and_signal_intensity["median_R("+condition+")"] = ratio_and_signal_intensity[ratio_condition_cols].median(axis=1)
    
    column_list = ratio_and_signal_intensity.columns.tolist()
    column_list = [colname for colname in column_list if "ratio" not in colname]
    
    for condition in conditions_list:
        conditions_cols_pos = [colname for colname in column_list if condition+"_"+treatment_labelling in colname]
        if len(conditions_cols_pos) == 0:
            conditions_cols_pos = [colname for colname in column_list if condition+treatment_labelling in colname]
        ratio_and_signal_intensity["median_SI("+condition+")"] = ratio_and_signal_intensity[conditions_cols_pos].median(axis=1)

    ratio_and_signal_intensity = ratio_and_signal_intensity.reset_index()
    
    #merge 
    df_all_list_columns_TPRFPR = df_all.filter(like="TPR-FPR").columns.tolist()
    df_all_list_columns_passcutoff = df_all.filter(like="pass_cutoff").columns.tolist()

    columns_list_df_all = df_all_list_columns_TPRFPR + df_all_list_columns_passcutoff
    columns_list_df_all.append("uniprot_id")
    
    #reduced_fasta_table = fasta_table.drop(fasta_table.columns[[0,-10,-9,-8,-7,-6]], axis = 1)
    
    ratio_and_signal_intensity = ratio_and_signal_intensity.rename(columns={"uniprot_id": "Entry"})
    ratio_and_signal_intensity.drop('description', axis=1, inplace=True)
    fasta_table = fasta_table.rename(columns={"uniprot_id": "alias_uniprot_id"})

    merged_with_metadata = pd.merge(ratio_and_signal_intensity, fasta_table, on=["Entry", "annotation"], how="left")
    merged_with_metadata = merged_with_metadata.rename(columns={"Entry": "uniprot_id"})
    merged_with_metadata = merged_with_metadata.drop_duplicates(subset='uniprot_id', keep='first')
    
    ratio_and_signal_intensity_merged = pd.merge(merged_with_metadata, df_all[columns_list_df_all],on=["uniprot_id"], how="outer")
    
    print(len(ratio_and_signal_intensity_merged))
    print(len(df_all))
    
    print("")
    ratio_and_signal_intensity_merged = add_pass_cutoff_analysis_to_df(ratio_and_signal_intensity_merged, df_all_list_columns_passcutoff)
    
    with pd.ExcelWriter(folder_path / (folder_name+'_results.xlsx')) as writer:
        ratio_and_signal_intensity_merged.to_excel(writer, sheet_name='ratio_raw_values', index=False)
        df_all.to_excel(writer, sheet_name='log2_norm_ratio', index=False)


# %%
#Analysis of Proteins that pass cutoff Filter (i.e. pass_cutoff_result=TRUE)
input_folder_path = output_folder_path / "05_cutoff_analysis" 
list_of_dir_paths = list(input_folder_path.iterdir())

list_of_result_file_paths = []

for directory in list_of_dir_paths:
    list_in_dir = list(directory.iterdir())
    for file_path in list_in_dir: 
        if "results" in file_path.stem:
            list_of_result_file_paths.append(file_path)

# %%
folder_path = output_folder_path / "06_results" 
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

pca_fig_list = []

for file_path in list_of_result_file_paths: 
    barplot_list = []
    if "~$" not in file_path.stem:
        sheet_1 = pd.read_excel(file_path, engine="openpyxl", sheet_name='ratio_raw_values') 
        sheet_2 = pd.read_excel(file_path, engine="openpyxl", sheet_name='log2_norm_ratio')
        
        #has to pass in at least 2 different cre+
        #list_columns_passcutoff = sheet_1.filter(like="pass_cutoff").columns.tolist()
        #list_columns_passcutoff = [e for e in sheet_1.columns if e.startswith('pass_cutoff')]
        #list_columns_passcutoff.pop()
        #sheet_1 = add_pass_cutoff_analysis_to_df(sheet_1, list_columns_passcutoff)
         
    
        for key, values in file_channel_dict.items():
            if file_path.stem.replace("_results", "") in key:
                col_list = list(file_channel_dict[key].values())
    
        pass_cutoff_true_df = sheet_1[sheet_1.pass_cutoff_result == True]
        pass_cutoff_true_df = pass_cutoff_true_df.set_index(["uniprot_id", "annotation", "Signal peptide"])
        pass_cutoff_true_df = pass_cutoff_true_df[col_list]
        
        title = file_path.stem.replace("_results", "") + " - min 2 cutoff/cond"
    
        pca_fig = get_pca_plot_after(pass_cutoff_true_df, title)
        pca_fig_list.append(pca_fig)
        
        pass_cutoff_true_df = pass_cutoff_true_df.reset_index()
        total = len(pass_cutoff_true_df)
        tp_num = len(pass_cutoff_true_df[pass_cutoff_true_df["annotation"] == "TP"])
        fp_num = len(pass_cutoff_true_df[pass_cutoff_true_df["annotation"] == "FP"])
        na_num = len(pass_cutoff_true_df[pass_cutoff_true_df["annotation"].isnull()])
        
        pass_cutoff_true_df_list = ["pass cutoff 2 \n ("+ str(total)+" proteins)", (tp_num/total)*100, (fp_num/total)*100, (na_num/total)*100]
        
        sheet_1 = pd.read_excel(file_path, engine="openpyxl", sheet_name='ratio_raw_values') 
        list_columns_ratio = sheet_1.filter(like="ratio_").columns.tolist()
        sheet_1 = sheet_1.set_index(["uniprot_id", "annotation"])
        sheet_1 = sheet_1[list_columns_ratio].dropna()
        sheet_1 = sheet_1.reset_index()
        #sheet_1 = sheet_1[["uniprot", "annotation"]]
        total = len(sheet_1)
        tp_num = len(sheet_1[sheet_1["annotation"] == "TP"])
        fp_num = len(sheet_1[sheet_1["annotation"] == "FP"])
        na_num = len(sheet_1[sheet_1["annotation"].isnull()])
        
        all_df = ["all \n ("+ str(total)+" proteins)", (tp_num/total)*100, (fp_num/total)*100, (na_num/total)*100]
        
        barplot_list.append(all_df)
        barplot_list.append(pass_cutoff_true_df_list)
        
        # Create the pandas DataFrame
        df = pd.DataFrame(barplot_list, columns=['state', 'True positives', 'False positives', 'Not annotated'])
        df = df.set_index("state")
        df.plot(kind='bar', stacked=True, color=['green','red','lightgrey'])
        plt.xlabel('')
        plt.ylabel("% Proteins")
        plt.title(title)
        plt.xticks(rotation = 0) 
        #plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        plt.savefig(folder_path / (file_path.stem.replace("_results", "")+'barplot.png'))
        plt.show()

with open(folder_path / 'pca_plots.html' , 'w') as f:
    for fig in pca_fig_list:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))

# %%
for file_path in list_of_result_file_paths: 
    if "~$" not in file_path.stem:
        sheet_2 = pd.read_excel(file_path, engine="openpyxl", sheet_name='log2_norm_ratio')
        sheet_2 = sheet_2.set_index(["uniprot_id", "annotation"])
        list_columns_passcutoff = sheet_2.filter(like="pass_cutoff").columns.tolist()
        
        folder_path = output_folder_path / "06_results" / file_path.stem
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)
            
        for col in list_columns_passcutoff:
            df_list = []
            col_clean = col.replace("pass_cutoff_", "")
    
            for file, cutoff_vals in cutoff_dict.items():
                if file_path.stem.replace("_results", "") in file:
                    cutoff_value = cutoff_vals[col_clean]
                    cutoff_value = round(cutoff_value, 3)   
        
            sub_df_cols = sheet_2.filter(like=col_clean).columns.tolist()
            sub_df = sheet_2[sub_df_cols]
            #print(len(sub_df))
            sub_df = sub_df.dropna(subset=["log2_norm_ratio_"+col_clean])
            #print(len(sub_df))
            sub_df = sub_df.reset_index()
            sub_df = sub_df.set_index("uniprot_id")
            sub_df_plot = sub_df[["log2_norm_ratio_"+col_clean, "annotation"]]
            TP = sub_df_plot[sub_df_plot["annotation"] == "TP"].drop("annotation", axis=1).iloc[:, 0].rename("TP")
            FP = sub_df_plot[sub_df_plot["annotation"] == "FP"].drop("annotation", axis=1).iloc[:, 0].rename("FP")
            no_val = sub_df_plot[sub_df_plot["annotation"].isnull()].drop("annotation", axis=1).iloc[:, 0].rename("NA")
    
            all_list = ["all \n ("+ str(len(sub_df_plot))+" proteins)", (len(TP)/len(sub_df_plot))*100, (len(FP)/len(sub_df_plot))*100, (len(no_val)/len(sub_df_plot))*100]
            df_list.append(all_list)
    
            df_to_plot = pd.concat([TP,FP, no_val],axis=1)
            gfg = sns.histplot(df_to_plot, palette=dict(TP="green", FP="red", NA="grey"))
            plt.axvline(cutoff_value, 0, color="black", linewidth=0.5, linestyle='dashed')
            plt.text(cutoff_value+1, 100, "Cutoff: "+str(cutoff_value))
            gfg.set(xlabel ="log2("+col_clean+")", title = file_path.stem+ "\n" + col_clean)
            plt.show()
            fig = gfg.get_figure()
            description = col_clean.replace("/", "_")
            fig.savefig(folder_path / (description+'hist_1.png')) 
    
            df_to_plot = pd.concat([TP,FP],axis=1)
            gfg_2 = sns.histplot(df_to_plot, palette=dict(TP="green", FP="red"))
            plt.axvline(cutoff_value, 0, color="black", linewidth=0.5, linestyle='dashed')
            plt.text(cutoff_value+1, 60, "Cutoff: "+str(cutoff_value))
            gfg_2.set(xlabel ="log2("+col_clean+")", title = file_path.stem+ "\n" + col_clean)
            plt.show()
            fig = gfg_2.get_figure()
            description = col_clean.replace("/", "_")
            fig.savefig(folder_path / (description+'hist_2.png')) 
 
    
            only_pass_cutoff = sub_df[sub_df[col] == 1]
            total = len(only_pass_cutoff)
            tp_num = len(only_pass_cutoff[only_pass_cutoff["annotation"] == "TP"])
            fp_num = len(only_pass_cutoff[only_pass_cutoff["annotation"] == "FP"])
            na_num = len(only_pass_cutoff[only_pass_cutoff["annotation"].isnull()])
            
            if total != 0:
                pass_cutoff_list = ["pass cutoff \n("+ str(total)+" proteins)", (tp_num/total)*100, (fp_num/total)*100, (na_num/total)*100]
            elif total == 0:
                pass_cutoff_list = ["pass cutoff \n("+ str(total)+" proteins)", 0, 0, 0]
            
            df_list.append(pass_cutoff_list)
    
            only_notpass_cutoff = sub_df[sub_df[col] == 0]
            total = len(only_notpass_cutoff)
            tp_num = len(only_notpass_cutoff[only_notpass_cutoff["annotation"] == "TP"])
            fp_num = len(only_notpass_cutoff[only_notpass_cutoff["annotation"] == "FP"])
            na_num = len(only_notpass_cutoff[only_notpass_cutoff["annotation"].isnull()])
            
            if total != 0:
                only_notpass_cutoff_list = ["not pass cutoff \n ("+ str(total)+" proteins)", (tp_num/total)*100, (fp_num/total)*100, (na_num/total)*100]
            elif total == 0:
                only_notpass_cutoff_list = ["not pass cutoff \n ("+ str(total)+" proteins)", 0, 0, 0]

                
            df_list.append(only_notpass_cutoff_list)
    
            # Create the pandas DataFrame
    
            df = pd.DataFrame(df_list, columns=['state', 'True positives', 'False positives', 'Not annotated'])
            df = df.set_index("state")
            df.plot(kind='bar', stacked=True, color=['green','red','lightgrey'])
            plt.xlabel('')
            plt.ylabel("% Proteins")
            plt.title(file_path.stem+ "\n" +col_clean )
            plt.xticks(rotation = 0) 
            #plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
            description = col_clean.replace("/", "_")
            plt.savefig(folder_path / (description+'barplot.png'))
            plt.show()

# %% [markdown]
# **Volcano plots**
# 
# Plan: 
# - mean of replicates per group 
# - fold changes of these groups 
# - log2 transform fold change 
# - stat.ttest_ind 
# - Bonferroni correction 
# - -log10 transform 
# - plot data 
# 
# https://thecodingbiologist.com/posts/Making-volcano-plots-in-python-in-Google-Colab
# 

# %%
folder_path = output_folder_path / "06_results" 
volcano_plot_list_10plex = []
volcano_plot_list_16plex = []

for file_path in list_of_result_file_paths: 
    barplot_list = []
    if "~$" not in file_path.stem:
        sheet_1 = pd.read_excel(file_path, engine="openpyxl", sheet_name='ratio_raw_values') 
        sheet_2 = pd.read_excel(file_path, engine="openpyxl", sheet_name='log2_norm_ratio')
        
        for key, values in file_channel_dict.items():
            if file_path.stem.replace("_results", "") in key:
                col_list = list(file_channel_dict[key].values())
        
        pass_cutoff_true_df = sheet_1[sheet_1.pass_cutoff_result == True]
        pass_cutoff_true_df = pass_cutoff_true_df.set_index(["uniprot_id", "annotation", "Signal peptide", "Entry Name"])
        cols_df = pass_cutoff_true_df[col_list]
        
        conditions_list = file_condition_dict["processed_census-out_"+file_path.stem.replace("_results", "")]["conditions"]
        control_labelling = file_condition_dict["processed_census-out_"+file_path.stem.replace("_results", "")]["control_labelling"]
        treatment_labelling = file_condition_dict["processed_census-out_"+file_path.stem.replace("_results", "")]["treatment_labelling"]
        
        if len(conditions_list) == 2: 
            list_columns_trt = cols_df.filter(like=treatment_labelling).columns.tolist()
            trt_col_df = cols_df[list_columns_trt]
            trt_col_df.shape
    
            list_cond_1 = trt_col_df.filter(like=conditions_list[0]).columns.tolist()
            list_cond_2 = trt_col_df.filter(like=conditions_list[1]).columns.tolist()
            idx_cond_1 = trt_col_df.columns.get_indexer(list_cond_1)
            idx_cond_2 = trt_col_df.columns.get_indexer(list_cond_2)
            trt_col_df["p_value"] = trt_col_df.apply(get_p_value, axis=1, args=(list_cond_1, idx_cond_2), )
            #-log10 
            #trt_col_df["-log10_pval"] = -1*np.log10(trt_col_df["p_value"])
            #Bonferroni correction
            #trt_col_df["-log10_pval"] = -1*np.log10(len(trt_col_df)*trt_col_df["p_value"])
            for condition in conditions_list:
                list_columns_cond = trt_col_df.filter(like=condition).columns.tolist()
                cond_df = trt_col_df[list_columns_cond]

                trt_col_df["mean_"+condition+"_"+treatment_labelling] = cond_df.mean(axis=1)
    
            trt_col_df["FC"] = trt_col_df["mean_"+conditions_list[0]+"_"+treatment_labelling] / trt_col_df["mean_"+conditions_list[1]+"_"+treatment_labelling]
            trt_col_df["log2_FC"] = np.log2(trt_col_df["FC"])

            volcano_df = trt_col_df[["p_value", "log2_FC"]]
            
            #print(volcano_df.shape)
            volcano_df = volcano_df.dropna() 
            #Bonferroni correction
            #volcano_df["pval_bonferroni_corr"] = volcano_df["p_value"] * len(volcano_df)
            #-log10
            #volcano_df["-log10_pval"] = -1*np.log10(volcano_df["pval_bonferroni_corr"])
            #volcano_df["-log10_pval"] = -1*np.log10(len(volcano_df)*volcano_df["p_value"])
            volcano_df["-log10_pval"] = -1*np.log10(volcano_df["p_value"])
            
            volcano_df["Regulation"] = volcano_df.apply(get_expr, axis=1)
            
            volcano_df["Regulation"] = volcano_df["Regulation"].astype('category')
            my_order = ["Down", "Stable", "Up"]
            volcano_df['Regulation'].cat.reorder_categories(my_order, inplace= True)
            
            
            #print(volcano_df.shape)
            volcano_df = volcano_df.reset_index()
            volcano_df = volcano_df.set_index(["annotation", "Signal peptide"])
            volcano_df["name"] = volcano_df['Entry Name'].str.split("_").str[0]
        
            x_axis_name = "log2(" + conditions_list[0] +"/" + conditions_list[1] +")"
            title_name = file_path.stem.replace("_results", "") + " - " + conditions_list[0] +" vs. " + conditions_list[1] + " ("+ str(len(volcano_df)) + " Proteins)"

            fig = px.scatter(volcano_df, 
                             x='log2_FC', 
                             y='-log10_pval',
                             color='Regulation',
                             color_discrete_sequence=["red", 'grey', 'green'],
                             hover_data=['name', 'uniprot_id'],
                             title = title_name,
                             labels = {"log2_FC": x_axis_name},
                             template = "simple_white",
                             category_orders={'Regulation': np.sort(volcano_df['Regulation'].unique())})
            fig.add_vline(x=0.58, line_width=2, line_dash="dash", line_color="grey")
            fig.add_vline(x=-0.58, line_width=2, line_dash="dash", line_color="grey")
            fig.add_hline(y=1.3, line_width=2, line_dash="dash", line_color="grey")
            fig.update_layout(legend=dict(title=""), title_x=0.5)
            fig.show()
            
            volcano_plot_list_10plex.append(fig)
            #print(volcano_df.head())
        
        elif len(conditions_list) == 3: 
            list_columns_trt = cols_df.filter(like=treatment_labelling).columns.tolist()
            trt_col_df = cols_df[list_columns_trt]
    
            list_cond_1 = trt_col_df.filter(like=conditions_list[0]).columns.tolist()
            list_cond_2 = trt_col_df.filter(like=conditions_list[1]).columns.tolist()
            list_cond_3 = trt_col_df.filter(like=conditions_list[2]).columns.tolist()
            
            idx_cond_1 = trt_col_df.columns.get_indexer(list_cond_1)
            idx_cond_2 = trt_col_df.columns.get_indexer(list_cond_2)
            idx_cond_3 = trt_col_df.columns.get_indexer(list_cond_3)
            
            trt_col_df["p_value_1"] = trt_col_df.apply(get_p_value, axis=1, args=(list_cond_1, idx_cond_2), )
            trt_col_df["p_value_2"] = trt_col_df.apply(get_p_value, axis=1, args=(list_cond_1, idx_cond_3), )
            trt_col_df["p_value_3"] = trt_col_df.apply(get_p_value, axis=1, args=(list_cond_2, idx_cond_3), )
            
            for condition in conditions_list:
                list_columns_cond = trt_col_df.filter(like=condition).columns.tolist()
                cond_df = trt_col_df[list_columns_cond]
                trt_col_df["mean_"+condition+"_"+treatment_labelling] = cond_df.mean(axis=1)
    
            trt_col_df["FC_1"] = trt_col_df["mean_"+conditions_list[1]+"_"+treatment_labelling] / trt_col_df["mean_"+conditions_list[0]+"_"+treatment_labelling]
            trt_col_df["FC_2"] = trt_col_df["mean_"+conditions_list[2]+"_"+treatment_labelling] / trt_col_df["mean_"+conditions_list[0]+"_"+treatment_labelling]
            trt_col_df["FC_3"] = trt_col_df["mean_"+conditions_list[1]+"_"+treatment_labelling] / trt_col_df["mean_"+conditions_list[2]+"_"+treatment_labelling]
            
            trt_col_df["log2_FC_1"] = np.log2(trt_col_df["FC_1"])
            trt_col_df["log2_FC_2"] = np.log2(trt_col_df["FC_2"])
            trt_col_df["log2_FC_3"] = np.log2(trt_col_df["FC_3"])
            
            df_list = []
            volcano_df_1 = trt_col_df[["p_value_1", "log2_FC_1"]]
            df_list.append(volcano_df_1)
            volcano_df_2 = trt_col_df[["p_value_2", "log2_FC_2"]]
            df_list.append(volcano_df_2)
            volcano_df_3 = trt_col_df[["p_value_3", "log2_FC_3"]]
            df_list.append(volcano_df_3)
            
            for volcano_df in df_list:
                if "p_value_1" in volcano_df.columns:
                    x_axis_name = "log2(" + conditions_list[1] +"/" + conditions_list[0] +")"
                    title_name = file_path.stem.replace("_results", "") + " - " + conditions_list[1] +" vs. " + conditions_list[0] + " ("+ str(len(volcano_df)) + " Proteins)"
                elif "p_value_2" in volcano_df.columns:
                    x_axis_name = "log2(" + conditions_list[2] +"/" + conditions_list[0] +")"
                    title_name = file_path.stem.replace("_results", "") + " - " + conditions_list[2] +" vs. " + conditions_list[0] + " ("+ str(len(volcano_df)) + " Proteins)"
                elif "p_value_3" in volcano_df.columns:
                    x_axis_name = "log2(" + conditions_list[1] +"/" + conditions_list[2] +")"
                    title_name = file_path.stem.replace("_results", "") + " - " + conditions_list[1] +" vs. " + conditions_list[2] + " ("+ str(len(volcano_df)) + " Proteins)"
            
                volcano_df.columns = ["p_value", "log2_FC"]
                volcano_df = volcano_df.dropna() 
                volcano_df["-log10_pval"] = -1*np.log10(volcano_df["p_value"])
                volcano_df["Regulation"] = volcano_df.apply(get_expr, axis=1)
            
                volcano_df["Regulation"] = volcano_df["Regulation"].astype('category')
                my_order = ["Down", "Stable", "Up"]
                #volcano_df['Regulation'].cat.reorder_categories(my_order, inplace= True)
            
                volcano_df = volcano_df.reset_index()
                volcano_df = volcano_df.set_index(["annotation", "Signal peptide"])
                volcano_df["name"] = volcano_df['Entry Name'].str.split("_").str[0]
                
                fig = px.scatter(volcano_df, 
                                 x='log2_FC', 
                                 y='-log10_pval',
                                 color='Regulation',
                                 color_discrete_sequence=["red", 'grey', 'green'],
                                 hover_data=['name', 'uniprot_id'],
                                 title = title_name,
                                 labels = {"log2_FC": x_axis_name},
                                 template = "simple_white",
                                 category_orders={'Regulation': np.sort(volcano_df['Regulation'].unique())})
                fig.add_vline(x=0.58, line_width=2, line_dash="dash", line_color="grey")
                fig.add_vline(x=-0.58, line_width=2, line_dash="dash", line_color="grey")
                fig.add_hline(y=1.3, line_width=2, line_dash="dash", line_color="grey")
                fig.update_layout(legend=dict(title=""), title_x=0.5)
                fig.show()
            
                volcano_plot_list_16plex.append(fig)

with open(folder_path / 'volcano_plots_10plex.html' , 'w') as f:
    for fig in volcano_plot_list_10plex:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
        
with open(folder_path / 'volcano_plots_16plex.html' , 'w') as f:
    for fig in volcano_plot_list_16plex:
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
        



