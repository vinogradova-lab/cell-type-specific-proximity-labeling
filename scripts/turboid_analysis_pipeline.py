# %% [markdown]
# **Import libraries**

# %%
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import scipy
from scipy import stats
import seaborn as sns
from functools import reduce
import pathlib
from csv import reader
%matplotlib inline
import warnings
import xlsxwriter
from sklearn.decomposition import PCA
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from scipy.stats import pearsonr
import scipy.stats as stat
import dash_bio
import plotly.express as px
warnings.filterwarnings('ignore')

# %% [markdown]
# **User input**

# %%
#provide folder paths
input_folder_path = "/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/TurboID_analysis/analysis_results/11_files_turbo_id_2pepperprotein_03062023/01_processed_files/"
output_path = '/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/TurboID_analysis/notebook_analysis/06_11files_030623/min_2_peptide_per_protein'

#get master list
master_list = pd.read_csv("/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/TurboID_analysis/mouse_lists/master_list/mouse_master_table_for_ROC_analysis.csv", index_col=[0])
TP_list = master_list[master_list["annotation"] == "TP"]
#TP_list = TP_list["uniprot_id"].tolist()
FP_list = master_list[master_list["annotation"] == "FP"]
#FP_list = FP_list["uniprot_id"].tolist()

print("Number of Uniport IDs in TP list:", len(TP_list)) #3436
print("Number of Uniport IDs in FP list:", len(FP_list)) #4353

# %% [markdown]
# **Change path into path object**

# %%
#turn path into object 
input_path_obj = pathlib.Path(input_folder_path)
#get absolute path
input_path_obj = input_path_obj.resolve()

output_path_obj = pathlib.Path(output_path)
#get absolute path
output_path_obj = output_path_obj.resolve()

# %% [markdown]
# **Get list of files in input folder**

# %%
list_of_file_paths = list(input_path_obj.iterdir())
list_of_file_paths = [x for x in list_of_file_paths if 'census-out' in x.stem]

list_of_file_names = []
for file_path in list_of_file_paths:
    file_name = file_path.stem
    list_of_file_names.append(file_name)
list_of_file_names

# %% [markdown]
# **Assign new column names accoriding to metadata_col.csv and remove keratins in each file and annotate**

# %%
meta_data_list = []
with open(input_path_obj / "metadata_col.csv", 'r') as read_obj:
    csv_reader = reader(read_obj)
    for row in csv_reader:
        meta_data_list.append(row)

# %%
file_channel_dict = {}
file_turbo_id_dict = {}
for file_name in list_of_file_names:
    channel_dict = {}
    turboid_dict = {}
    for meta_data_item in meta_data_list[1:]:
        if file_name in meta_data_item:
            channel_dict[meta_data_item[1]] = meta_data_item[2]
            turboid_dict[meta_data_item[2]] = float(meta_data_item[4])
    file_channel_dict[file_name] = channel_dict
    file_turbo_id_dict[file_name] = turboid_dict

# %%
#get list of keratins uniprot ids 
keratins = master_list[master_list.keratin == True]
keratins_list = keratins.uniprot_id.tolist()

# %%
#save new dfs in folder 
folder_path = output_path_obj / "01_raw_files_with_correct_channel_names_keratins_removed_TP_FP_annotated"
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
    merged = df.merge(master_list[["annotation", "Entry"]], on="Entry", how="left")
    merged = merged.drop_duplicates()
    merged = merged.rename(columns={"Entry": "uniprot_id"})
    #merged["annotation"].value_counts()
    print(len(merged))
    print("")
    merged_df = merged.set_index(["uniprot_id", "description", "pep_num", "annotation"])
    dfs_dict[file_name] = merged_df
    merged_df.to_csv(folder_path / (file.stem + ".csv"))

# %% [markdown]
# **Get condition information from conditions_metadata.csv**

# %%
condition_list = []
with open(input_path_obj / "conditions_metadata.csv", 'r') as read_obj:
    csv_reader = reader(read_obj)
    for row in csv_reader:
        condition_list.append(row)

# %%
file_condition_dict = {}
for file_name in list_of_file_names:
    condition_dict = {}
    for file_row in condition_list[1:]:
        if file_name in file_row:
            condition_dict["conditions"] = file_row[1].replace(" ", "").split(",")
            condition_dict["control_labelling"] = file_row[2]
            condition_dict["treatment_labelling"] = file_row[3]
    file_condition_dict[file_name] = condition_dict

# %% [markdown]
# **Filter based on condition**

# %%
#Functions for Filtering
def get_condition_df(df, condition_name):
    sub_df = df.filter(regex=condition_name)
    return(sub_df)

def filter_condition_df(sub_df, condition_cols_contain, control_cols_contain):
    cond_columns = sub_df.filter(like=condition_cols_contain).columns.tolist()
    ctrl_columns = sub_df.filter(like=control_cols_contain).columns.tolist()
    
    if len(cond_columns) == 3:
        column_pairs = [[cond_columns[0], cond_columns[1]], 
                        [cond_columns[0], cond_columns[2]], 
                        [cond_columns[1], cond_columns[2]]]
    elif len(cond_columns) == 4:
        column_pairs = [[cond_columns[0], cond_columns[1]], 
                        [cond_columns[0], cond_columns[2]], 
                        [cond_columns[0], cond_columns[3]], 
                        [cond_columns[1], cond_columns[2]], 
                        [cond_columns[1], cond_columns[3]], 
                        [cond_columns[2], cond_columns[3]]]
    elif len(cond_columns) == 5:
        column_pairs = [[cond_columns[0], cond_columns[1]], 
                        [cond_columns[0], cond_columns[2]], 
                        [cond_columns[0], cond_columns[3]], 
                        [cond_columns[0], cond_columns[4]], 
                        [cond_columns[1], cond_columns[2]], 
                        [cond_columns[1], cond_columns[3]], 
                        [cond_columns[1], cond_columns[4]], 
                        [cond_columns[2], cond_columns[3]],
                        [cond_columns[2], cond_columns[4]],
                        [cond_columns[3], cond_columns[4]]]
        
    list_of_uniprots = []
    for columnpair in column_pairs: 
        filter_df = sub_df[columnpair]
        filter_df["sum"] = filter_df[columnpair].sum(axis=1)
        filter_df["cv"] = filter_df[columnpair].std(axis=1).div(filter_df[columnpair].mean(axis=1))
        filter_df = filter_df[(filter_df["sum"] >= 10000) & (filter_df["cv"] <= 0.5)]
        filter_df = filter_df.reset_index()
        filter_df_ids = filter_df.uniprot_id.tolist()
        [list_of_uniprots.append(item) for item in filter_df_ids if item not in list_of_uniprots]
            
    sub_df = sub_df.reset_index()
    sub_df = sub_df[sub_df["uniprot_id"].isin(list_of_uniprots)]
    sub_df = sub_df.set_index(["uniprot_id", 'description', 'pep_num', 'annotation'])
    
    return(sub_df, cond_columns, ctrl_columns)


# %%
#save new dfs in folder 
folder_path = output_path_obj / "02_filtering_per_cond_per_file"
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

# %% [markdown]
# **Normalisation**
# 
# (1) generate the output put with raw SI and annotation TP/FP; (2) look at TP proteins with 2+ peptides; (3) calculate median SI for each cre+ channel and median of sums; (4) calculate normalization ratios by dividing (median of median)/(median cre+ channel), and (5) use those normalization values for each channel;

# %%
#Normalisation Functions/pca func
def get_pca_plot(df, title_string): #df needs to be without log
    #remove index from data 
    df.head()
    df = df.reset_index()
    df = df.drop("description", 1)
    df = df.drop("pep_num", 1)
    df = df.drop("uniprot_id", 1)
    df = df.drop("annotation", 1)
    #drop na
    df = df.dropna()
    #number_of_proteins_in_common = df.shape[0]

    #log2 transform
    df_log = np.log2(df)
    df_log = df_log.replace(-np.inf, np.nan)
    df_log = df_log.dropna()
    number_of_proteins_in_common = df_log.shape[0]
    
    #transpose table
    df_log_t = df_log.transpose()
    df_log_t.reset_index(inplace=True)
    df_log_t = df_log_t.rename(columns={"index":"channel_name"})

    #get X and y 
    X = df_log_t.drop('channel_name', 1)
    y = df_log_t['channel_name']

    #get PCA 2, fit transform, get df 
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(X)
    principalDf = pd.DataFrame(data = principalComponents, 
                               columns = ['principal component 1', 'principal component 2'])
    #add channel name
    finalDf = pd.concat([principalDf, df_log_t[['channel_name']]], axis = 1)

    #get PCA %
    pca_1_percent = round((pca.explained_variance_ratio_[0] * 100),2)
    pca_2_percent = round((pca.explained_variance_ratio_[1] * 100),2)

    #add condition
    finalDf["condition"] = finalDf["channel_name"]
    finalDf["condition"] = finalDf["condition"].str[:-2]
    
    #get range
    range_list = []

    range_list.append(finalDf["principal component 1"].max())
    range_list.append(finalDf["principal component 1"].min())

    range_list.append(finalDf["principal component 2"].max())
    range_list.append(finalDf["principal component 2"].min())

    max_number = max(range_list)+10
    min_number = min(range_list)-10

    fig = px.scatter(finalDf, 
                     x='principal component 1', 
                     y='principal component 2', 
                     hover_data=['channel_name'], 
                     color=finalDf['condition'],
                     color_discrete_sequence=["orange", "green", "blue", "purple", "gold", "lime", "dodgerblue", "rosybrown", "black"],  #, , , "lightblue", , "black", , "peru", , , "indigo", "teal", , "cyan", "fuchsia"
                     labels={
                         "principal component 1": "PC1 ({}%)".format(pca_1_percent),
                         "principal component 2": "PC2 ({}%)".format(pca_2_percent),
                         "condition": "Condition"},
                     title=title_string+" ("+str(number_of_proteins_in_common)+" Proteins)")

    fig.update_xaxes(dtick=10, range=[min_number, max_number])
    fig.update_yaxes(dtick=10, range=[min_number, max_number])

    fig.update_layout(plot_bgcolor='rgb(243, 243, 243)', 
                      height=500, width=600, showlegend=True)
    return(fig)

def pearson_corr_channels(df, file_name, folder_path, title_str):
    pearsoncorr = df.corr(method='pearson')

    fig, ax = plt.subplots(figsize=(15,15))  
    sns.heatmap(pearsoncorr, 
                xticklabels=pearsoncorr.columns,
                yticklabels=pearsoncorr.columns,
                cmap='RdBu_r',
                annot=True,
                linewidth=0.5,
                center=0.)
    ax.set_title(title_str+'- Pearson Correlation of Channels')
    plt.savefig(folder_path / (file_name+"_"+title_str+"_cor.png"))
    return("complete")

# %%
normalized_dict = {}

folder_path = output_path_obj / "03_normalisation_plots"
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
folder_path = output_path_obj / "04_normalised"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)

for file_name in list_of_file_names:
    df = normalized_dict[file_name]
    print(len(df))
    print("")
    df.to_csv(folder_path / (file_name + ".csv"))

# %% [markdown]
# **FUNCTIONS for cutoff analysis**

# %%
def get_ratio_condition_df(sub_df, cond_columns, ctrl_columns):
    for condition in cond_columns:
        for control in ctrl_columns:
            sub_df["ratio_"+condition+"/"+control] = sub_df[condition] / sub_df[control]
        
    ratio_df = sub_df.filter(regex='ratio_')
    return(ratio_df)

def get_cutoffs(ratio_dfs_dict, folder_path, all_cond_cutoff_tables, sub_cutoff_dict):
    
    folder_path = folder_path / ('plots')
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
    
    for condition, table in ratio_dfs_dict.items():     
        this_condition_full_tables = []
        this_condition_cutoff_tables = []
    
        columns_list = table.columns.tolist()
        table = table.reset_index()
    
        for column in columns_list:
            
            #get replicate
            df_rep = table[["uniprot_id", "annotation", column]]
        
            #f.write("Analyzing: " + column + "\n")
        
            #normalize ratios
            median_FP_ratio = df_rep.loc[df_rep['annotation'] == 'FP', column].median()
            #print(column)
            #print(median_FP_ratio)
            #df_rep[df_rep['annotation'] == 'FP'].to_csv("A.csv")
            #break
            
            df_rep[column] = df_rep[column].div(median_FP_ratio)
            df_rep["norm_"+column] = df_rep[column]
            
            df_rep[column] = np.log2(df_rep[column]) 
            log2_norm_ratio_column = column.replace("ratio_", "log2_norm_ratio_")
            df_rep.rename(columns={column: log2_norm_ratio_column}, inplace=True)
            column = column.replace("ratio_", "")

            df_rep = df_rep[df_rep[log2_norm_ratio_column].notna()]
        
            #get total number of annotated TP and FP
            total_TP = len(df_rep[df_rep["annotation"] == "TP"])
            total_FP = len(df_rep[df_rep["annotation"] == "FP"])
        
            #rank ratio in descending order
            df_rep = df_rep.sort_values(log2_norm_ratio_column, ascending=False)
        
            #clean up index
            df_rep = df_rep.reset_index(drop=True)
            #get index as list
            index_list = df_rep.index.values.tolist()
        
            #for each row calculate number of TP and FP in all the rows before
            for row_number in index_list:
                if row_number == 0:
                    df_rep.loc[df_rep.index[row_number], 'FP'] = 0
                    df_rep.loc[df_rep.index[row_number], 'TP'] = 0
            
                else:
                    results_dict = {}
                    subset_df = df_rep.loc[0:row_number-1]
                    for idx, name in enumerate(subset_df.annotation.value_counts().index.tolist()):
                        results_dict[name] = subset_df.annotation.value_counts()[idx]

                    if 'FP' in results_dict:
                        df_rep.loc[df_rep.index[row_number], 'FP'] = results_dict["FP"]
                    else:
                        df_rep.loc[df_rep.index[row_number], 'FP'] = 0
                    if 'TP' in results_dict:
                        df_rep.loc[df_rep.index[row_number], 'TP'] = results_dict["TP"]
                    else:
                        df_rep.loc[df_rep.index[row_number], 'TP'] = 0

            #calculate TPR and FPR
            df_rep["TPR"] = df_rep["TP"] / total_TP
            df_rep["FPR"] = df_rep["FP"] / total_FP
        
            #plot TPR and FPR
            plt.figure()
            plt.plot(df_rep["FPR"], df_rep["TPR"], color='red', label='ROC')
            plt.plot([0, 1], [0, 1], color='green', linestyle='--')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(column)
            plt.legend()
            plt.savefig(folder_path / (column.replace("/","-")+"_roc.png"))
        
            #calculate TPR-FPR, get cutoff value = maximum TPR-FPR value
            df_rep["TPR-FPR"] = df_rep["TPR"] - df_rep["FPR"]
            cutoff_value = df_rep.loc[df_rep["TPR-FPR"].idxmax()][log2_norm_ratio_column]
            #f.write("Cutoff value:" + str(cutoff_value) + "\n")
            sub_cutoff_dict[column] = cutoff_value
    
            #plot log2 ratio and TPR-FPR
            plt.figure()
            ax = sns.lineplot(df_rep[log2_norm_ratio_column], df_rep["TPR-FPR"])
    
            # Setting the X and Y Label
            xlabel_string = "Log2(" + column + ")"
            ax.set_xlabel(xlabel_string)
            ax.set_ylabel('TPR-FPR')
            ax.set_title("Cutoff:"+str(cutoff_value))
            ax.axvline(cutoff_value, linewidth=1, color="black", linestyle = "--")
            ax.figure.savefig(folder_path / (column.replace("/","-")+'_lineplot.png'))
    
            #Retain all proteins with log2 ratios higher than that of the determined cutoff
            #label proteins with pass cutoff or did not pass cutoff 
            #df_rep_cutoff["pass_cutoff"] = df_rep[log2_norm_ratio_column].where(df_rep[log2_norm_ratio_column] > cutoff_value, other='0')
            df_rep["pass_cutoff"] = np.where(df_rep[log2_norm_ratio_column] > cutoff_value, 1, 0)
            #df_rep_cutoff = df_rep[df_rep["FPR"] < 0.1] FDR 10 %
        
            this_condition_full_table = []
            this_condition_cutoff_table = []
        
            this_condition_full_tables.append(df_rep)
            this_condition_cutoff_tables.append(df_rep)
    
        all_cond_cutoff_tables[condition] = this_condition_cutoff_tables
    return(all_cond_cutoff_tables)

def add_pass_cutoff_analysis_to_df(ratio_and_signal_intensity_merged, df_all_list_columns_passcutoff):
    new_list = df_all_list_columns_passcutoff+["uniprot_id"]
    pass_cutoff_sub_df = ratio_and_signal_intensity_merged[new_list]
    pass_cutoff_sub_df.columns = pass_cutoff_sub_df.columns.str.split('/').str[0]
    pass_cutoff_sub_df.columns = pass_cutoff_sub_df.columns.str.replace('pass_cutoff_', "")
    sum_pass_cutoff = pass_cutoff_sub_df.transpose().reset_index().groupby("index").max().transpose() 
    sum_pass_cutoff.columns = sum_pass_cutoff.columns.str.rsplit("_", 1).str[0]
    sum_pass_cutoff = sum_pass_cutoff.reset_index()
    sum_pass_cutoff = sum_pass_cutoff.drop("index",1)
    new_df = sum_pass_cutoff.transpose().reset_index().groupby("index").sum().transpose()
    new_df = new_df.set_index("uniprot")
    new_df = new_df.add_suffix('_pass_cutoff_sum')
    col_list = new_df.columns.tolist()
    for col in col_list:
        new_df[col+"_min_2_cutoffs"] = new_df[col] >= 2
    new_df["pass_cutoff_result"] = (new_df == True).any(axis=1)
    new_df = new_df.reset_index()
    new_df = new_df.rename(columns={"uniprot":"uniprot_id"})
    ratio_and_signal_intensity_merged = ratio_and_signal_intensity_merged.reset_index()
    with_cutoff_df = ratio_and_signal_intensity_merged.merge(new_df, on="uniprot_id") 
    return(with_cutoff_df)

# %% [markdown]
# **Cutoff analysis of all files**

# %%
folder_path = output_path_obj / "05_cutoff_analysis"
if not os.path.exists(folder_path):
    os.mkdir(folder_path)
    
cutoff_dict = {}
for file_name in list_of_file_names:
    print(file_name)
    
    #create output folder
    folder_name = file_name.replace("processed_census-out_", "")
    folder_path = output_path_obj / "05_cutoff_analysis" / folder_name
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
    
    #reduced_master_list = master_list.drop(master_list.columns[[0,-10,-9,-8,-7,-6]], axis = 1)
    
    ratio_and_signal_intensity = ratio_and_signal_intensity.rename(columns={"uniprot_id": "Entry"})
    ratio_and_signal_intensity.drop('description', axis=1, inplace=True)
    master_list = master_list.rename(columns={"uniprot_id": "alias_uniprot_id"})

    merged_with_metadata = pd.merge(ratio_and_signal_intensity, master_list, on=["Entry", "annotation"], how="left")
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

# %% [markdown]
# **Analysis of Proteins that pass cutoff Filter (i.e. pass_cutoff_result=TRUE)**

# %%
def get_pca_plot_after(df, title_string): #df needs to be without log
    #remove index from data 
    df.head()
    df = df.reset_index()
    df = df.drop("Signal peptide", 1)
    df = df.drop("uniprot_id", 1)
    df = df.drop("annotation", 1)
    #drop na
    df = df.dropna()
    #number_of_proteins_in_common = df.shape[0]

    #log2 transform
    df_log = np.log2(df)
    df_log = df_log.replace(-np.inf, np.nan)
    df_log = df_log.dropna()
    number_of_proteins_in_common = df_log.shape[0]
    #transpose table
    df_log_t = df_log.transpose()
    df_log_t.reset_index(inplace=True)
    df_log_t = df_log_t.rename(columns={"index":"channel_name"})

    #get X and y 
    X = df_log_t.drop('channel_name', 1)
    y = df_log_t['channel_name']

    #get PCA 2, fit transform, get df 
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(X)
    principalDf = pd.DataFrame(data = principalComponents, 
                               columns = ['principal component 1', 'principal component 2'])
    #add channel name
    finalDf = pd.concat([principalDf, df_log_t[['channel_name']]], axis = 1)

    #get PCA %
    pca_1_percent = round((pca.explained_variance_ratio_[0] * 100),2)
    pca_2_percent = round((pca.explained_variance_ratio_[1] * 100),2)

    #add condition
    finalDf["condition"] = finalDf["channel_name"]
    finalDf["condition"] = finalDf["condition"].str[:-2]
    
    #get range
    range_list = []

    range_list.append(finalDf["principal component 1"].max())
    range_list.append(finalDf["principal component 1"].min())

    range_list.append(finalDf["principal component 2"].max())
    range_list.append(finalDf["principal component 2"].min())

    max_number = max(range_list)+10
    min_number = min(range_list)-10

    fig = px.scatter(finalDf, 
                     x='principal component 1', 
                     y='principal component 2', 
                     hover_data=['channel_name'], 
                     color=finalDf['condition'],
                     color_discrete_sequence=["orange", "green", "blue", "purple", "gold", "lime", "dodgerblue", "rosybrown", "black"],  #, , , "lightblue", , "black", , "peru", , , "indigo", "teal", , "cyan", "fuchsia"
                     labels={
                         "principal component 1": "PC1 ({}%)".format(pca_1_percent),
                         "principal component 2": "PC2 ({}%)".format(pca_2_percent),
                         "condition": "Condition"},
                     title=title_string+" ("+str(number_of_proteins_in_common)+" Proteins)")

    fig.update_xaxes(dtick=10, range=[min_number, max_number])
    fig.update_yaxes(dtick=10, range=[min_number, max_number])

    fig.update_layout(plot_bgcolor='rgb(243, 243, 243)', 
                      height=500, width=600, showlegend=True)
    return(fig)

# %%
input_folder_path = output_path_obj / "05_cutoff_analysis" 
list_of_dir_paths = list(input_folder_path.iterdir())

list_of_result_file_paths = []

for directory in list_of_dir_paths:
    list_in_dir = list(directory.iterdir())
    for file_path in list_in_dir: 
        if "results" in file_path.stem:
            list_of_result_file_paths.append(file_path)

# %%
folder_path = output_path_obj / "06_results" 
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
        
        folder_path = output_path_obj / "06_results" / file_path.stem
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
#FUNCTIONS
def get_p_value(row, cond_1, cond_2):
    ttest_result = stat.ttest_ind(row[cond_1],  row[cond_2])
    return(ttest_result[1])

def get_expr(row):
    if (row['log2_FC'] > 0.58) & (row['-log10_pval'] > 1.3):
        return("Up")
    elif (row['log2_FC'] < -0.58) & (row['-log10_pval'] > 1.3):
        return("Down")
    else:
        return("Stable")

# %%
folder_path = output_path_obj / "06_results" 
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
        

# %%



