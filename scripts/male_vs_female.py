# %%
import pandas as pd 
from pathlib import Path
import seaborn as sns

# %%
filt_male_path = Path("/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/052423_results/min_1_pep_per_protein/02_filtering_per_cond_per_file/20230512_EV2-28A/passed_filter_processed_census-out_20230512_EV2-28A.csv")
filt_female_path = Path("/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/052423_results/min_1_pep_per_protein/02_filtering_per_cond_per_file/20230513_EV2-28B/passed_filter_processed_census-out_20230513_EV2-28B.csv")

filt_norm_male_path = Path("/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/052423_results/min_1_pep_per_protein/04_normalized/20230512_EV2-28A.csv")
filt_norm_female_path = Path("/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/052423_results/min_1_pep_per_protein/04_normalized/20230513_EV2-28B.csv")

filt_norm_roc_male_path = "/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/052423_dontuse/min_1_pep_per_protein/05_results/01_tissue/20230512_EV2-28A/final_protein_table_20230512_EV2-28A.csv"
filt_norm_roc_female_path = "/Users/nropek/Dropbox (Dropbox @RU)/TurboID manuscript/Mass-spectrometry datasets/03_results/03_downstream_analysis/052423_dontuse/min_1_pep_per_protein/05_results/01_tissue/20230513_EV2-28B/final_protein_table_20230513_EV2-28B.csv"

# %%
def channel_ratio_suffix(file_path, suffix, only_cre_plus=False):
    df = pd.read_csv(file_path) 
    df = df.set_index(["uniprot_id", "description", "pep_num", "annotation"])
    if only_cre_plus==True:
        df = df.filter(like="(+)", axis=1)

    if df.shape[1] > 16:
        df = df.iloc[:, 20:36]
        print(df.columns)
    #get channel ratio
    df_cols = df.columns.tolist()
    df["sum_per_row"] = df[df_cols].sum(axis='columns')
    df = df[df_cols].div(df.sum_per_row, axis=0)
    #add suffix 
    df = df.add_suffix("_" + suffix)
    return df

def get_pca_plot(df, title_string, color_list): #df needs to be without log
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
    finalDf["gender"] = finalDf["condition"].str.rsplit("_", 1).str[1]
    finalDf["condition"] = finalDf["condition"].str.rsplit("_", 2).str[0] #[:-2]
    
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
                     symbol=finalDf['gender'],
                     color_discrete_sequence=color_list, #["lightgreen", "green", "lightblue", "blue", "gold", "orange", "dodgerblue", "rosybrown", "black"],  #, , , "lightblue", , "black", , "peru", , , "indigo", "teal", , "cyan", "fuchsia"
                     labels={
                         "principal component 1": "PC1 ({}%)".format(pca_1_percent),
                         "principal component 2": "PC2 ({}%)".format(pca_2_percent),
                         "condition": "Condition"},
                     title=title_string+" ("+str(number_of_proteins_in_common)+" Proteins)")

    fig.update_traces(marker=dict(size=12),
                      selector=dict(mode='markers'))

    fig.update_xaxes(dtick=10, range=[min_number, max_number])
    fig.update_yaxes(dtick=10, range=[min_number, max_number])

    fig.update_layout(plot_bgcolor='rgb(243, 243, 243)', 
                      height=500, width=600, showlegend=True)
    return(fig)

def get_heatmap(df, treatment_labelling, file_name):
    subset_for_heatmap = df.filter(like=treatment_labelling, axis=1)
    cols_list = subset_for_heatmap.columns.tolist()
    subset_for_heatmap = subset_for_heatmap.reset_index().drop(["pep_num", "annotation", "uniprot_id"], axis=1) #.set_index("description")
    subset_for_heatmap.description = subset_for_heatmap.description.str.split("GN=").str[1]
    subset_for_heatmap.replace([np.inf, -np.inf], np.nan, inplace=True)
    subset_for_heatmap.dropna(inplace=True)
    subset_for_heatmap = subset_for_heatmap.set_index("description")

    g = sns.clustermap(subset_for_heatmap, z_score=0, cmap=sns.diverging_palette(15, 145, s=60, as_cmap=True), center=0, figsize=(12,int(len(subset_for_heatmap) / 5 )), yticklabels=True)
    #g.ax_row_dendrogram.set_visible(False)
    #g.ax_col_dendrogram.set_visible(False)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 10)
    g.ax_heatmap.set_title(file_name) 
    g.ax_cbar.set_ylabel("z-score",size=15)
    g.ax_cbar.set_position((0.1, .2, .03, .4))

    #g.savefig("heatmap_malevsfemale.pdf", dpi=400)
    #subset_for_heatmap_merged.to_csv(folder_path / "heatmap_data.csv") 
    return "done"

# %% 
color_list = ["lightgreen", "green", "lightblue", "blue", "gold", "orange"]
filt_male_df = channel_ratio_suffix(filt_male_path, "male")
filt_female_df = channel_ratio_suffix(filt_female_path, "female")
filt_merged = filt_male_df.join(filt_female_df, how="outer")
get_pca_plot(filt_merged, "PCA", color_list)
# %%
filt_male_df = channel_ratio_suffix(filt_norm_male_path, "male")
filt_female_df = channel_ratio_suffix(filt_norm_female_path, "female")
filt_merged = filt_male_df.join(filt_female_df, how="outer")
get_pca_plot(filt_merged, "PCA", color_list)
filt_merged.to_csv("filt_merged_after_roc.csv")

# %% 
filt_male_df = channel_ratio_suffix(filt_norm_roc_male_path, "male")
filt_female_df = channel_ratio_suffix(filt_norm_roc_female_path, "female")
filt_merged = filt_male_df.join(filt_female_df, how="outer")
get_pca_plot(filt_merged, "PCA", color_list)
get_heatmap(filt_merged, '(+)', "male vs female")

# %%
color_list = ["green", "blue", "orange"]
filt_male_df = channel_ratio_suffix(filt_male_path, "male", True)
filt_female_df = channel_ratio_suffix(filt_female_path, "female", True)
filt_merged = filt_male_df.join(filt_female_df, how="outer")
get_pca_plot(filt_merged, "PCA", color_list)
# %%
filt_male_df = channel_ratio_suffix(filt_norm_male_path, "male", True)
filt_female_df = channel_ratio_suffix(filt_norm_female_path, "female", True)
filt_merged = filt_male_df.join(filt_female_df, how="outer")
get_pca_plot(filt_merged, "PCA", color_list)
# %%
filt_male_df = channel_ratio_suffix(filt_norm_roc_male_path, "male", True)
filt_female_df = channel_ratio_suffix(filt_norm_roc_female_path, "female", True)
filt_merged = filt_male_df.join(filt_female_df, how="outer")
get_pca_plot(filt_merged, "PCA", color_list)
# %%
