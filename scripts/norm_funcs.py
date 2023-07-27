import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import plotly.express as px

def get_condition_df(df, condition_name):
    sub_df = df.filter(regex=condition_name)
    return(sub_df)

def normalization_cre_groups(df, conditions_list, treatment_labelling, control_labelling):
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
    #logging.info("Normalization factors")
    #logging.info(norm_factor)

    normalized_df = df[norm_factor.index] * norm_factor
    list_of_ctrl_cols = [column for column in df.columns if column not in norm_factor.index]
    untouched_ctrl_columns = df[list_of_ctrl_cols]
    
    norm_df = normalized_df.reset_index().merge(untouched_ctrl_columns.reset_index(), on=["uniprot_id", "description", "pep_num", "annotation"])
    norm_df = norm_df.set_index(["uniprot_id", "description", "pep_num", "annotation"])
    norm_df = norm_df[df.columns] # get the same order of columns 
    return norm_df

def normalization_allcre_channels(df, treatment_labelling):
    norm_factor_series = []
    
    cond_columns = df.filter(like=treatment_labelling).columns.tolist()
    sub_df = df[cond_columns]
        
    norm_factors = sub_df.median().median() / sub_df.median()
    norm_factor_series.append(norm_factors)

    norm_factor = pd.concat(norm_factor_series, axis=1).sum(1)
    print(norm_factor)
    #logging.info("Normalization factors")
    #logging.info(norm_factor)

    normalized_df = df[norm_factor.index] * norm_factor
    list_of_ctrl_cols = [column for column in df.columns if column not in norm_factor.index]
    untouched_ctrl_columns = df[list_of_ctrl_cols]
    
    norm_df = normalized_df.reset_index().merge(untouched_ctrl_columns.reset_index(), on=["uniprot_id", "description", "pep_num", "annotation"])
    norm_df = norm_df.set_index(["uniprot_id", "description", "pep_num", "annotation"])
    norm_df = norm_df[df.columns] # get the same order of columns 
    return norm_df

def normalization_all_channels(df):
    norm_factor_series = []

    norm_factors = df.median().median() / df.median()
    norm_factor_series.append(norm_factors)

    norm_factor = pd.concat(norm_factor_series, axis=1).sum(1)
    #logging.info("Normalization factors")
    print(norm_factor)
    #logging.info(norm_factor)

    normalized_df = df[norm_factor.index] * norm_factor
    #list_of_ctrl_cols = [column for column in df.columns if column not in norm_factor.index]
    #untouched_ctrl_columns = df[list_of_ctrl_cols]
    
    #norm_df = normalized_df.reset_index().merge(untouched_ctrl_columns.reset_index(), on=["uniprot_id", "description", "pep_num", "annotation"])
    #norm_df = norm_df.set_index(["uniprot_id", "description", "pep_num", "annotation"])
    normalized_df = normalized_df[df.columns] # get the same order of columns 
    return normalized_df

# Median normalization is a global normalization method correcting for differential sample amounts in a robust manner. 
# It centers the sample data to its corresponding median. 
# First, median protein abundance ratio is calculated for each sample. 
# Next, each individual protein abundance in a channel is divided by its corresponding median

def normalization_all_channels_median(df):
    norm_factor_series = []

    norm_factors = df.median() #.median() / df.median()
    norm_factor_series.append(norm_factors)

    norm_factor = pd.concat(norm_factor_series, axis=1).sum(1)
    #logging.info("Normalization factors")
    print(norm_factor)
    #logging.info(norm_factor)

    normalized_df = df[norm_factor.index] / norm_factor
    #list_of_ctrl_cols = [column for column in df.columns if column not in norm_factor.index]
    #untouched_ctrl_columns = df[list_of_ctrl_cols]
    
    #norm_df = normalized_df.reset_index().merge(untouched_ctrl_columns.reset_index(), on=["uniprot_id", "description", "pep_num", "annotation"])
    #norm_df = norm_df.set_index(["uniprot_id", "description", "pep_num", "annotation"])
    normalized_df = normalized_df[df.columns] # get the same order of columns 
    return normalized_df



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

    fig.update_traces(marker=dict(size=12),
                      selector=dict(mode='markers'))

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