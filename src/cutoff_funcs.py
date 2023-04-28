import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import plotly.express as px
import scipy.stats as stat
from functools import reduce
from src.filter_funcs import get_condition_df


def get_detailed_protein_annotation(df, volcano_df, fasta_table): 
    df = df.reset_index().set_index("uniprot_id")
    df = df.drop(["annotation", "description"], axis = 1)

    annotated_df = pd.concat([df, volcano_df, fasta_table], axis=1)
    return annotated_df


def get_ratio_condition_df(sub_df, cond_columns, ctrl_columns):
    for condition in cond_columns:
        for control in ctrl_columns:
            sub_df["ratio_"+condition+"/"+control] = sub_df[condition] / sub_df[control]    
    ratio_df = sub_df.filter(regex='ratio_')
    return ratio_df


def get_cutoffs(ratio_dfs_dict, folder_path, all_cond_cutoff_tables, sub_cutoff_dict):
    folder_path = folder_path / ('cutoff_roc_plots')
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)
    for condition, table in ratio_dfs_dict.items():     
        this_condition_full_tables = []
        this_condition_cutoff_tables = []
        columns_list = table.columns.tolist()
        table = table.reset_index()
        for column in columns_list:
            # get replicate
            df_rep = table[["uniprot_id", "annotation", column]]
            # normalize ratios
            median_FP_ratio = df_rep.loc[df_rep['annotation'] == 'FP', column].median()
            df_rep[column] = df_rep[column].div(median_FP_ratio)
            df_rep["norm_"+column] = df_rep[column]
            df_rep[column] = np.log2(df_rep[column]) 
            log2_norm_ratio_column = column.replace("ratio_", "log2_norm_ratio_")
            df_rep.rename(columns={column: log2_norm_ratio_column}, inplace=True)
            column = column.replace("ratio_", "")
            df_rep = df_rep[df_rep[log2_norm_ratio_column].notna()]
            # get total number of annotated TP and FP
            total_TP = len(df_rep[df_rep["annotation"] == "TP"])
            total_FP = len(df_rep[df_rep["annotation"] == "FP"])
            # rank ratio in descending order
            df_rep = df_rep.sort_values(log2_norm_ratio_column, ascending=False)
            # clean up index
            df_rep = df_rep.reset_index(drop=True)
            # get index as list
            index_list = df_rep.index.values.tolist()
            # for each row calculate number of TP and FP in all the rows before
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

            # calculate TPR and FPR
            df_rep["TPR"] = df_rep["TP"] / total_TP
            df_rep["FPR"] = df_rep["FP"] / total_FP
            # plot TPR and FPR
            plt.figure()
            plt.plot(df_rep["FPR"], df_rep["TPR"], color='red', label='ROC')
            plt.plot([0, 1], [0, 1], color='green', linestyle='--')
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(column)
            plt.legend()
            plt.savefig(folder_path / (column.replace("/","-")+"_roc.png"))

            # calculate TPR-FPR, get cutoff value = maximum TPR-FPR value
            df_rep["TPR-FPR"] = df_rep["TPR"] - df_rep["FPR"]
            cutoff_value = df_rep.loc[df_rep["TPR-FPR"].idxmax()][log2_norm_ratio_column]
            sub_cutoff_dict[column] = cutoff_value

            # plot log2 ratio and TPR-FPR
            plt.figure()
            ax = sns.lineplot(df_rep[log2_norm_ratio_column], df_rep["TPR-FPR"])

            # Setting the X and Y Label
            xlabel_string = "Log2(" + column + ")"
            ax.set_xlabel(xlabel_string)
            ax.set_ylabel('TPR-FPR')
            ax.set_title("Cutoff:"+str(cutoff_value))
            ax.axvline(cutoff_value, linewidth=1, color="black", linestyle = "--")
            ax.figure.savefig(folder_path / (column.replace("/","-")+'_lineplot.png'))
    
            # Retain all proteins with log2 ratios higher than that of the determined cutoff
            # label proteins with pass cutoff or did not pass cutoff 
            # df_rep_cutoff["pass_cutoff"] = df_rep[log2_norm_ratio_column].where(df_rep[log2_norm_ratio_column] > cutoff_value, other='0')
            df_rep["pass_cutoff"] = np.where(df_rep[log2_norm_ratio_column] > cutoff_value, 1, 0)
            # df_rep_cutoff = df_rep[df_rep["FPR"] < 0.1] FDR 10 %
        
            this_condition_full_table = []
            this_condition_cutoff_table = []
            this_condition_full_tables.append(df_rep)
            this_condition_cutoff_tables.append(df_rep)
        all_cond_cutoff_tables[condition] = this_condition_cutoff_tables
    return all_cond_cutoff_tables 


def add_pass_cutoff_analysis_to_df(ratio_and_signal_intensity_merged, df_all_list_columns_passcutoff):
    new_list = df_all_list_columns_passcutoff+["uniprot_id"]
    pass_cutoff_sub_df = ratio_and_signal_intensity_merged[new_list]
    pass_cutoff_sub_df.columns = pass_cutoff_sub_df.columns.str.split('/').str[0]
    pass_cutoff_sub_df.columns = pass_cutoff_sub_df.columns.str.replace('pass_cutoff_', "")
    sum_pass_cutoff = pass_cutoff_sub_df.transpose().reset_index().groupby("index").max().transpose() 
    sum_pass_cutoff.columns = sum_pass_cutoff.columns.str.rsplit("_", 1).str[0]
    sum_pass_cutoff = sum_pass_cutoff.reset_index()
    sum_pass_cutoff = sum_pass_cutoff.drop("index", 1)
    new_df = sum_pass_cutoff.transpose().reset_index().groupby("index").sum().transpose()
    new_df = new_df.set_index("uniprot")
    new_df = new_df.add_suffix('_pass_cutoff_sum')
    col_list = new_df.columns.tolist()
    for col in col_list:
        new_df[col+"_min_2_cutoffs"] = new_df[col] >= 2
    
    min_2_cutoffs_list = new_df.filter(like="_min_2_cutoffs").columns.tolist()
    new_df["pass_cutoff_result"] = (new_df[min_2_cutoffs_list] == True).any(axis=1)
    new_df = new_df.reset_index()
    new_df = new_df.rename(columns={"uniprot":"uniprot_id"})
    ratio_and_signal_intensity_merged = ratio_and_signal_intensity_merged.reset_index()
    with_cutoff_df = ratio_and_signal_intensity_merged.merge(new_df, on="uniprot_id") 
    return with_cutoff_df 


def get_pca_plot_after(df, title_string): #df needs to be without log
    # remove index from data 
    df.head()
    df = df.reset_index()
    df = df.drop("Signal peptide", 1)
    df = df.drop("uniprot_id", 1)
    df = df.drop("annotation", 1)
    df = df.dropna()
    # number_of_proteins_in_common = df.shape[0]

    # log2 transform
    df_log = np.log2(df)
    df_log = df_log.replace(-np.inf, np.nan)
    df_log = df_log.dropna()
    number_of_proteins_in_common = df_log.shape[0]
    # transpose table
    df_log_t = df_log.transpose()
    df_log_t.reset_index(inplace=True)
    df_log_t = df_log_t.rename(columns={"index":"channel_name"})

    # get X and y 
    X = df_log_t.drop('channel_name', 1)
    y = df_log_t['channel_name']

    # get PCA 2, fit transform, get df 
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(X)
    principalDf = pd.DataFrame(data = principalComponents, 
                               columns = ['principal component 1', 'principal component 2'])
    # add channel name
    finalDf = pd.concat([principalDf, df_log_t[['channel_name']]], axis = 1)

    # get PCA %
    pca_1_percent = round((pca.explained_variance_ratio_[0] * 100),2)
    pca_2_percent = round((pca.explained_variance_ratio_[1] * 100),2)

    # add condition
    finalDf["condition"] = finalDf["channel_name"]
    finalDf["condition"] = finalDf["condition"].str[:-2]
    
    # get range
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
    return fig


def get_p_value(row, cond_1, cond_2):
    ttest_result = stat.ttest_ind(row[cond_1],  row[cond_2])
    return ttest_result[1]


def get_expr(row):
    if (row['log2_FC'] > 0.58) & (row['-log10_pval'] > 1.3):
        return "Up"
    elif (row['log2_FC'] < -0.58) & (row['-log10_pval'] > 1.3):
        return "Down"
    else:
        return "Stable"


def get_volcano_plot(conditions_list, control_labelling, treatment_labelling, df, file_name, folder_path):
    volcano_plot_list = []
    volcano_df_list = []
    if len(conditions_list) == 2:
        list_columns_trt = df.filter(like=treatment_labelling).columns.tolist()
        trt_col_df = df[list_columns_trt]
        list_cond_1 = trt_col_df.filter(like=conditions_list[0]).columns.tolist()
        list_cond_2 = trt_col_df.filter(like=conditions_list[1]).columns.tolist()

        for condition in conditions_list:
            list_columns_cond = trt_col_df.filter(like=condition).columns.tolist()
            cond_df = trt_col_df[list_columns_cond]
            trt_col_df["mean_"+condition+"_"+treatment_labelling] = cond_df.mean(axis=1)
            
        trt_col_df["FC"] = trt_col_df["mean_"+conditions_list[0]+"_"+treatment_labelling] / trt_col_df["mean_"+conditions_list[1]+"_"+treatment_labelling]
        x_axis_name = "log2(" + conditions_list[0] +"/" + conditions_list[1] +")"
    
        idx_cond_1 = trt_col_df.columns.get_indexer(list_cond_1)
        idx_cond_2 = trt_col_df.columns.get_indexer(list_cond_2)
        trt_col_df["p_value"] = trt_col_df.apply(get_p_value, axis=1, args=(idx_cond_1, idx_cond_2), )
        
        trt_col_df["log2_FC"] = np.log2(trt_col_df["FC"])
        volcano_df = trt_col_df[["p_value", "log2_FC"]]
            
        volcano_df = volcano_df.dropna() 
        volcano_df["-log10_pval"] = -1*np.log10(volcano_df["p_value"])
        volcano_df["Regulation"] = volcano_df.apply(get_expr, axis=1)
            
        volcano_df["Regulation"] = volcano_df["Regulation"].astype('category')
        
        volcano_df = volcano_df.reset_index()
        volcano_df = volcano_df.set_index(["annotation"])
        volcano_df["name"] = volcano_df['description'].str.split(" ").str[0]

        title_name = file_name + " - " + conditions_list[0] +" vs. " + conditions_list[1] + " ("+ str(len(volcano_df)) + " Proteins)"
        my_order = ["Down", "Stable", "Up"]
        volcano_df['Regulation'].cat.reorder_categories(my_order, inplace= True)
        colors_volcano = ["red", 'grey', 'green']

        fig = px.scatter(volcano_df, 
                         x='log2_FC', 
                         y='-log10_pval',
                         color='Regulation',
                         color_discrete_sequence=colors_volcano,
                         hover_data=['name', 'uniprot_id'],
                         title = title_name,
                         labels = {"log2_FC": x_axis_name},
                         template = "simple_white",
                         category_orders={'Regulation': np.sort(volcano_df['Regulation'].unique())})
        fig.add_vline(x=0.58, line_width=2, line_dash="dash", line_color="grey")
        fig.add_vline(x=-0.58, line_width=2, line_dash="dash", line_color="grey")
        fig.add_hline(y=1.3, line_width=2, line_dash="dash", line_color="grey")
        fig.update_layout(legend=dict(title=""), title_x=0.5)
        volcano_plot_list.append(fig)
        volcano_df = volcano_df.reset_index().set_index("uniprot_id")
        volcano_df = volcano_df.drop(["annotation", "description", "pep_num", "name"], axis=1)
        volcano_df = volcano_df.add_suffix("_"+title_name)
        volcano_df_list.append(volcano_df)
        
    elif len(conditions_list) == 3: 
        list_columns_trt = df.filter(like=treatment_labelling).columns.tolist()
        trt_col_df = df[list_columns_trt]
    
        list_cond_1 = trt_col_df.filter(like=conditions_list[0]).columns.tolist()
        list_cond_2 = trt_col_df.filter(like=conditions_list[1]).columns.tolist()
        list_cond_3 = trt_col_df.filter(like=conditions_list[2]).columns.tolist()
            
        idx_cond_1 = trt_col_df.columns.get_indexer(list_cond_1)
        idx_cond_2 = trt_col_df.columns.get_indexer(list_cond_2)
        idx_cond_3 = trt_col_df.columns.get_indexer(list_cond_3)
            
        trt_col_df["p_value_1"] = trt_col_df.apply(get_p_value, axis=1, args=(idx_cond_1, idx_cond_2), )
        trt_col_df["p_value_2"] = trt_col_df.apply(get_p_value, axis=1, args=(idx_cond_1, idx_cond_3), )
        trt_col_df["p_value_3"] = trt_col_df.apply(get_p_value, axis=1, args=(idx_cond_2, idx_cond_3), )
            
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
                title_name = file_name.split("processed_census-out_")[1] + " - " + conditions_list[1] +" vs. " + conditions_list[0] + " ("+ str(len(volcano_df)) + " Proteins)"
            elif "p_value_2" in volcano_df.columns:
                x_axis_name = "log2(" + conditions_list[2] +"/" + conditions_list[0] +")"
                title_name = file_name.split("processed_census-out_")[1] + " - " + conditions_list[2] +" vs. " + conditions_list[0] + " ("+ str(len(volcano_df)) + " Proteins)"
            elif "p_value_3" in volcano_df.columns:
                x_axis_name = "log2(" + conditions_list[1] +"/" + conditions_list[2] +")"
                title_name = file_name.split("processed_census-out_")[1] + " - " + conditions_list[1] +" vs. " + conditions_list[2] + " ("+ str(len(volcano_df)) + " Proteins)"
            
            volcano_df.columns = ["p_value", "log2_FC"]
            volcano_df = volcano_df.dropna() 
            volcano_df["-log10_pval"] = -1*np.log10(volcano_df["p_value"])
            volcano_df["Regulation"] = volcano_df.apply(get_expr, axis=1)
            volcano_df["Regulation"] = volcano_df["Regulation"].astype('category')
            my_order = ["Down", "Stable", "Up"]
            #volcano_df['Regulation'].cat.reorder_categories(my_order, inplace= True)
            
            volcano_df = volcano_df.reset_index()
            volcano_df = volcano_df.set_index(["annotation"])
            volcano_df["name"] = volcano_df['description'].str.split(" ").str[0]
                
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
            volcano_plot_list.append(fig)
            volcano_df = volcano_df.reset_index().set_index("uniprot_id")
            volcano_df = volcano_df.drop(["annotation", "description", "pep_num", "name"], axis=1)
            volcano_df = volcano_df.add_suffix("_"+title_name)
            volcano_df_list.append(volcano_df)

    with open(folder_path / 'volcano_plots.html' , 'w') as f:
        for fig in volcano_plot_list:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))

    volcano_df = pd.concat(volcano_df_list, axis=1)
            
    return volcano_df

def get_volcano_plot_treatment_vs_control(conditions_list, control_labelling, treatment_labelling, df, file_name, folder_path):
    volcano_plot_list = []
    volcano_df_list = []

    for condition in conditions_list:
        cond_cols_df = df.filter(like=condition)
        list_cond_1 = cond_cols_df.filter(like=treatment_labelling).columns.tolist()
        list_cond_2 = cond_cols_df.filter(like=control_labelling).columns.tolist()
        if len(list_cond_1) <= 1 or len(list_cond_2) <= 1:
            print(f"Condition {condition} does not have enough replicates to be shown in volcano plot!")
            print(list_cond_1)
            print(list_cond_2)
            continue

        for labelling in [control_labelling, treatment_labelling]:
            list_columns_labelling = cond_cols_df.filter(like=labelling).columns.tolist()
            labelling_df = cond_cols_df[list_columns_labelling]
            cond_cols_df["mean_"+condition+"_"+labelling] = labelling_df.mean(axis=1)
            
        cond_cols_df["FC"] = cond_cols_df["mean_"+condition+"_"+treatment_labelling] / cond_cols_df["mean_"+condition+"_"+control_labelling]
        x_axis_name = "log2(" + condition + "_" + treatment_labelling + "/" + condition + "_" + control_labelling + ")"

        idx_cond_1 = cond_cols_df.columns.get_indexer(list_cond_1)
        idx_cond_2 = cond_cols_df.columns.get_indexer(list_cond_2)
        cond_cols_df["p_value"] = cond_cols_df.apply(get_p_value, axis=1, args=(idx_cond_1, idx_cond_2), )
        
        cond_cols_df["log2_FC"] = np.log2(cond_cols_df["FC"])
        volcano_df = cond_cols_df[["p_value", "log2_FC"]]
            
        volcano_df = volcano_df.dropna() 
        volcano_df["-log10_pval"] = -1*np.log10(volcano_df["p_value"])
        volcano_df["Regulation"] = volcano_df.apply(get_expr, axis=1)
            
        volcano_df["Regulation"] = volcano_df["Regulation"].astype('category')
        
        volcano_df = volcano_df.reset_index()
        volcano_df = volcano_df.set_index(["annotation"])
        volcano_df["name"] = volcano_df['description'].str.split(" ").str[0]

        title_name = file_name + " - " + condition + "_" +treatment_labelling + " vs. " + condition + "_" + control_labelling + " (" + str(len(volcano_df)) + " Proteins)"
        #my_order = ["Down", "Stable", "Up"]
        my_order = ["Stable", "Up"]
        volcano_df['Regulation'].cat.reorder_categories(my_order, inplace= True)
        colors_volcano = ['grey', 'green']
        #colors_volcano = ["red", 'grey', 'green']
    
        fig = px.scatter(volcano_df, 
                         x='log2_FC', 
                         y='-log10_pval',
                         color='Regulation',
                         color_discrete_sequence=colors_volcano,
                         hover_data=['name', 'uniprot_id'],
                         title = title_name,
                         labels = {"log2_FC": x_axis_name},
                         template = "simple_white",
                         category_orders={'Regulation': np.sort(volcano_df['Regulation'].unique())})
        fig.add_vline(x=0.58, line_width=2, line_dash="dash", line_color="grey")
        fig.add_vline(x=-0.58, line_width=2, line_dash="dash", line_color="grey")
        fig.add_hline(y=1.3, line_width=2, line_dash="dash", line_color="grey")
        fig.update_layout(legend=dict(title=""), title_x=0.5)
        volcano_plot_list.append(fig)
        volcano_df = volcano_df.reset_index().set_index("uniprot_id")
        volcano_df = volcano_df.drop(["annotation", "description", "pep_num", "name"], axis=1)
        volcano_df = volcano_df.add_suffix("_"+title_name)
        volcano_df_list.append(volcano_df)
        
    with open(folder_path / 'volcano_plots.html' , 'w') as f:
        for fig in volcano_plot_list:
            f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))

    volcano_df = pd.concat(volcano_df_list, axis=1)
            
    return volcano_df


def get_ratios_and_cutoffs(df, 
                           conditions_list, 
                           control_labelling, 
                           treatment_labelling, 
                           tissue_file_folder_path, 
                           cutoff_dict, 
                           file_name, 
                           fasta_table):
    # get all ratios
    ratio_dfs_dict = {}

    for condition in conditions_list:
        sub_df = get_condition_df(df, condition)
        sub_df = sub_df.dropna() 
        cond_columns = sub_df.filter(like=treatment_labelling).columns.tolist()
        ctrl_columns = sub_df.filter(like=control_labelling).columns.tolist()
        ratio_df = get_ratio_condition_df(sub_df, cond_columns, ctrl_columns)
        ratio_dfs_dict[condition] = ratio_df
        
    ratio_tables = list(ratio_dfs_dict.values())
    ratio_tables = [table.reset_index() for table in ratio_tables]
    col_list = ['uniprot_id', 'description', 'pep_num', 'annotation']
    ratio_tables_merged = reduce(lambda df1, df2: pd.merge(df1, df2, on = col_list, how="outer"), ratio_tables)
    
    # merge with normalized SI
    ratio_and_signal_intensity = ratio_tables_merged.merge(df, on=['uniprot_id', 'description', 'pep_num'], how='left')
        
    all_cond_cutoff_tables = {}
    sub_cutoff_dict = {}
    all_cond_cutoff_tables = get_cutoffs(ratio_dfs_dict, tissue_file_folder_path, all_cond_cutoff_tables, sub_cutoff_dict)
    cutoff_dict[file_name] = sub_cutoff_dict
    
    # For each replicate get number of proteins that pass cutoff
    for condition, table_list in all_cond_cutoff_tables.items():
        for table in table_list:
            col_name = table.filter(like='log2_norm_ratio_').columns.tolist()
            col_name = col_name[0]
            col_name = col_name.replace("log2_norm_ratio_", "")
            table.drop('TP', axis=1, inplace=True)
            table.drop('FP', axis=1, inplace=True)
            table.drop('TPR', axis=1, inplace=True)
            table.drop('FPR', axis=1, inplace=True)
            table.rename(columns={"pass_cutoff": "pass_cutoff_"+col_name,
                                  "TPR-FPR": "TPR-FPR_"+col_name}, inplace=True)
            
    # Merge replicates of the same condition into one table and get everything in one table (outer merge)
    merged_replicate_tables = []
    for condition, table_list in all_cond_cutoff_tables.items():
        df = reduce(lambda df1,df2: pd.merge(df1,df2,on=['uniprot_id', 'annotation'], how="outer"), table_list)
        merged_replicate_tables.append(df)

    # Merge cutoff ratio from all conditions in one file
    df_all = reduce(lambda df1,df2: pd.merge(df1,df2,on=['uniprot_id', 'annotation'], how="outer"), merged_replicate_tables)
    
    # merge ratio and signal intensity with cutoff values
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
    
    # merge 
    #df_all_list_columns_TPRFPR = df_all.filter(like="TPR-FPR").columns.tolist()
    df_all_list_columns_passcutoff = df_all.filter(like="pass_cutoff").columns.tolist()

    columns_list_df_all = df_all_list_columns_passcutoff #df_all_list_columns_TPRFPR
    columns_list_df_all.append("uniprot_id")
        
    ratio_and_signal_intensity = ratio_and_signal_intensity.rename(columns={"uniprot_id": "Entry"})
    ratio_and_signal_intensity.drop('description', axis=1, inplace=True)
    fasta_table = fasta_table.rename(columns={"uniprot_id": "alias_uniprot_id"})

    merged_with_metadata = pd.merge(ratio_and_signal_intensity, fasta_table, on=["Entry", "annotation"], how="left")
    merged_with_metadata = merged_with_metadata.rename(columns={"Entry": "uniprot_id"})
    merged_with_metadata = merged_with_metadata.drop_duplicates(subset='uniprot_id', keep='first')
    
    ratio_and_signal_intensity_merged = pd.merge(merged_with_metadata, df_all[columns_list_df_all],on=["uniprot_id"], how="outer")
    
    #print(len(ratio_and_signal_intensity_merged))
    #print(len(df_all))
    
    #print("")
    ratio_and_signal_intensity_merged = add_pass_cutoff_analysis_to_df(ratio_and_signal_intensity_merged, df_all_list_columns_passcutoff)
    
    with pd.ExcelWriter(tissue_file_folder_path / ('decision_table' + file_name.split("processed_census-out")[1] + '.xlsx')) as writer:
        ratio_and_signal_intensity_merged.to_excel(writer, sheet_name='ratio_normSI_annot', index=False)
        df_all.to_excel(writer, sheet_name='cutoff_plots_data', index=False)
    
    return ratio_and_signal_intensity_merged, df_all


def get_before_after_cutoff_barplots(decision_table, pass_cutoff_path, file_name):
    pass_cutoff_true_df = decision_table[decision_table.pass_cutoff_result == True]   
    pass_cutoff_cols = pass_cutoff_true_df.filter(like="pass_cutoff").columns.tolist()
    pass_cutoff_true_df = pass_cutoff_true_df.drop(pass_cutoff_cols, axis = 1)
    #pass_cutoff_true_df = pass_cutoff_true_df.reset_index()

    barplot_list = []
    total = len(pass_cutoff_true_df)
    tp_num = len(pass_cutoff_true_df[pass_cutoff_true_df["annotation"] == "TP"])
    fp_num = len(pass_cutoff_true_df[pass_cutoff_true_df["annotation"] == "FP"])
    na_num = len(pass_cutoff_true_df[pass_cutoff_true_df["annotation"].isnull()])
    pass_cutoff_true_df_list = ["pass cutoff 2 \n ("+ str(total)+" proteins)", (tp_num/total)*100, (fp_num/total)*100, (na_num/total)*100]
        
    list_columns_ratio = decision_table.filter(like="ratio_").columns.tolist()
    decision_table = decision_table.set_index(["uniprot_id", "annotation"])
    decision_table = decision_table[list_columns_ratio].dropna()
    decision_table = decision_table.reset_index()
    #sheet_1 = sheet_1[["uniprot", "annotation"]]
    total = len(decision_table)
    tp_num = len(decision_table[decision_table["annotation"] == "TP"])
    fp_num = len(decision_table[decision_table["annotation"] == "FP"])
    na_num = len(decision_table[decision_table["annotation"].isnull()])
    all_proteins = ["all \n ("+ str(total)+" proteins)", (tp_num/total)*100, (fp_num/total)*100, (na_num/total)*100]
        
    barplot_list.append(all_proteins)
    barplot_list.append(pass_cutoff_true_df_list)

    # Create the pandas DataFrame
    df = pd.DataFrame(barplot_list, columns=['state', 'True positives', 'False positives', 'Not annotated'])
    df = df.set_index("state")
    df.plot(kind='bar', stacked=True, color=['green','red','lightgrey'])
    plt.xlabel('')
    plt.ylabel("% Proteins")
    plt.title(file_name.split("processed_census-out_")[1] + " 2 pass cutoff")
    plt.xticks(rotation = 0) 
    #plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.savefig(pass_cutoff_path / (file_name.split("processed_census-out_")[1] +'_barplot.png'))
    plt.show()
        
    return pass_cutoff_true_df

def get_tp_fp_cutoff_plots(cutoff_plots_table, cutoff_roc_path, file_name, cutoff_dict):
    sheet_2 = cutoff_plots_table.set_index(["uniprot_id", "annotation"])
    list_columns_passcutoff = sheet_2.filter(like="pass_cutoff").columns.tolist()
        
    folder_path = cutoff_roc_path / "tp_fp_plots"
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)

    for col in list_columns_passcutoff:
        df_list = []
        col_clean = col.replace("pass_cutoff_", "")
        
        for file, cutoff_vals in cutoff_dict.items():
            if file_name in file:
                cutoff_value = cutoff_vals[col_clean]
                cutoff_value = round(cutoff_value, 3)   
            
        sub_df_cols = sheet_2.filter(like=col_clean).columns.tolist()
        sub_df = sheet_2[sub_df_cols]
        sub_df = sub_df.dropna(subset=["log2_norm_ratio_"+col_clean])
        sub_df = sub_df.reset_index()
        sub_df = sub_df.set_index("uniprot_id")
        sub_df_plot = sub_df[["log2_norm_ratio_"+col_clean, "annotation"]]
        TP = sub_df_plot[sub_df_plot["annotation"] == "TP"].drop("annotation", axis=1).iloc[:, 0].rename("TP")
        FP = sub_df_plot[sub_df_plot["annotation"] == "FP"].drop("annotation", axis=1).iloc[:, 0].rename("FP")
        no_val = sub_df_plot[sub_df_plot["annotation"].isnull()].drop("annotation", axis=1).iloc[:, 0].rename("NA")
        all_list = ["all \n ("+ str(len(sub_df_plot))+" proteins)", (len(TP)/len(sub_df_plot))*100, (len(FP)/len(sub_df_plot))*100, (len(no_val)/len(sub_df_plot))*100]
        df_list.append(all_list)
        
        df_to_plot = pd.concat([TP,FP, no_val],axis=1)
        gfg = sns.histplot(df_to_plot,
                           palette=dict(TP="green", FP="red", NA="grey"))
        plt.axvline(cutoff_value, 0, color="black", linewidth=0.5, linestyle='dashed')
        #plt.text(cutoff_value+1, 100, "Cutoff: "+str(cutoff_value))
        gfg.set(xlabel ="log2("+col_clean+")", title = file_name + "\n" + col_clean + " Cutoff: "+str(cutoff_value))
        plt.show()
        fig = gfg.get_figure()
        description = col_clean.replace("/", "_")
        fig.savefig(folder_path / (description+'hist_1.png')) 
        
        df_to_plot = pd.concat([TP,FP],axis=1)
        gfg_2 = sns.histplot(df_to_plot, 
                             palette=dict(TP="green", FP="red"))
        plt.axvline(cutoff_value, 0, color="black", linewidth=0.5, linestyle='dashed')
        #plt.text(cutoff_value+1, 60, "Cutoff: "+str(cutoff_value))
        gfg_2.set(xlabel ="log2("+col_clean+")", title = file_name + "\n" + col_clean + " Cutoff: "+str(cutoff_value))
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
        plt.title(file_name + "\n" + col_clean)
        plt.xticks(rotation = 0) 
        #plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        description = col_clean.replace("/", "_")
        plt.savefig(folder_path / (description+'barplot.png'))
        plt.show()   
    
    return "done"            