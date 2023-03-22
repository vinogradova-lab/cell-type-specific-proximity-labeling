#FUNCTIONS for cutoff analysis
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