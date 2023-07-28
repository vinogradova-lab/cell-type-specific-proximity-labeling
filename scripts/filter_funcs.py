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