import pandas as pd
import numpy as np
from Bio import SeqIO

def remove_duplicates_from_aggregate(row): 
    list_of_cols = row.index.tolist()
    for col in list_of_cols: 
        row[col] = ', '.join(map(str, set(row[col])))
    return row

def annotate_TP(row):
    tp_term_list = ["Endoplasmic Reticulum", "Endoplasmic reticulum", "Secreted", "secreted", "Rough endoplasmic reticulum", "endoplasmic reticulum"]
    subcelullar_location = str(row["Subcellular location [CC]"])
    #"Golgi", "Extracellular space", "extracellular space", 
    subcelullar_location = str(subcelullar_location)
    if "Note:" in subcelullar_location:
        subcelullar_location = subcelullar_location.split("Note:")[0]
    for item in tp_term_list: 
        if item in subcelullar_location:
            return("TP")

def annotate_TP_signalp(row):
    tp_term_list = ["Endoplasmic Reticulum", "Endoplasmic reticulum", "Secreted", "secreted", "Rough endoplasmic reticulum", "endoplasmic reticulum"]
    # "Golgi", "Extracellular space", "extracellular space", 
    subcelullar_location = str(row["Subcellular location [CC]"])
    sp_score = row["Score"]
    ups_score = row["UPS_Score"]
    if "Note:" in subcelullar_location:
        subcelullar_location = subcelullar_location.split("Note:")[0]
    for item in tp_term_list: 
        if item in subcelullar_location:
            if sp_score > 0.7 and ups_score == 0: 
                return "TP"

def annotate_FP_subcellular_loc(subcelullar_location):
    fp_term_list = ["Mitochondrion", "mitochondrion"] # removed: 'cytoskeleton' "Cytoskeleton", "Nucleus", "nucleus", 
    subcelullar_location = str(subcelullar_location)
    if "Note:" in subcelullar_location:
        subcelullar_location = subcelullar_location.split("Note:")[0]
    for item in fp_term_list: 
        if item in subcelullar_location:
            return "FP" 
            

def annotate_FP_mitomatrix(gene_names, mitomatrix_list):
    gene_name_list = str(gene_names).split(' ')
    for gene in gene_name_list:
        if gene in mitomatrix_list:
            return "FP"

def conclude_annotation(df):   
    if df["TP"] == "TP":
        return "TP" 
    # if signal peptide is present then even though annotation says FP, change annotation to nan
    elif df["FP"] == "FP" and pd.notna(df["Signal peptide"]):
        return np.nan
    elif df["FP"] == "FP" and pd.isnull(df["Signal peptide"]):
        return "FP"

def conclude_annotation_signalp(df):   
    sp_score = df["Score"]
    ups_score = df["UPS_Score"]

    if df["TP"] == "TP":
        return "TP" 
    elif df["FP"] == "FP":
        if sp_score > 0.7 and ups_score == 0:  
            return np.nan
        else: 
            return "FP"

def check_uniprot(row, sp_list): 
    if row.uniprot_id in sp_list:
        return row.uniprot_id
    elif "," in row.Entry: 
        for id in row.Entry.split(", "): 
            if id in sp_list: 
                return id
    else:
        return np.nan

# create main fasta table (dataframe) from fasta entries, remove reverse entries, set uniprot id as index 
def create_df_from_fasta(mouse_fasta_path):
    # read in original massspec fasta file as dict 
    fasta_dict = {record.id : record.description for record in SeqIO.parse(mouse_fasta_path, "fasta")}
    main_fasta_table = pd.DataFrame.from_dict(fasta_dict, orient='index')
    main_fasta_table = main_fasta_table.reset_index()
    main_fasta_table.columns = ["uniprot_id", "description"]
    main_fasta_table["id"] = main_fasta_table["uniprot_id"]
    main_fasta_table[["db", "uniprot_id", "entry_name"]] = main_fasta_table["uniprot_id"].str.split("|", expand=True)
    # filter out reverse entires
    main_fasta_table = main_fasta_table[~main_fasta_table["id"].str.contains("Reverse_")]
    main_fasta_table = main_fasta_table.set_index("uniprot_id")
    return main_fasta_table

def import_uniprot_annot(uniprot_id_mapping_path):
    # read in uniprot id mapping results, annotate keratins
    uniprot_mapping = pd.read_csv(uniprot_id_mapping_path)
    uniprot_mapping["keratin"] = uniprot_mapping["Protein names"].str.contains("keratin", case=False)
    uniprot_mapping.rename({"From":"uniprot_id"}, axis=1, inplace=True)
    uniprot_mapping = uniprot_mapping.set_index("uniprot_id")
    return uniprot_mapping

def clean_up_fasta_table_merge(merged_main_fasta_table):
    # change dtypes so that we can concat later 
    merged_main_fasta_table = merged_main_fasta_table.astype({'Length': 'object', 'Mass': 'object'})
    # clean
    duplicates_list = merged_main_fasta_table[merged_main_fasta_table.index.duplicated()].index.tolist()
    duplicates_df = merged_main_fasta_table[merged_main_fasta_table.index.isin(duplicates_list)]
    deduplicated_entries = duplicates_df.groupby("uniprot_id").agg(list) #.reset_index()
    deduplicated_entries = deduplicated_entries.apply(remove_duplicates_from_aggregate, axis=1)
    deduplicated_entries["Gene Names"] = deduplicated_entries["Gene Names"].str.replace(",", "")

    merged_main_fasta_table = merged_main_fasta_table.drop(index=duplicates_list)
    merged_main_fasta_table = pd.concat([merged_main_fasta_table, deduplicated_entries])
    return merged_main_fasta_table

def add_secretion_prediction_to_fasta_table(merged_main_fasta_table, secretion_prediction_resource_path):
    spr_df = pd.read_csv(secretion_prediction_resource_path)
    spr_df = spr_df.rename(columns={"Accession":"uniprot_sp_key", "Organism":"organism"})
    sp_list = spr_df.reset_index().uniprot_sp_key.tolist()
    merged_main_fasta_table = merged_main_fasta_table.reset_index()
    merged_main_fasta_table['Entry'] = merged_main_fasta_table['Entry'].astype(str)
    merged_main_fasta_table['uniprot_id'] = merged_main_fasta_table['uniprot_id'].astype(str)
    merged_main_fasta_table["uniprot_sp_key"] = merged_main_fasta_table.apply(check_uniprot, args=(sp_list,), axis=1)
    merged_main_fasta_table = merged_main_fasta_table.merge(spr_df, on="uniprot_sp_key", how="left")
    merged_main_fasta_table = merged_main_fasta_table.set_index('uniprot_id')
    return merged_main_fasta_table

def get_human_mitomatrix_mouse_ortholog(mitomatrix_protein_path, mart_export_path): 
    mitomatrix_df = pd.read_csv(mitomatrix_protein_path)
    mitomatrix_df['Gene Names'] = mitomatrix_df['Gene Names'].apply(lambda x: x.split(";"))
    mitomatrix_df = mitomatrix_df.explode("Gene Names")
    mitomatrix_df = mitomatrix_df.rename(columns={"Gene Names":"Gene name"})

    # we need to be careful when crossreferencing by gene names
    # https://www.biostars.org/p/149115/
    # https://pypi.org/project/mousipy/
    # https://bioinformatics.stackexchange.com/questions/17486/converting-mouse-genes-to-human-genes

    # follow https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/ turotial and read in results
    mart_human_mouse_orthologs = pd.read_csv(mart_export_path)

    # merge mart results to get mouse gene list 
    mitomatrix_mart_merge_df = mitomatrix_df.merge(mart_human_mouse_orthologs, how="left", on="Gene name")
    mouse_mitomatrix_genes_list = mitomatrix_mart_merge_df["Mouse gene name"].tolist()
    return mouse_mitomatrix_genes_list

def annotate_TP_FP(with_signalp, merged_main_fasta_table, mouse_mitomatrix_genes_list):
    if with_signalp: 
        merged_main_fasta_table["TP"] = merged_main_fasta_table.apply(annotate_TP_signalp, axis=1) 
    else:
        merged_main_fasta_table["TP"] = merged_main_fasta_table.apply(annotate_TP, axis=1) 

    merged_main_fasta_table["FP"] = merged_main_fasta_table["Gene Names"].apply(annotate_FP_mitomatrix, args=(mouse_mitomatrix_genes_list, )) 
    # print(merged_main_fasta_table["TP"].value_counts())
    # print(merged_main_fasta_table["FP"].value_counts())
    # clean up annotation 
    if with_signalp: 
        merged_main_fasta_table["annotation"] = merged_main_fasta_table.apply(conclude_annotation_signalp, axis=1) 
    else: 
        merged_main_fasta_table["annotation"] = merged_main_fasta_table.apply(conclude_annotation, axis=1) 
    # merged_main_fasta_table["annotation"].value_counts()
    return merged_main_fasta_table

def add_tissue_annotation_from_ken(tissue_localisation_path, merged_main_fasta_table):
    # check spleen, adipose tissue, etc file names
    files_in_folder = list(tissue_localisation_path.iterdir())
    #[file_path.stem for file_path in files_in_folder]
    for file_path in files_in_folder:
        column_name = file_path.stem.split("_", 1)[0]
        df = pd.read_excel(file_path, engine='openpyxl')
        entry_list = df.Entry.tolist()
        merged_main_fasta_table[column_name] = merged_main_fasta_table.index.isin(entry_list)
    return merged_main_fasta_table

def add_gene_id_mapping(gprofiler_mouse_geneid_mapping_path, merged_main_fasta_table):
    gene_id_annot_df = pd.read_csv(gprofiler_mouse_geneid_mapping_path)
    gene_id_annot_df.columns = ['uniprot_id', 'gene_id', 'gene_name', 'description', 'namespace']
    gene_id_annot_df = gene_id_annot_df[['uniprot_id', 'gene_id']].drop_duplicates()
    gene_id_annot_df = gene_id_annot_df.groupby("uniprot_id").agg(list)
    for col in gene_id_annot_df.columns: 
        gene_id_annot_df[col] = [','.join(map(str, l)) for l in gene_id_annot_df[col]]
        gene_id_annot_df[col] = gene_id_annot_df[col].str.split(",").str[0]
        gene_id_annot_df[col] = gene_id_annot_df[col].replace("None", np.nan)

    merged_main_fasta_table = merged_main_fasta_table.join(gene_id_annot_df)
    return merged_main_fasta_table

def add_human_uniprot_mapping(human_mouse_homology_mapping_path):
    human_mouse_homology_mapping_df = pd.read_csv(human_mouse_homology_mapping_path, delimiter="\t")
    human_mouse_homology_mapping_df = human_mouse_homology_mapping_df[["DB Class Key", "Common Organism Name", "SWISS_PROT IDs", "Symbol"]].drop_duplicates()
    mouse_df = human_mouse_homology_mapping_df.loc[human_mouse_homology_mapping_df["Common Organism Name"] == 'mouse, laboratory'].set_index("DB Class Key").add_suffix("_mouse")
    human_df = human_mouse_homology_mapping_df.loc[human_mouse_homology_mapping_df["Common Organism Name"] == 'human'].set_index("DB Class Key").add_suffix("_human")
    human_df = human_df.reset_index().groupby("DB Class Key").agg(list)
    for col in human_df.columns: 
        human_df[col] = [','.join(map(str, l)) for l in human_df[col]]

    human_mouse_homology_mapping_rearranged_df = mouse_df.join(human_df)
    human_mouse_homology_mapping_rearranged_df = human_mouse_homology_mapping_rearranged_df[["SWISS_PROT IDs_mouse", "Symbol_mouse", "SWISS_PROT IDs_human", "Symbol_human"]].reset_index().add_suffix("_jackson_homology_db")

    maping_dict = human_mouse_homology_mapping_rearranged_df[["Symbol_mouse_jackson_homology_db", "DB Class Key_jackson_homology_db", "SWISS_PROT IDs_mouse_jackson_homology_db"]]
    gene_name_mapping_dict = pd.Series(maping_dict["DB Class Key_jackson_homology_db"].values,index=maping_dict["Symbol_mouse_jackson_homology_db"]).to_dict() 

    maping_dict["SWISS_PROT IDs_mouse_jackson_homology_db"] = maping_dict["SWISS_PROT IDs_mouse_jackson_homology_db"].str.split(",")
    maping_dict = maping_dict.explode("SWISS_PROT IDs_mouse_jackson_homology_db")
    maping_dict = maping_dict.drop_duplicates().dropna()
    uniprot_mapping_dict = pd.Series(maping_dict["DB Class Key_jackson_homology_db"].values,index=maping_dict["SWISS_PROT IDs_mouse_jackson_homology_db"]).to_dict() 
    return gene_name_mapping_dict, uniprot_mapping_dict, human_mouse_homology_mapping_rearranged_df

def add_db_key(df, uniprot_mapping_dict, gene_name_mapping_dict):
    gene_name = df["gene_name"] 
    uniprot_id = df["uniprot_id"]
    if uniprot_id in uniprot_mapping_dict.keys():
        return uniprot_mapping_dict[uniprot_id]
    elif gene_name in gene_name_mapping_dict: 
        return gene_name_mapping_dict[gene_name]

def clean_human_uniprot(str_of_uniprots):
    uniprot_list = str(str_of_uniprots).split(",")
    clean_list = []
    for uniprot in uniprot_list: 
        if uniprot != 'nan':
            if uniprot not in clean_list:
                clean_list.append(uniprot)

    return(", ".join(clean_list))