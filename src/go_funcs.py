import pandas as pd
import numpy as np
from goatools.cli.ncbi_gene_results_to_python import ncbi_tsv_to_py
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
#from __future__ import print_function
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import seaborn as sns
import textwrap
from pathlib import Path
import sys
import logging

def read_in_ncbi_go_associations_data():
    # NCBI Mouse data and GO data was Downloaded on May 17 2023
    # read in ncbi data into py 
    ncbi_txt = 'gene_result.txt'
    output_py = 'genes_ncbi_10090_proteincoding.py'

    if not Path("./" + ncbi_txt).exists():
        print('You need to download NCBI data and save it as gene_result.txt')
        print('Go to:')
        print("https://www.ncbi.nlm.nih.gov/gene")
        print("Enter:")
        print('10090"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties]')
        print("if you need human change the number to 9606")
        print("send to file")
        print("save in current dir and try again!")
        sys.exit()
    
    if not Path("./" + output_py).exists():
        ncbi_tsv_to_py(ncbi_txt, output_py)

    # import py as gneid2nt_mus
    from genes_ncbi_10090_proteincoding import GENEID2NT as GeneID2nt_mus

    # map gene ids to gene name
    mapper = {}
    for key in GeneID2nt_mus:
        mapper[GeneID2nt_mus[key].Symbol] = GeneID2nt_mus[key].GeneID 
    inv_map = {v: k for k, v in mapper.items()}
        
    # download go data
    try:
        obo_fname
    except NameError:
        obo_fname = download_go_basic_obo()
        
    # download ncbi associations
    try:
        file_gene2go
    except NameError:
        file_gene2go = download_ncbi_associations()
    
    return GeneID2nt_mus, inv_map, file_gene2go


def initialize_godag_obj(file_gene2go):
    # Initialize GODag object
    obodag = GODag("go-basic.obo")
    logging.info("Loading Ontologies...")
    logging.info("go-basic.obo from Gene Ontology Consortium website: %s GO terms", len(obodag))

    # Read NCBI's gene2go. Store annotations in a list of namedtuples
    objanno = Gene2GoReader(file_gene2go, taxids=[10090])

    # Get associations for each branch of the GO DAG (BP, MF, CC)
    ns2assoc = objanno.get_ns2assc()

    logging.info("Loading Associations...")
    for nspc, id2gos in ns2assoc.items():
        logging.info("{NS} {N:,} annotated mouse genes".format(NS=nspc, N=len(id2gos)))

    # subset mouse genome to only proteins/genes in experiment 
    return obodag, ns2assoc


def create_godag_obj(obodag, ns2assoc, GeneID2nt_mus, reference_list=None):
    logging.info("Loading Background gene set...")
    if reference_list: 
        reference_list_wona = [item for item in reference_list if str(item) != 'nan']
        GeneID2nt_mus = dict(filter(lambda item: item[0] in reference_list_wona, GeneID2nt_mus.items()))
        logging.info('Number of items in background gene set (experiment specific set): %s', len(reference_list_wona))
    else:
        logging.info('Number of items in background gene set (NCBI all mouse protein coding genes): %s', len(GeneID2nt_mus))

    # create go object this needs to be done only once if whole mouse genome is used
    # if individual reference list are used then it needs to be done for every file
    logging.info("Initializing GOEA object which holds the Ontologies, Associations, and background...")
    goeaobj = GOEnrichmentStudyNS(GeneID2nt_mus, # List of mouse protein-coding genes
                                  ns2assoc, # geneid/GO associations
                                  obodag, # Ontologies
                                  propagate_counts = False,
                                  alpha = 0.05, # default significance cut-off
                                  methods = ['fdr_bh']) # defult multipletest correction method
                         
    return goeaobj

def get_all_goterms(goeaobj):
    # create list that holds all GO terms
    # we will need this to count the number of times this go term appears in the mouse genome 
    GO_items = []

    temp = goeaobj.ns2objgoea['BP'].assoc
    for item in temp:
        GO_items += temp[item]
        

    temp = goeaobj.ns2objgoea['CC'].assoc
    for item in temp:
        GO_items += temp[item]
        

    temp = goeaobj.ns2objgoea['MF'].assoc
    for item in temp:
        GO_items += temp[item]
    
    return GO_items

def go_it(test_genes, goeaobj, GO_items, inv_map):
    logging.info(f'Number of input GeneIDs for GO analysis: {len(test_genes)}')
    
    goea_results_all = goeaobj.run_study(test_genes)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    
    df = pd.DataFrame(list(map(lambda x: [x.GO, 
                                          x.goterm.name, 
                                          x.goterm.namespace, 
                                          x.p_uncorrected, 
                                          x.p_fdr_bh,
                                          x.ratio_in_study[0], 
                                          x.ratio_in_study[1], 
                                          GO_items.count(x.GO), 
                                          list(map(lambda y: inv_map[y], x.study_items)),], 
                                          goea_results_sig)), 
                                          columns = ['GO', 
                                                     'term', 
                                                     'class', 
                                                     'p-value', 
                                                     'p_fdr_bh', 
                                                     'n_genes',
                                                     'n_study', 
                                                     'n_go', 
                                                     'study_genes'])

    df = df[df.n_genes > 1]
    df["-log10(FDR)"] = -1*np.log10(df["p_fdr_bh"])
    df = df.sort_values(by='-log10(FDR)', ascending=False)
    df['fold_enrichment'] = (df.n_genes/df.n_go) * 100 
    
    return df

#sns.set(font_scale=1.5)
#sns.set_style("ticks")

def create_go_plots(df):
    class_list = list(set(df['class'].tolist()))
    class_list_sorted = sorted(class_list)
    
    fig, axs = plt.subplots(nrows=3, figsize=(5,10))
    fig.tight_layout(pad=4.0)
    count = 0 

    for go_class in class_list_sorted:
        sub_df = df.loc[df['class'] == go_class]
        sub_df = sub_df.head(10)
    
        # = mpl.cm.Blues
        #min_val, max_val, n = 0.3, 1.0, 10
        #colors = cmap(np.linspace(min_val, max_val, n))
        #cmap = mpl.colors.LinearSegmentedColormap.from_list("cmap", colors)

        #norm = mpl.colors.Normalize(vmin = sub_df['-log10(FDR)'].min(), vmax = sub_df['-log10(FDR)'].max())
        #mapper = cm.ScalarMappable(norm = norm, cmap = cmap)
        #mapper.set_array([]) 

        #minsize = min(sub_df['n_genes'])
        #maxsize = max(sub_df['n_genes'])
        bp = sns.scatterplot(data = sub_df, 
                             x = '-log10(FDR)', 
                             y = 'term',
                             size = 'n_genes',
                             sizes=(10, 100),
                             hue = "n_genes",
                             hue_norm=(0, 50),
                             #palette = mapper.to_rgba(sub_df["-log10(FDR)"].values),
                             ax=axs[count],
                             zorder=1)
        bp.set_yticklabels([textwrap.fill(e, 40) for e in sub_df['term']])
        bp.set_title(go_class.replace("_", " ").title())
        bp.set_ylabel("")
        bp.legend(loc='upper left', bbox_to_anchor=(1, 1), labelspacing=1, title='N. of genes')
        xlabels = ['{:,.1f}'.format(x) for x in bp.get_xticks()]
        bp.set_xticklabels(xlabels)
        
        #cbar = plt.colorbar(mapper, shrink=0.50, ax=axs[count], location='right')
        #cbar.set_label('-log10(FDR)', rotation=270,labelpad=30)
        count += 1
        
    return fig    


def get_up_down_goterm(pass_cutoff_true_df, goeaobj, GO_items, inv_map, folder_path, reference):
    reg_cols = pass_cutoff_true_df.filter(regex="Regulation_").columns.tolist()
    cols_to_select = ["gene_id"] + reg_cols
    subset_df = pass_cutoff_true_df[cols_to_select]
    subset_df = subset_df.reset_index().set_index(["uniprot_id", "gene_id"])
    cols_list = [col.split(" - ")[1].split(" (")[0] for col in subset_df]
    unique_col_list = list(set(cols_list))

    for comparison in unique_col_list: 
        logging.info(comparison)
        comparison_df = subset_df.filter(like=comparison)
        reg_col = [col for col in comparison_df if col.startswith('Regulation_')]

        # significant up 
        sign_up = comparison_df.loc[comparison_df[reg_col[0]] == "Significant Up"]
        up_gene_ids_list = sign_up.reset_index().gene_id.tolist()
        up_gene_ids_list_wona = [item for item in up_gene_ids_list if str(item) != 'nan']
        logging.info("nas in UP reg list: %s", len(up_gene_ids_list) - len(up_gene_ids_list_wona))

        df = go_it(up_gene_ids_list_wona, goeaobj, GO_items, inv_map)
        fig = create_go_plots(df)
        df.to_csv(folder_path / ("goterm_" + reference + '_' + comparison + "_up_go_term_results.csv"))
        fig.savefig(folder_path / ("goterm_" + reference + '_' + comparison + "_up_go_term_plot.pdf"), bbox_inches='tight')

        # significant down
        sign_down = comparison_df.loc[comparison_df[reg_col[0]] == "Significant Down"]
        down_gene_ids_list = sign_down.reset_index().gene_id.tolist()
        down_gene_ids_list_wona = [item for item in down_gene_ids_list if str(item) != 'nan']
        logging.info("nas in DOWN reg list: %s", len(down_gene_ids_list) - len(down_gene_ids_list_wona))

        df = go_it(down_gene_ids_list_wona, goeaobj, GO_items, inv_map)
        fig = create_go_plots(df)
        df.to_csv(folder_path / ("goterm_" + reference + '_' + comparison + "_down_go_term_results.csv"))
        fig.savefig(folder_path / ("goterm_" + reference + '_' + comparison + "_down_go_term_plot.pdf"), bbox_inches='tight')

        # significant up and down 
        sign_up_down = comparison_df.loc[comparison_df[reg_col[0]].isin(["Significant Up", "Significant Down"])]
        up_down_gene_ids_list = sign_up_down.reset_index().gene_id.tolist()
        up_down_gene_ids_list_wona = [item for item in up_down_gene_ids_list if str(item) != 'nan']
        logging.info("nas in UP and DOWN reg list: %s", len(up_down_gene_ids_list) - len(up_down_gene_ids_list_wona))

        df = go_it(up_down_gene_ids_list_wona, goeaobj, GO_items, inv_map)
        fig = create_go_plots(df)
        df.to_csv(folder_path / ("goterm_" + reference + '_' + comparison + "_up_and_down_go_term_results.csv"))
        fig.savefig(folder_path / ("goterm_" + reference + '_' + comparison + "_up_and_down_go_term_plot.pdf"), bbox_inches='tight')
    
    return "done"