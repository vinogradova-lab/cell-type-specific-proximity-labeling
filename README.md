# TurboID


### **Step 1 (processing of raw census-out files)**

Files were processed with the Vinogradova Lab pipeline and were treated as “whole proteome” files – i.e., we quantify protein expression changes and corresponding peptide values were summed to receive the respective raw signal intensity per protein

--> Processing details
- Files were processed twice: 
  - with the restriction that at least 1 peptide must be present per protein
  - with the restriction that at least 2 peptides must be present per protein 
- Additionally following filters were used:
  - Non unique peptides were removed
  - Peptides that are not fully tryptic were removed
  - Peptides with more than 1 missed cleavage site were removed

Note: make sure to add (copy from old folders) both conditions_metadata.csv and metadata_col.csv - these files are needed for renaming the channels as well as specifying which are cre- and cre+ conditions

### **Step 2a (mouse_fasta_to_FP_TP_tables_without_signalp.py)**
- create main table (based off of our mouse fasta file used for IP2 search)
  - read in fasta file to get all uniprot IDs
  - use uniprot.org to get all annotations for these proteins
  - label proteins as TP based on following criteria in subcellular location column: 
    - "Endoplasmic Reticulum", "Secreted", "secreted", "Endoplasmic reticulum", "Rough endoplasmic reticulum", "endoplasmic reticulum"
  - label proteins as FP if uniprot id is found in human mitomatrix list from supplementary table S1 [Rhee et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3916822/)
    - obtain human-mouse gene orthologs by following tutorial: [How to get all the orthologous genes between two species](https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/) 
    - this allows us to get the correct mouse gene names which correspond to the mitomatrix protein list, we then crossreference our fasta table and label them as FP
  - additionally make sure that if there is a SignalP annotation (ie if SignalP column is not empty) for a protein marked as FP then remove the FP annotation, also if protein is marked as both TP and FP, keep it in the TP list and remove it from FP list
  - crossreference table with secretion prediction resource shared by Corey 
  - crossreference main list with spleen, adipose tissue etc. (files provided by Ken)
  - crossreference table with Crapome database
  - import NCBI gene annotations (using goatools package) and crossreference gene symbols to add gene_ids into our final list (for go term analysis)
  - save table
  - Number of TP proteins: 2806
  - Number of FP proteins: 435

### **Step 2b (mouse_fasta_to_FP_TP_tables_with_signalp.py)**
- create main table (based off of our mouse fasta file used for IP2 search)
  - read in fasta file to get all uniprot IDs
  - use uniprot.org to get all annotations for these proteins
  - label proteins as TP based on following criteria in subcellular location column: 
    - "Endoplasmic Reticulum", "Secreted", "secreted", "Endoplasmic reticulum", "Rough endoplasmic reticulum", "endoplasmic reticulum"
  - label proteins as FP if uniprot id is found in human mitomatrix list from supplementary table S1 [Rhee et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3916822/)
    - obtain human-mouse gene orthologs by following tutorial: [How to get all the orthologous genes between two species](https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/) 
    - this allows us to get the correct mouse gene names which correspond to the mitomatrix protein list, we then crossreference our fasta table and label them as FP
  - add table shared by Corey "secretion prediction resource" to main table and check if if sp_score > 0.7 and ups_score == 0 is True otherwise exclude the proteins from TP list
  - if protein is marked as both TP and FP, keep it in the TP list and remove it from FP list
  - crossreference table with secretion prediction resource shared by Corey 
  - crossreference main list with spleen, adipose tissue etc. (files provided by Ken)
  - crossreference table with Crapome database
  - import NCBI gene annotations (using goatools package) and crossreference gene symbols to add gene_ids into our final list (for go term analysis)
  - save table
  - Number of TP proteins: 1602
  - Number of FP proteins: 438


### **Step 3 (turboid_analysis_pipeline.py)**

  - `conda activate turbo_id_notebook`
  - `python turboid_analysis_pipeline.py` (in tuboid/scripts folder)

  - read in data and metadata 
  
**01_raw_files_with_correct_channel_names_keratins_removed_TP_FP_annotated**
  - read in each experiment, rename channels, remove keratins, annotate if TP or FP
  - save experiment table in 01_raw_files_with_correct_channel_names_keratins_removed_TP_FP_annotated

**02_filtering_per_cond_per_file**
   - filter each experiment based on cre+ channels for each condition
  - At least one pair (any 2 channels) out of all cre+ channels should have a min raw signal intensity sum >= 10000 and a coefficient of variation <= 0.5
  - if all conditions did not pass filter, the corresponding protein is removed from list, but if at least one condition for this protein passed this filter then the protein passes and is reported
  
**03_normalization_plots and 04_normlized**
  - Normalization is performed for cre+ channels for each condition
  - Per cre+ channel TP proteins with 2+ peptides were used to calculate the median per channel 
  - normalization factors were calculated by dividing (median of cre+ channel medians)/(median cre+ channel)
  - each cre+ channel was then multiplied with the corresponding normalization factor 

**05_results**
1. Tissue samples 
   - each sample has a 01_cutoff_roc_filter folder 
   - this folder contains a huge excel table with two sheets (for more details on Cutoff ROC analysis please see below)
   - sheet "ratio_normSI_anot" holds all ratios, the normalized SIs from the previous step, detailed annotation and pass cutoff columns 
   - sheet "cutoff_plots_data" contains normalized ratios (median of FP proteins), log2 of that normalized ratio as well as TPR-FPR columns - this data corresponds to plots in folder cutof_ROC_plots, based on this data we decide if a protein passes or not (highest TPR-FPR value is our treshold in the log2 ratio column), this sheet will be of help if anyone is interested to see which other proteins were very close to the treshold but didn't pass etc. 
   - the folder called tp_fp_plots contains visualization of number of proteins that pass and do not pass cutoff and if they are in the TP or FP or no group
   - a second filtering step determines our final protein list which is labelled as final_protein_list_xx
   - this final list contains proteins that pass cutoffs with at least 2 replicates per condition in at least one condition of this file (see filtering in 01_cutoff_roc_filter, excel table, sheet "ratio_normSI_annot" - pass_cutoff_result column)
   - based on this final list volcano plots are created, volcano plot data can be found in the final table
   - and you can also find a barplot indicating the abundance of TP and FP proteins in that final list compared to a barplot without any filtering (ie with what we started initially)
2. Serum samples
   - also here each sample has its own folder 
   - but with serum samples we do not perform any cutoff roc analysis we only add fold change columns (median cre+ / median cre-) per condition for easier filtering in final table
   - there is no enrichment filtering appplied to this data, all identified proteins are shown in volcano plot 
   - the final table is annotated and volcano plots are created based on that table, volcano plot data can be found in the final table

Both serum and tissue sample folders contain heatmap_plot, volcano_plot and goterm analysis folders. 

**GO term analysis**

Go term analysis with significant up/down/up and down regulated proteins - implement in python using goatools
   - Using full mouse genome as background
   - Using experiment specific background (before ROC)
 
Each file folder has a separate goterm_analysis folder which holds the plots and tables. The names of the files follow this scheme:
goterm_BACKGROUND_CONDITION _DOWN/UP/DOWN&UP_go_term_plot.csv or pdf
 
More details: GOEA analysis was performed on the significant up or down or up & down regulated proteins (mapped to gene IDs). Enriched ontologies with a Benjamini-Hochberg-corrected p-value less than 0.05 were retained and top 10 terms were plotted. I am happy to improve the plots as I think the color and shape definitely need improvement but as far as information, they show number of genes and FDR which I hope will be helpful in interpreting them. I also provide corresponding tables which also show the names of the genes for each term. I also want to note that there are some proteins that I could not map to gene ids, which means if they are significantly up or down regulated they are not included in the goterm analysis because they don’t have a gene id  – I will need to look more into this, and will keep you all updated but for now I just wanted to share what I have. Some plots are empty which means that nothing passed the p-value 0.05 threshold.

**Volcano plots**

For every protein the log2 fold-change of the ratio between condtions was calculated. p-values were calculated with T-test for the means of two independent samples for all replicates and volcano plots show uncorrected p-values. Significant proteins which show a > 1.5 fold change are highlighted in either red or green. 

**Heatmap**

A heatmap of top proteins is provided as pdf (can be searched cmd+F) and can be compared to volcano plots

### Details Cutoff ROC analysis 
Following the approach of: [Cho et al. 2020](https://www.nature.com/articles/s41596-020-0399-0)

1. Calculate Ratios per condition (example of one condition) 
   - Each condition within each file was analyzed on its own, therefore each condition will have its own cutoff in the end. 
   - Example ratios of file EV1-100A for condition iWAT_HF
        - ratio_iWAT_HF_cre+_1/iWAT_HF_cre-1
        - ratio_iWAT_HF_cre+_1/iWAT_HF_cre-_2 	
        - ratio_iWAT_HF_cre+_2/iWAT_HF_cre-_1 	
        - ratio_iWAT_HF_cre+_2/iWAT_HF_cre-_2 	
        - ratio_iWAT_HF_cre+_3/iWAT_HF_cre-_1 	
        - ratio_iWAT_HF_cre+_3/iWAT_HF_cre-_2 
   - Each replicate will be treated on its own (all steps according to Cho et al - Nature Protocol)
   - Example: ratio_ iWAT_HF_cre+_1/iWAT_HF_cre-1
   - Normalize ratio_ iWAT_HF_cre+_1/iWAT_HF_cre-1 by the median ratio of proteins that are labelled as FP
   - Rank ratios in descending order 
   - For each row in table calculate number of TP and FP in all the rows before
   - Calculate TPR and FPR 
   - Plot TPR and FPR (see ROC plot) 
   - calculate TPR-FPR 
   - find cutoff ratio which has the highest value for TPR-FPR
   - Plot log2 ratio and TPR-FPR (see line plot – cutoff will be the highest point) 

2. Now that cutoffs for each ratio are identified we merge this information with our original processed raw intensity values file to keep an overview per file and add additional annotation columns from uniprot
   - Sheet “ratio_normSI_annot”:
     - “ratio_” columns: ratios of replicates 
     - iWAT_HF_cre+_1, iWAT_HF_cre+_2 etc.: normalized signal intensity values
     - “pass_cutoff” columns: 1/0  = passed cutoff/did not pass cutoff
     - “median_R” columns: average (median) of ratios for a particular condition (cre+)
     - “median_SI” columns: average (median) of signal intensity for a particular condition (cre+)
     - “condition_pass_cutoff_sum_min_2_cutoffs” columns: evaluates if this condition has at least two different cre+ replicates that pass the ratio cutoff
     - “pass_cutoff_result” columns: proteins that are marked as True represent our final protein list (volcano plots are created using this list). This column is True if any of the conditions are True for the “condition_pass_cutoff_sum_min_2_cutoffs” columns 
   - Sheet “cutoff_plots_data”:
     - “log2_norm_ratio_” columns: log2 normalized ratio (x axis in line plot) 
     - “norm_ratio_” columns: normalized ratio (same as above but without log2) 
     - “TPR-FPR” columns: True positive rate – False positive rate for this protein (y axis in line plot)
     - “pass_cutoff” columns: 1/0  = passed cutoff/did not pass cutoff


**Notes for file: processed_census-out_04172023_CRW_A-5_16pl_M**
- condition scheme does not match other files so needed to make follwing changes:
  - usually we compare the two conditions in a file in the end this is not the case here so skipping final pass_cutoff filter - which checks if either condition has passed with min 2 replicates - for this file specifically we treat GFAP and BAT individually thats why each condition has its own final table (BAT - min 2 replicates must have passed cutoff to be in final protein list (-> volcano plot), GFAP: at least one of the 3 ratios must have passed to be in final protein list (->scatterplots)) 
  - for GFAP cre+ filtering is done on all 3 cre+ channels (ctrl, fast, LPS) if at least 2 pass sum and cv then protein passes, Cre(-)_Ctrl_W6 is the control for all Cre(+) samples 
  - since there are no replicates for GFAP - creating scatterplots to compare fold change for GFAP  


**Notes for comparing files 20230512_EV2-28A (male) vs. 20230513_EV2-28B (female):**
- First we processed both files individually 
- Using the filtered and normalized files perform normalization to sum per protein (channel ratio)
  
 **1**

- Then male and female files were merged 
- And PCA plot was created 
- Proteins detected in male file: 1231 
- Proteins detected in female file: 1297
- Only ~300 overlap

**2**

- Using filtered, normalized and roc result
- Normalization to sum per protein (channel ratio)
- Then male and female files were merged 
- And PCA plot was created 
- Proteins detected in male file: 745 
- Proteins detected in female file: 1011
- And only 170 overlap after ROC 

Todo: investigate overlaping proteins

