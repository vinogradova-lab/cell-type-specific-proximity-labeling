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

### **Step 2 (mouse_fasta_to_FP_TP_tables.py)**
- create main table (based off of our mouse fasta file used for IP2 search)
  - read in fasta file to get all uniprot IDs
  - use uniprot.org to get all annotations for these proteins
  - label proteins as TP based on following criteria in subcellular location column: 
    - "Endoplasmic Reticulum", "Secreted", "secreted", "Endoplasmic reticulum", "Rough endoplasmic reticulum", "endoplasmic reticulum"
  - label proteins as FP if uniprot id is found in human mitomatrix list from supplementary table S1 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3916822/)
    - obtain human-mouse gene orthologs by following tutorial: https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
    - this allows us to get the correct mouse gene names which correspond to the mitomatrix protein list, we then crossreference our fasta table and label them as FP
  - additionally make sure that if there is a SignalP annotation (ie if SignalP column is not empty) for a protein marked as FP then remove the FP annotation, also if protein is marked as both TP and FP, keep it in the TP list and remove it from FP list
  - crossreference main list with spleen, adipose tissue etc. (files provided by Ken)
  - save table
  - Number of TP proteins: 2806
  - Number of FP proteins: 261


### **Step 3 (turboid_analysis_pipeline.py)**

  - read in data and metadata 
  
**01_raw_files_with_correct_channel_names_keratins_removed_TP_FP_annotated**
  - read in each experiment, rename channels, remove keratins, annotate if TP or FP
  - save experiment table in 01_raw_files_with_correct_channel_names_keratins_removed_TP_FP_annotated

**02_filtering_per_cond_per_file**
   - read in each experiment and filter based on cre+ channels for each condition 
  - At least one pair (any 2 channels) out of the 3 or 4 cre+ channels should have a min raw signal intensity sum >= 10000 and a coefficient of variation <= 0.5
  - if all conditions did not pass filter, the corresponding protein is removed from list, but if at least one condition for this protein passed this filter then the protein passes and is reported
  
**03_normalization_plots and 04_normlized**
  - Normalization is performed for cre+ channels. 
  - Per cre+ channel TP proteins with 2+ peptides were used to calculate the median per channel 
  - normalization factors were calculated by dividing (median of cre+ channel medians)/(median cre+ channel)
  - each cre+ channel was then multiplied with the corresponding normalization factor 

**05_results**
1. Tissue samples 
   - each sample has a 01_cutoff_roc_filter folder 
   - this folder contains a huge excel table with two sheets 
   - sheet "ratio_normSI_anot" holds all ratios, the normalized SIs from the previous step, detailed annotation and pass cutoff columns 
   - sheet "cutoff_plots_data" contains normalized ratios (median of FP proteins), log2 of that normalized ratio as well as TPR-FPR columns - this data corresponds to plots in folder cutof_ROC_plots, based on this data we decide if a Protein passes or not (highest TPR-FPR value is our treshold in the log2 ratio column), this sheet will be of help if anyone is interested to see which other proteins were very close to the treshold but didn't pass etc. 
   - the folder called tp_fp_plots contains visualization of number of proteins that pass and do not pass cutoff and if they are in the TP or FP or no group
   - a second filtering step determines our final protein list which is labelled as final_protein_list_xx
   - this final list contains proteins that pass cutoffs with at least 2 replicates per condition in at least one condition of this file (see filtering in 01_cutoff_roc_filter, excel table, sheet "ratio_normSI_annot" - pass_cutoff_result column)
   - based on this final list volcano plots are created, volcano plot data can be found in the final table
   - and you can also find a barplot indicating the abundance of TP and FP proteins in that final list compared to a barplot without any filtering (ie with what we started initially)
2. Serum samples
   - also here each sample has its own folder 
   - but with serum samples we do not perform any cutoff roc analysis
   - TODO: before providing the final table for serum, I could filter if based on enrichment in min 2 replicates, for now there is no filtering step  
   - the final table is annotated and volcano plots are created based on that table, volcano plot data can be found in the final table

Ad volcano plots: p-values were calculated with T-test for the means of two independent samples and volcano plots show uncorrected p-values. Significant proteins which show a > 1.5 fold change are highlighted in either red or green. 