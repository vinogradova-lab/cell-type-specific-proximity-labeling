# Cell type-specific proximity labeling analysis

Analysis code accompanying the 2026 manuscript

## Repository Structure

| Folder | Description | Author |
|---|---|---|
| `scripts/` | Processing of raw mass spectrometry searches | [Nathalie Ropek](nropek@protonmail.com) |
| `notebooks/` | Figure generation and post-processing | [Henry Sanford](hsanford@rockefeller.edu) |
| `templates/` | Examples of formatted metadata for mass spec processing | [Nathalie Ropek](nropek@protonmail.com) |
| `data/` | Reference datasets for annotating protein secretion | [Nathalie Ropek](nropek@protonmail.com) |

## Guide to processing raw mass spectrometry searches

### Step 1 (mouse_fasta_to_FP_TP_tables_without_signalp.py) - for reference in case fasta file need to be modified
- create main table containing protein metadata (based off of our mouse fasta file used for proteomics search)
  - read in fasta file to get all uniprot IDs
  - use uniprot.org to get all annotations for these proteins
  - label proteins as TP based on following criteria in subcellular location column: 
    - "Endoplasmic Reticulum", "Secreted", "secreted", "Endoplasmic reticulum", "Rough endoplasmic reticulum", "endoplasmic reticulum"
  - label proteins as FP if uniprot id is found in human mitomatrix list from supplementary table S1 [Rhee et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3916822/)
    - obtain human-mouse gene orthologs by following tutorial: [How to get all the orthologous genes between two species](https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/) 
    - this allows us to get the correct mouse gene names which correspond to the mitomatrix protein list, we then crossreference our fasta table and label them as FP
  - additionally make sure that if there is a SignalP annotation (ie if SignalP column is not empty) for a protein marked as FP then remove the FP annotation, also if protein is marked as both TP and FP, keep it in the TP list and remove it from FP list
  - crossreference table with secretion prediction resource
  - crossreference main list with spleen, adipose tissue etc.
  - crossreference table with Crapome database - cutoff > 200
  - import NCBI gene annotations (using goatools package) and crossreference gene symbols to add gene_ids into our final list (for go term analysis)
  - save table


### Step 2 (turboid_analysis_pipeline.py)

---------------------------------------------------------------------------------------------------------------------------------------------

Mac with Docker installed: 
- go to data turboid folder (where start_turbo_id_downstream.sh is located)
- `docker build -t turboid_downstream .` 
- `./start_turbo_id_downstream.sh process_data`

Then: 
- follow instructions in templates/metadata_turbo_id.xlsx, Analysis details tab
- In your 01_processed_files/conditions_metadata.csv, order your conditoins so that the condition you want to measure regulation against should be at the 
start of your condition list
---------------------------------------------------------------------------------------------------------------------------------------------
  
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
 
More details: GOEA analysis was performed on the significant up or down or up & down regulated proteins (mapped to gene IDs). Enriched ontologies with a Benjamini-Hochberg-corrected p-value less than 0.05 were retained and top 10 terms were plotted. 

**Volcano plots**

For every protein the log2 fold-change of the ratio between condtions was calculated. p-values were calculated with T-test for the means of two independent samples for all replicates and volcano plots show uncorrected p-values. Significant proteins which show a > 1.5 fold change are highlighted in either red or green. 

**Heatmap**

A heatmap of top proteins is provided as pdf (can be searched cmd+F) and can be compared to volcano plots

### Details of ROC cutoff analysis 
Following the approach of: [Cho et al. 2020](https://www.nature.com/articles/s41596-020-0399-0)

1. Calculate Ratios per condition (example of one condition) 
   - Each condition within each file was analyzed on its own, therefore each condition will have its own cutoff in the end. 
   - Example ratios of file CRW_B17A for condition iWAT_HFD
        - ratio_iWAT_HFD_cre+_1/iWAT_HFD_cre-1
        - ratio_iWAT_HFD_cre+_1/iWAT_HFD_cre-_2 	
        - ratio_iWAT_HFD_cre+_2/iWAT_HFD_cre-_1 	
        - ratio_iWAT_HFD_cre+_2/iWAT_HFD_cre-_2 	
        - ratio_iWAT_HFD_cre+_3/iWAT_HFD_cre-_1 	
        - ratio_iWAT_HFD_cre+_3/iWAT_HFD_cre-_2 
   - Each replicate will be treated on its own (all steps according to Cho et al - Nature Protocol)
   - Example: ratio_ iWAT_HFD_cre+_1/iWAT_HFD_cre-1
   - Normalize ratio_ iWAT_HFD_cre+_1/iWAT_HFD_cre-1 by the median ratio of proteins that are labelled as FP
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
     - iWAT_HFD_cre+_1, iWAT_HFD_cre+_2 etc.: normalized signal intensity values
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


