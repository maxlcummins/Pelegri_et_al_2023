# Pelegri et al, 2023
A repository containing supplementary scripts for data processing and visualisation in Pelegri et al, 2023

## Workflow overview

1. LFQ analyst

After running LFQ analyst, the resulting CSV file is processed using the script `scripts/proteomics_workflow_1.5_fold.Rmd`. This script determines which proteins were significantly over or under expressed within each treatments, determines their intersection between treatments. It writes to file the genes i) unique to each treatment, ii) found in pairwise combinations of two treatments and iii) found in all treatments. It also generates the venn diagriams visualised in Figure 5.

`LFQ analyst Results` contains two files - one with edited gene names (flagged as 'edited') and the original. Some gene names were manually changed to gene alisases due to differences in the name outputs between String and LFQ Analyst to simplify the relinkage of fold-change values back to the String output. Additionally, LFQ Analyst erroneously reformats some gene names to date format, for example SEPT7 is changed to Sep-7. This prevents analysis in String as this is an invalid alias for this gene. This and other examples of SEPT genes in our LFQ Analyst output were also manually corrected.

2. String

Outputs from step 1 are run in String (see paper methods) and functional associations of proteins are assigned. We next filter this string data to include only proteins associated with pathways of interest (e.g. Neural pathways) using scripts titled `filter_string_data_N2.R`, etc. Run order for these files doesn't matter. This script then binds the fold changes for each protein back in as well.

3. Data visualisation

Next we run `Figure_6.R` and `Figure_7_N2.R`, `Figure_7_C1.R`, and `Figure_7_B27.R` to generate our chord diagrams.
