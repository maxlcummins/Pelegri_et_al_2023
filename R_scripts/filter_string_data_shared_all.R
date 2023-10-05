library(readxl)
library(readr)
library(tidyverse)

#Change here
treatment_name <- "Shared_all"

#Change here
#Define string file path
string_file_path <- "String/Shared_All/Shared_all_enrichment.all.tsv"

#Read in string data
string_data <- read_delim(string_file_path, delim = "\t")

#Define proteomic file path
file_path <- "LFQ analyst Results/RS2_RS3_ALL_1.5FOLD.csv"

#Read in proteomic data
df <- read_csv(file_path)

#Define list of descriptors to search for, pipe delimited
list_of_descriptions <- 'brain|nervous|synap|axon|neural|cerebral|ganglion|temporal lobe|synap|neuro|filopodium|amygdala|dendrit'

#Filter the string data to only contin things from our list of descriptions in the term description column
#THIS WILL ALLOW PARTIAL MATCHES AND MIGHT BE CASE SENSITIVE
neural <- string_data %>% filter(str_detect(tolower(`term description`), list_of_descriptions))

#Make a new dataframe for neural hits
neural <- neural %>%
  #Select columns of intereset
  select(`term ID`, `term description`, `matching proteins in your network (labels)`) %>%
  #Split comma delimited column into rows
  mutate(`Gene Name` = strsplit(as.character(`matching proteins in your network (labels)`), ",")) %>% 
  unnest(`Gene Name`) %>%
  #Remove old column with comma delimited protein names
  select(-`matching proteins in your network (labels)`) %>%
  #Join to our proteomic dataframe with fold changes etc
  left_join(df, by  = "Gene Name")

#Change here
#Select columns of interest - note we need to change to the appropriate treatment's fold changes and significant columns
neural <- neural %>%
  select(`term ID`, `term description`, `Gene Name`, `RS2_3D_B27_vs_RS2_3D_DMEM_log2 fold change`, `RS2_3D_B27_vs_RS2_3D_DMEM_significant`, `RS2_3D_C1_vs_RS2_3D_DMEM_log2 fold change`, `RS2_3D_C1_vs_RS2_3D_DMEM_significant`, `RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`, `RS3_3D_N2_vs_RS3_3D_DMEM_significant`) %>%
  mutate(category = "Neural")




list_of_descriptions <- 'actin'

actin <- string_data %>% filter(str_detect(tolower(`term description`), list_of_descriptions))

actin <- actin %>%
  select(`term ID`, `term description`, `matching proteins in your network (labels)`) %>%
  mutate(`Gene Name` = strsplit(as.character(`matching proteins in your network (labels)`), ",")) %>% 
  unnest(`Gene Name`) %>%
  select(-`matching proteins in your network (labels)`) %>%
  left_join(df, by  = "Gene Name")

#Change here
actin <- actin %>%
  select(`term ID`, `term description`, `Gene Name`, `RS2_3D_B27_vs_RS2_3D_DMEM_log2 fold change`, `RS2_3D_B27_vs_RS2_3D_DMEM_significant`, `RS2_3D_C1_vs_RS2_3D_DMEM_log2 fold change`, `RS2_3D_C1_vs_RS2_3D_DMEM_significant`, `RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`, `RS3_3D_N2_vs_RS3_3D_DMEM_significant`) %>%
  mutate(category = "Actin")



list_of_descriptions <- 'ribo'

ribosomal <- string_data %>% filter(str_detect(tolower(`term description`), list_of_descriptions))

ribosomal <- ribosomal %>%
  select(`term ID`, `term description`, `matching proteins in your network (labels)`) %>%
  mutate(`Gene Name` = strsplit(as.character(`matching proteins in your network (labels)`), ",")) %>% 
  unnest(`Gene Name`) %>%
  select(-`matching proteins in your network (labels)`) %>%
  left_join(df, by  = "Gene Name")


#Change here
ribosomal <- ribosomal %>%
  select(`term ID`, `term description`, `Gene Name`, `RS2_3D_B27_vs_RS2_3D_DMEM_log2 fold change`, `RS2_3D_B27_vs_RS2_3D_DMEM_significant`, `RS2_3D_C1_vs_RS2_3D_DMEM_log2 fold change`, `RS2_3D_C1_vs_RS2_3D_DMEM_significant`, `RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`, `RS3_3D_N2_vs_RS3_3D_DMEM_significant`) %>%
  mutate(category = "Ribosomal")

#Reorder columns
out_df <- bind_rows(ribosomal, neural, actin) %>% select(category, everything())


#Generalised for other scripts
write_delim(x = out_df , file = paste0("String/", treatment_name, "/", treatment_name, "_string_categorised.tsv"), delim = "\t")


#Heatmap
#out_df %>% filter(category == "Neural") %>% select(`Gene Name`, `RS2_3D_B27_vs_RS2_2D_DMEM_log2 fold change`, `RS2_3D_C1_vs_RS2_2D_DMEM_log2 fold change`, `RS3_3D_N2_vs_RS3_2D_DMEM_log2 fold change`) %>% unique() %>% column_to_rownames('Gene Name') %>% mutate(across(everything(), ~ replace_na(., 0))) %>% pheatmap(fontsize_row = 5, fontsize_col = 8)