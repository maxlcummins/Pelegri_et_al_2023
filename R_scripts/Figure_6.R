library(circlize)

library(ggplot2)
library(ComplexUpset)
library(tidyverse)
library(readxl)

B27 <- read_delim("Treatments_3D_DMEM_3D/String/B27/B27_string_categorised.tsv", delim = "\t")
B27_C1 <- read_delim("Treatments_3D_DMEM_3D/String/B27_C1/B27_C1_string_categorised.tsv", delim = "\t")
B27_N2 <- read_delim("Treatments_3D_DMEM_3D/String/B27_N2/B27_N2_string_categorised.tsv", delim = "\t")
C1 <- read_delim("Treatments_3D_DMEM_3D/String/C1/C1_string_categorised.tsv", delim = "\t")
C1_N2 <- read_delim("Treatments_3D_DMEM_3D/String/C1_N2/C1_N2_string_categorised.tsv", delim = "\t")
N2 <- read_delim("Treatments_3D_DMEM_3D/String/N2/N2_string_categorised.tsv", delim = "\t")
shared_all <- read_delim("Treatments_3D_DMEM_3D/String/Shared_All/Shared_all_string_categorised.tsv", delim = "\t")

B27 <- B27 %>% mutate(treatment = 'B27') %>% rename(fold_change = `RS2_3D_B27_vs_RS2_3D_DMEM_log2 fold change`)
C1 <- C1 %>% mutate(treatment = 'C1') %>% rename(fold_change = `RS2_3D_C1_vs_RS2_3D_DMEM_log2 fold change`)
N2 <- N2 %>% mutate(treatment = 'N2') %>% rename(fold_change = `RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`)
B27_C1 <- B27_C1 %>% mutate(treatment = 'B27_C1') %>% rename(fold_change = `RS2_3D_B27_vs_RS2_3D_DMEM_log2 fold change`)

all <- bind_rows(B27, C1) %>% bind_rows(N2) %>% bind_rows(B27_C1)  %>% filter(fold_change > 0)

mat <- all %>% filter(category == 'Neural') %>% select(treatment, `term description`, `Gene Name`) %>% group_by(treatment, `Gene Name`, `term description`) %>% summarise(count = n())

mat1 <- mat %>% as.data.frame() %>% select(-`Gene Name`)

mat2 <- mat %>% as.data.frame() %>% select(-`term description`)

circos.clear()

#Change number to change colours
set.seed(4)

chordDiagram(mat1, annotationTrack = "grid", preAllocateTracks = 1, grid.col = c("Brain" = "#913640", "Nervous system" = '#E8C167', 'Axon guidance' = '#D67500', 'Brain cell line' = '#FCFFC9', 'Forebrain' = '#1D0B14'))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1] + cm_h(1.5), sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex =1)
        circos.axis(h = "top", labels.cex =0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

#Down regulated genes now

B27 <- read_delim("Treatments_3D_DMEM_3D/String/B27/B27_string_categorised.tsv", delim = "\t")
B27_C1 <- read_delim("Treatments_3D_DMEM_3D/String/B27_C1/B27_C1_string_categorised.tsv", delim = "\t")
B27_N2 <- read_delim("Treatments_3D_DMEM_3D/String/B27_N2/B27_N2_string_categorised.tsv", delim = "\t")
C1 <- read_delim("Treatments_3D_DMEM_3D/String/C1/C1_string_categorised.tsv", delim = "\t")
C1_N2 <- read_delim("Treatments_3D_DMEM_3D/String/C1_N2/C1_N2_string_categorised.tsv", delim = "\t")
N2 <- read_delim("Treatments_3D_DMEM_3D/String/N2/N2_string_categorised.tsv", delim = "\t")
shared_all <- read_delim("Treatments_3D_DMEM_3D/String/Shared_All/Shared_all_string_categorised.tsv", delim = "\t")

B27 <- B27 %>% mutate(treatment = 'B27') %>% rename(fold_change = `RS2_3D_B27_vs_RS2_3D_DMEM_log2 fold change`)
C1 <- C1 %>% mutate(treatment = 'C1') %>% rename(fold_change = `RS2_3D_C1_vs_RS2_3D_DMEM_log2 fold change`)
N2 <- N2 %>% mutate(treatment = 'N2') %>% rename(fold_change = `RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`)
B27_C1 <- B27_C1 %>% mutate(treatment = 'B27_C1') %>% rename(fold_change = `RS2_3D_B27_vs_RS2_3D_DMEM_log2 fold change`)

all <- bind_rows(B27, C1) %>% bind_rows(N2) %>% bind_rows(B27_C1)  %>% filter(fold_change < 0)

mat <- all %>% filter(category == 'Neural') %>% select(treatment, `term description`, `Gene Name`) %>% group_by(treatment, `Gene Name`, `term description`) %>% summarise(count = n())

mat1 <- mat %>% as.data.frame() %>% select(-`Gene Name`)

mat2 <- mat %>% as.data.frame() %>% select(-`term description`)

circos.clear()

#Change number to change colours
set.seed(4)

chordDiagram(mat1, annotationTrack = "grid", preAllocateTracks = 1, grid.col = c("Brain" = "#913640", "Nervous system" = '#E8C167', 'Axon guidance' = '#D67500', 'Brain cell line' = '#FCFFC9', 'Forebrain' = '#1D0B14'))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        circos.text(mean(xlim), ylim[1] + cm_h(1.5), sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex =1)
        circos.axis(h = "top", labels.cex =0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)