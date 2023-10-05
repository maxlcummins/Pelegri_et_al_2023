library(circlize)
library(pheatmap)
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

#change treatment here
N2 <- N2 %>% mutate(treatment = 'N2')

#change treatment here
mat <- N2 %>% filter(category == 'Neural') %>% select(treatment, `term description`, `Gene Name`) %>% group_by(`Gene Name`, `term description`) %>% summarise(count = n())

#up_genes <- C1 %>% filter(category == 'Neural' & `RS2_3D_C1_vs_RS2_3D_DMEM_log2 fold change` > 0) %>% pull(`Gene Name`)

#down_genes <- C1 %>% filter(category == 'Neural' & `RS2_3D_C1_vs_RS2_3D_DMEM_log2 fold change` < 0) %>% pull(`Gene Name`)

circos.clear()

#Change number to change colours
set.seed(1)

#group1 <- rep("A", length(up_genes))
#names(group1) <- up_genes

#group2 <- rep("B", length(down_genes))
#names(group2) <- down_genes

test = chordDiagram(mat, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        #font size for labels = cex
        circos.text(mean(xlim), ylim[1] + cm_h(1.5), sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex =0.7)
        #font size for ticks = cex
        circos.axis(h = "top", labels.cex =0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)

mat2 <- N2 %>% filter(category == 'Neural') %>% select(treatment, `term description`, `Gene Name`, `RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`)

mat2 <- mat2 %>% as.data.frame() %>% select(`Gene Name`, `RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`)

mat2 <- mat %>% inner_join(mat2, by = "Gene Name") %>% mutate(`RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change` = as.numeric(`RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`))

mat2 <- mat2 %>% unique()

mat2[nrow(mat2)+1,] <- NA
mat2[nrow(mat2)+1,] <- NA

mat2[nrow(mat2),4] <- 4.05
mat2[nrow(mat2)-1,4] <- -4.05

#Normalized Data
mat2$normalized = (mat2$`RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`-min(mat2$`RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`, na.rm = T))/(max(mat2$`RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`, na.rm = T)-min(mat2$`RS3_3D_N2_vs_RS3_3D_DMEM_log2 fold change`, na.rm = T))

mat2 <- mat2 %>% mutate(normalized = replace_na(normalized, 0.5))

f <- colorRamp(colors = c("red", '#FFBFBF', "#BFBFFF", "blue"))

values <- mat2$normalized

mat2$col <- rgb(f(values)/255)

grid_cols <- mat2 %>% select(`Gene Name`)

grid_cols <- unique(grid_cols)

grid_cols <- left_join(grid_cols, mat2, by = "Gene Name") %>% select(`Gene Name`, col)

grid_cols <- grid_cols %>% unique()

grid.cols <- grid_cols$col

names(grid.cols) <- grid_cols$`Gene Name`

circos.clear()

mat2 <- mat2 %>% arrange(desc(normalized))

mat3 <- mat2 %>% select(all_of(colnames(mat)))

set.seed(1)

pdf("Treatments_3D_DMEM_3D/new_figures/N2_chord.pdf") 

chordDiagram(mat, annotationTrack = "grid", preAllocateTracks = 1, grid.col = c("Brain" = "#913640", "Nervous system" = '#E8C167', 'Axon guidance' = '#D67500', 'Brain cell line' = '#FCFFC9', 'Forebrain' = '#1D0B14', grid.cols), order = c(mat3$`Gene Name`, unique(mat3$`term description`)))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        #font size for labels = cex
        circos.text(mean(xlim), ylim[1] + cm_h(1.5), sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex =0.7)
        #font size for ticks = cex
        circos.axis(h = "top", labels.cex =0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)


dev.off()

pdf("Treatments_3D_DMEM_3D/new_figures/chord_legend.pdf") 

pheatmap(seq(4, -4, by = -0.25), color = colorRampPalette(c("red", '#FFBFBF', "#BFBFFF", "blue"))(100), cluster_rows = F, cluster_cols = F, display_numbers = T)

dev.off()
