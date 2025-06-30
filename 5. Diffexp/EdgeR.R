## load the libraries here for the analyses
library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(ggrepel)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(EnhancedVolcano)
library(forcats)

## read the inputs here
countfile     = "combined_counts.tsv" # this is the gene to quantification profile for the samples that we generated in Step/Folder 4
control       = "P-MyC-CaP" # this tells the reference condition for the analyses
treatment     = "SMC-NEPT" # this tells the condition that is to be compared to the reference in the analyses
control_rep   = "3" # this tells the number of samples that belong to the reference condition
treatment_rep = "3" # this tells the number of samples that belong to the contrasting/treatment condition
path          = "mart_export.txt" # this is the annotation file that maps the gene-identifiers to gene-symbols that are readily recognized, 
                                         #and adds a few columns that says the gene function, location in the genome etc

## define the output prefixes here, this uses the control and treatment variables that we define to create filenames that are easy to identify
outprefix <- paste(treatment, control, sep="_vs_")

######
# Read the files and prepare for analyses
######
#setwd("/depot/tlratlif/data/2024_4_Gada_SMC_RNA_Bulk/30-994315926/P1/output/Diffexp/")
Anno <- read.table(path, sep = "\t", stringsAsFactors = F, header = TRUE, quote = "")

#setwd("/depot/tlratlif/data/2024_4_Gada_SMC_RNA_Bulk/30-994315926/P1/output/counts/")

raw.data = read.table(countfile, row.names="Gene_ID", header = T,sep = "\t", as.is = T, check.names = F, quote = "")

group <- factor(rep(c(control, treatment),
                    times=c(control_rep, treatment_rep)),
                levels = c(control, treatment))

control_mat<- select(raw.data,starts_with(control))
treatment_mat<- select(raw.data,starts_with(treatment))
raw.counts<- merge(control_mat,treatment_mat,by="row.names",all.x=TRUE)

raw.counts<- tibble::column_to_rownames(raw.counts, var="Row.names")


###################################################################################
#                          edgeR Analysis
###################################################################################
#### start using edgeR here to get DE. 
# start by prepping an assay object that edgeR can use and later export out the CPM matrices
y <- DGEList(counts=as.matrix(raw.counts), group=group)
x <- calcNormFactors(y)
xlog <- cpm(x, log=T)
write.table(xlog,file = paste0(outprefix,".edgeR_logCPM.txt"), sep = "\t", quote = F, row.names = T)
xcpm <- cpm(x, log=F)
write.table(xcpm,file = paste0(outprefix,".edgeR_CPM.txt"), sep = "\t", quote = F, row.names = T)


### filter to remove the low expressed genes. filtering keeps genes that have count-per-million (CPM) above a minimum number of samples 
#in the smallest group or the group that is specified. by default, we pick the smallest group. and then we normalize the data
keep <- filterByExpr(y)
filtered.data <- y[keep,keep.lib.sizes=FALSE]
y <- calcNormFactors(filtered.data)



### create a design matrix without the intercept, edgeR basically gives the same results, this just makes the design matrix look simpler to troubleshoot.
# rows in the designmatrix are samples, and the columns are conditions
design <- model.matrix(~0+group)
colnames(design)<- levels(y$samples$group)


###################################################################################
#                          PCAPlot
###################################################################################
# We use the filtered gene-counts to generate the pca
# get logCPMs again for the filtered genes and calculate row-wise variances, we then select the top 500 genes based on the largest variations that occur
logCPMs<- cpm(y, log=T)
rowvar <- apply(logCPMs, 1, var)
ordered_genes <- order(rowvar, decreasing=TRUE)
top5000 <- head(ordered_genes, 5000)
logCPM_top5000 <- logCPMs[top5000,]

# use the base-R prcomp function to get the principal components, and calculate the percent variations, we then select the top 2 PCs to plot

pca <- prcomp(t(logCPM_top5000))
to_plot <- data.frame(pca$x, y$samples)
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
use.pcs <- c(1,2)

group_colors <- c("pink", "royalblue")
to_plot_group_color <- unique(to_plot$group)
color_mapping <- setNames(group_colors, to_plot_group_color)


p1 <- ggplot(to_plot, aes(PC1, PC2,color= group,
                        theme(legend.position = "none"),
                        label = rownames(to_plot))) +
  geom_point(size=1) + scale_color_manual(values = color_mapping) +
  xlab(paste0("PC1: ",round(percentVar[1],2),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2],2),"% variance"))

p2<- p1 + geom_point(size = 2) + geom_text_repel(hjust = 0.1, vjust = -0.1, nudge_x = 0.1, size = 4) +
  theme(legend.position = "none",
        axis.title=element_text(size=14))

png(file=paste0(outprefix,".PCA_plot.png") ,width = 2000, height = 2000, bg = "white", pointsize = 10, res=300)
p2
dev.off()
###################################################################################
#                          Expression Analyses
###################################################################################

### we then estimate the dispersions, which initially figures out the biological variation between the samples, 
#then fits a dispersion which quantifies the degree of variability relative to a poisson distribution
# we then do a Generalized Linear Models with Quasi-likelihood Tests to figure out the DEs and then annotate the gene-lists

y <- estimateDisp(filtered.data,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,contrast=c(-1,1))


edgeR_DE = as.data.frame(topTags(qlf, sort.by = "PValue", n = Inf))
edgeR_DE = tibble::rownames_to_column(edgeR_DE, "Gene_ID")
names(Anno)<-c("Gene_ID","Gene.stable.ID.version","Transcript.stable.ID","Transcript.stable.ID.version","Gene.description","Gene.name","Gene.type")
edgeR_DE = left_join(edgeR_DE, Anno, by="Gene_ID")
edgeR_DE <- edgeR_DE[!duplicated(edgeR_DE$Gene_ID), ]

##To export edgeR_DE into excel sheet

library(xlsx)
write.xlsx(edgeR_DE, "/depot/tlratlif/data/2024_4_Gada_SMC_RNA_Bulk/30-994315926/P1/output/Diffexp/edgeR_DE.xlsx", sheetName = "Edger_Top5000", col.names = TRUE, row.names = TRUE)
#saveRDS(edgeR_DE, file = paste0(outprefix,".edgeR_DE.rds"))

significant_genes <- subset(edgeR_DE, FDR < 0.05)
write.xlsx(significant_genes, "/depot/tlratlif/data/2024_4_Gada_SMC_RNA_Bulk/30-994315926/P1/output/Diffexp/edgeR_DE_significant_genes-fdr-5percent.xlsx", sheetName = "Edger_FDR_5P", col.names = TRUE, row.names = TRUE)
#saveRDS(edgeR_DE, file = paste0(outprefix,".edgeR_DE.rds"))
###################################################################################
#                          VolcanoPlot
###################################################################################
### we use the logFC metric to determine and color the direction of regulation. Significant genes are then determined 
#if the FDR threshold of 0.05 is being met or not and colored accordingly

edgeR_DE.show = edgeR_DE %>%
  dplyr::select(Gene.name, logFC, PValue, FDR)

# add color to logFC direction mapping and assign what genes are up, neutral or down

keyvals <- ifelse(
  edgeR_DE.show$logFC < 0 & edgeR_DE.show$FDR < 0.05, 'royalblue',
  ifelse(edgeR_DE.show$logFC > 0 & edgeR_DE.show$FDR < 0.05, 'red',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'Down-regulated'
names(keyvals)[keyvals == 'grey'] <- 'No change'
names(keyvals)[keyvals == 'red'] <- 'Up-regulated'

# plot here



volc<-EnhancedVolcano(edgeR_DE.show,
                      lab = edgeR_DE.show$Gene.name,
                      x = 'logFC',
                      y = 'FDR', legendLabels=c('Up-regulated','No change','Down-regulated'),
                      legendPosition = 'right',
                      legendLabSize = 16,
                      legendIconSize = 5.0, drawConnectors = TRUE,
                      widthConnectors = 0.75, colCustom = keyvals, title = paste0("Volcano plot of the comparison: ", outprefix), subtitle="", cutoffLineType = 'blank', FCcutoff = 4.0, ylab = bquote(~-Log[10]~italic(FDR)))

png(file=paste0(outprefix,".Volcano_plot.png") ,width = 3000, height = 3000, bg = "white", pointsize = 1, res=300)
volc
dev.off()
###################################################################################
#                          Volcano plot for highlighting genes of interest
###################################################################################

# Define the significance thresholds
log2fc_threshold <- 1.0
fdr_threshold <- 0.05

# Add a column to your data frame to categorize each gene based on the criteria
edgeR_DE.show$color <- ifelse(edgeR_DE.show$FDR < fdr_threshold & abs(edgeR_DE.show$logFC) > log2fc_threshold, 'red',ifelse(edgeR_DE.show$FDR < fdr_threshold, 'royalblue',ifelse(abs(edgeR_DE.show$logFC) > log2fc_threshold, 'green', 'grey')))

# Plot the volcano plot
EnhancedVolcano(edgeR_DE.show, lab =edgeR_DE.show$Gene.name,
                selectLab = c("Krt8","Notch1","Ar","Psca","Tmprss2","Etv3","Etv4","Etv5","Dht","Folh1","Rb1","Pten","Ascl2","Foxa1","Cxcl15","Cxcl5","Cxcl1","Cxcl14","Cxcl16","Cxcl12","Mapk1","Chga","Chgb","Syp","Ncam1","Ncam2","Notch2","Notch3","Notch4","Twist2","Vim","Snail1","Mycn","Onecut2","Akt1","Akt2","Akt3"),
                x = 'logFC',
                y = 'FDR',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'adjusted p-value'),
                pCutoff = fdr_threshold,
                FCcutoff = log2fc_threshold,
                pointSize = 3.0,
                #legendLabels = c("Not Significant", "Log2FC Significant", "FDR Significant", "Both Significant"),  # Customize legend labels
)


## Highlighting the genes of interest with threshold 

# Assuming edgeR_DE.show is your data frame
#highlight_rows <- c(50, 100)  # Adjust based on your specific row numbers

# Create a custom color vector
#keyvals <- ifelse(
#  edgeR_DE.show$logFC <= -0.4 & edgeR_DE.show$FDR <-0.3, 'royalblue',
#  ifelse(edgeR_DE.show$logFC > 0 & edgeR_DE.show$FDR < 0.05, 'red', 
#     ifelse(edgeR_DE.show$logFC < 4 & edgeR_DE.show$FDR > 0.05, 'forestgreen',
#     ifelse(edgeR_DE.show$logFC < 0 & edgeR_DE.show$FDR > 0.05, 'forestgreen', 'grey')))
#)




#keyvals[is.na(keyvals)] <- 'grey'
#names(keyvals)[keyvals == 'royalblue'] <- 'Down-regulated'
#names(keyvals)[keyvals == 'grey'] <- 'No change'
#names(keyvals)[keyvals == 'red'] <- 'Up-regulated'
#names(keyvals)[keyvals == 'forestgreen'] <- 'Log2fc'

#volc <- EnhancedVolcano(edgeR_DE.show,
#                        lab = ifelse(1:nrow(edgeR_DE.show) %in% highlight_rows, edgeR_DE.show$Gene.name, ""),
#                        x = 'logFC',
#                        y = 'FDR',
#                        legendLabels = c('Up-regulated', 'No change', 'Down-regulated','Log2fc'),
#                        legendPosition = 'right',
#                        legendLabSize = 16,
#                        legendIconSize = 5.0,
#                        drawConnectors = TRUE,
#                        widthConnectors = 0.75,
#                        colCustom = keyvals,
#                        title = paste0("Volcano plot of the comparison: ", outprefix),
#                        subtitle = "",
#                       pCutoff = 0.03,   # Adjusted p-value cutoff
#                        FCcutoff = 0.4,   # Adjusted fold change cutoff
#                        cutoffLineType = 'solid',
#                        cutoffLineWidth = 0.5,
#                        cutoffLineCol = 'black',
#                        ylab = bquote(~-Log[10]~italic(FDR)))

# Print the plot
#print(volc)

###################################################################################
#                          Heatmaps
###################################################################################
## get the CPM information first

CPM <- cpm(y, log=F)
CPM <- tibble::rownames_to_column(data.frame(CPM),"Gene_ID")

## subset to get the gene_id to symbol mapping
gene_sym <- Anno %>% select(Gene_ID, Gene.name)
head(gene_sym)
## This is a quick way to select the top 10 and bottom 10 genes. This ensures that we have both Up and Down regulated genes.
gene_order <- select(edgeR_DE, "Gene_ID", "logFC", "PValue","Gene.name")
gene_order$Rank <- sign(gene_order$logFC) * -log10(gene_order$PValue)
gene_order_export <- select(gene_order, "Gene_ID", "Rank","Gene.name")
gene_order_export <- dplyr::arrange(gene_order_export, desc(Rank)) %>% drop_na()
gene_order_export <- select(gene_order_export, "Gene_ID", "Rank","Gene.name")
top_60_genes <- rbind(head(gene_order_export,30), tail(gene_order_export,30)) 


top_60_cpms <- inner_join(top_60_genes, CPM, by="Gene_ID")
#top_20_cpms <- inner_join(gene_sym, top_20_cpms, by="Gene_ID")
top_60_cpms <- top_60_cpms %>% select(-Gene_ID,-Rank) %>% tibble::column_to_rownames(., "Gene.name")
#top_20_cpms <- top_20_cpms[!duplicated(top_20_cpms$Gene_ID), ]


png(file=paste0(outprefix,".Heatmap.png") ,width = 3000, height = 3000, bg = "white", pointsize = 1, res=300)
draw(Heatmap(t(scale(t(top_60_cpms))), name = "CPM", show_column_dend=F, show_row_dend=T,
cluster_columns=F, cluster_rows=T, row_names_gp=gpar(fontsize=8), row_title_rot=0,
column_split=group, heatmap_legend_param = list(direction = "horizontal"), border = TRUE),
heatmap_legend_side = "bottom", annotation_legend_side = "bottom",  merge_legend = TRUE)
dev.off()


## To pull out ceratin genes for heatmap

library(circlize)

# Define your list of genes

##Ar related genes
genes_list1 <- c("Ar", "Psca", "Tmprss6","Tmprss9","Tmprss11e","Etv3","Etv4","Etv5","Dhtkd1","Folh1") 


##Lineage plasticity related genes
genes_list2 <- c("Rb1", "Pten", "Ascl2","Foxa1","Cxcl15","Cxcl5","Cxcl1","Cxcl114","Cxcl16","Cxcl12","Mapk1","Mapk14","Mapk10","Mapk11","Mapk15","Mapk12","Mapk13") 



## neuroendocrine differentaiation genes'
genes_list3 <- c("Chga","Chgb","Syp","Ncam1","Ncam2")



## Emt changes
genes_list4 <- c("Krt8","Notch1","Notch2","Notch3","Twist2","Vim","Snail1")



## Myc family
genes_list5 <- c("Mycn","Onecut2","Akt1","Akt2","Akt3")

# Subset the CPM matrix to include only the genes in your list

gene_order <- select(edgeR_DE, "Gene_ID", "logFC", "PValue","Gene.name")
gene_order$Rank <- sign(gene_order$logFC) * -log10(gene_order$PValue)
gene_order_export <- select(gene_order, "Gene_ID", "Rank","Gene.name")
gene_order_export <- dplyr::arrange(gene_order_export, desc(Rank)) %>% drop_na()
gene_order_export <- select(gene_order_export, "Gene_ID", "Rank","Gene.name")
View(gene_order_export)
top_genes <- gene_order_export
top_cpms <- inner_join(top_genes, CPM, by="Gene_ID")
subset_cpms <- top_cpms[rownames(top_cpms) %in% genes_list5, ]
View(top_cpms)
colnames(top_cpms)
[1] "Gene_ID"     "Rank"        "Gene.name"   "P.MyC.CaP.1" "P.MyC.CaP.2" "P.MyC.CaP.3"
[7] "SMC.NEPT.1"  "SMC.NEPT.2"  "SMC.NEPT.3" 
filtered_top_cpms <- top_cpms %>% filter(Gene.name %in% genes_list5)
filtered_top_cpms
Gene_ID      Rank Gene.name P.MyC.CaP.1 P.MyC.CaP.2 P.MyC.CaP.3 SMC.NEPT.1
1 ENSMUSG00000046532  5.149788        Ar   1826.6951   1844.8252   1919.0093  2535.1142
2 ENSMUSG00000030541 -6.383452      Idh2    436.5612    434.3761    421.1909   315.6164
SMC.NEPT.2 SMC.NEPT.3
1  2327.5614  2321.4906
2   324.2846   328.0963
filtered_top_cpms <- filtered_top_cpms %>% select(-Gene_ID,-Rank) %>% tibble::column_to_rownames(., "Gene.name")


# Draw the heatmap using the subsetted CPM matrix


draw(Heatmap(t(scale(t(filtered_top_cpms))), name = "CPM", show_column_dend = FALSE, show_row_dend = TRUE,
             cluster_columns = FALSE, cluster_rows = TRUE, row_names_gp = gpar(fontsize = 8), row_title_rot = 0,
             column_split = group, heatmap_legend_param = list(direction = "horizontal"), border = TRUE),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)

###################################################################################
#                          Filtering
###################################################################################
## we export the genes that pass a 5% FDR threshold and/or various logFC cutoffs
# EdgeR - FDR filtered
edgeR_DE_FDR5P = dplyr::filter(edgeR_DE, FDR < 0.05)

edgeR_DE_FDR5P_logFC1 = dplyr::filter(edgeR_DE, FDR < 0.05 &  (logFC < -1 | logFC > 1)  )

edgeR_DE_FDR5P_logFC2 = dplyr::filter(edgeR_DE, FDR < 0.05 &  (logFC < -2 | logFC > 2)  )

cutoff_5P = c(nrow(edgeR_DE_FDR5P), nrow(edgeR_DE_FDR5P_logFC1), nrow(edgeR_DE_FDR5P_logFC2))

DE_stats = data.frame(cutoff_5P)
rownames(DE_stats) = c("No cutoff", "log2FC cutoff 1", "log2FC cutoff 2")
DE_stats = tibble::rownames_to_column(DE_stats, "log2FC_cutoff")
names(DE_stats) = c("log2FC_cutoff","FDR < 0.05")

###################################################################################
#                          Write Table
###################################################################################

write.table(DE_stats,        file = paste0(outprefix,".Summary_stats.txt"), sep = "\t", quote = F, row.names = F)
write.table(edgeR_DE,        file = paste0(outprefix,".DE_edgeR_All.txt"), sep = "\t", quote = F, row.names = F)
write.table(edgeR_DE_FDR5P,  file = paste0(outprefix,".DE_edgeR_5PCT.txt"), sep = "\t", quote = F, row.names = F)


###################################################################################
#                          Generate Input file for GO Term/Pathway Analyses
###################################################################################
# generate the input files at FDR 5P and at 2 Fold and 4 Fold thresholds to perform pathway analyses
my_list = list(
               "edgeR_DE_FDR5P" = edgeR_DE_FDR5P,
               "edgeR_DE_FDR5P_logFC1" = edgeR_DE_FDR5P_logFC1,
               "edgeR_DE_FDR5P_logFC2" = edgeR_DE_FDR5P_logFC2)


list_names = names(my_list)

for (name in names(my_list)) {
  
  l = my_list[[name]]
  v_CP = dplyr::select(l, "Gene_ID")
  outfile = paste0(outprefix,".CP_", name, ".txt")
  write.table(v_CP, outfile, row.names = F, quote = F, sep = "\t")
  
}

###################################################################################
#                          Pathway Analyses
###################################################################################
#
library(edgeR)
library(clusterProfiler)
library(org.Mm.eg.db)  # Use 'org.Mm.eg.db' for mouse data
library(pathview)

# Example edgeR results table
# head(edgeR_DE.show)

# Define thresholds
fdr_threshold <- 0.05
log2fc_threshold <- 1.0

# Filter significant genes
sig_genes <- edgeR_DE.show[edgeR_DE.show$FDR < fdr_threshold & abs(edgeR_DE.show$logFC) > log2fc_threshold, ]

# Convert gene symbols to Entrez IDs
gene_symbols <- sig_genes$Gene.name  # Assuming 'Gene.name' is the column with gene symbols
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)  
# Perform KEGG pathway enrichment analysis
kegg_enrich <- enrichKEGG(gene = entrez_ids$ENTREZID, 
                          organism = 'mmu',  # Use 'mmu' for mouse
                          pvalueCutoff = 0.05)

# View the enrichment results
head(kegg_enrich,n=3)

# Plot the first pathway in the result
if (length(kegg_enrich$ID) > 0) {
  pathview(gene.data = sig_genes$logFC, 
           pathway.id = kegg_enrich$ID[1], 
           species = "mmu")  # Use 'mmu' for mouse
}

kegg_df <- as.data.frame(kegg_enrich)

# Plot using ggplot2
ggplot(kegg_df, aes(x = Count, y = reorder(Description, Count), size = -log10(p.adjust), color = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(size = expression(-log[10](adjusted~p-value)),
       color = expression(-log[10](adjusted~p-value)),
       x = "Count",
       y = "ID",
       title = "KEGG Pathway Enrichment")
#boxplot
top_kegg <- kegg_df[1:20, ]
# Plot using ggplot2
ggplot(top_kegg, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = Count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(x = "Pathway", 
       y = expression(-log[10](adjusted~p-value)), 
       fill = "Gene Count",
       title = "Top 20 KEGG Pathway Enrichment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


######################## Other R functions to capture the session for reproducibility#########
## Document the various libraries loaded and capture the stderr
header = paste(rep("#", 50), collapse = "")

sink(file = paste0(outprefix,".sessioninfo.txt"))

cat(paste(header, "#Version Information", header, sep = '\n'))
cat('\n')
version
cat('\n')

cat(paste(header, "#Session Information", header, sep = '\n'))
cat('\n')
sessionInfo()
sink()

#################### Move all the analyzed files into a directory for cleaner hoarding ########
files<- list.files(pattern=outprefix)
dir.create(outprefix)
invisible(file.copy(files, outprefix))
invisible(file.remove(files))
