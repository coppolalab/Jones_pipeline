#For WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()
library(biomaRt)

#For baseline processing
library(limma)
library(sva)
library(R.utils)
library(Biobase)

#Functional programming
library(magrittr)
library(purrr)
library(functional)

#Data arrangement
library(dplyr)
library(tidyr)
library(broom)
library(PMCMR)
library(BayesFactor)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(igraph)

#Reading and writing tables
library(readr)
library(openxlsx)
library(parallel)

source("../../FRDA project/common_functions.R")

PCAPlot <- function(filename, dataset, facet.bool, size.height, size.width) {
    colnames(dataset)[2] <- "Module"
    dataset$Module %<>% str_replace("ME", "") 
    p <- ggplot(dataset, aes(x = as.numeric(x), y = value, fill = Module, color = Module)) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("First Principal Component")
    p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(dataset$x)))
    p <- p + scale_color_manual(values = sort(unique(dataset$Module)))
    if (facet.bool == TRUE)
    {
        p <- p + facet_wrap(~ Module)
        p <- p + theme(legend.position = "none")
    } 
    CairoPDF(filename, height = size.height, width = size.width)
    print(p)
    dev.off()
}

Heatmap <- function(dataset, ME.genes) {
    color <- as.character(unique(dataset$module.colors))
    dataset %<>% select(-module.colors) %>% scale
    max.dataset <- max(abs(dataset))
    print(dim(dataset))
    CairoPDF(paste("./modules/", color, sep = ""), width = 21, height = 12)
    par(mar = c(3.5,3,2,3))
    par(oma = c(4,0,2,0))
    plotMat(dataset, zlim = c(-max.dataset, max.dataset), main = paste(color, " (", nrow(dataset), ")", sep = ""))

    ME.genes.plot <- select(ME.genes, Sample.ID, matches(color))
    p <- ggplot(ME.genes.plot, aes_string(x = "Sample.ID", y = color))
    p <- p + geom_bar(stat = "identity") + xlab("Eigengene Expression")#+ ylim(c(-6, 16)) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.x = element_text(angle = 90, size = 2))  
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    print(p)
    dev.off()
}

EnrichrSubmit <- function(index, full.df, enrichr.terms, use.weights = FALSE) {
    dataset <- filter(full.df, Module == index)
    dir.create(file.path("./enrichr", index), recursive = TRUE, showWarnings = FALSE)
    enrichr.data <- map(enrichr.terms, get.enrichrdata, dataset, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), EnrichrWorkbook, enrichr.data, index)
}

EnrichrWorkbook <- function(subindex, full.df, index) {
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    filename = paste("./enrichr/", index, "/", index, "_", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

EigenPlot <- function(module.traits.pval, status.col, status.vector) {
    sig.col <- paste(status.col, ".p.value", sep = "")
    cor.status.labeled <- data.frame(Color = rownames(module.traits.pval), select_(data.frame(module.traits.pval), sig.col))
    filter.cond <- paste(sig.col, "< 0.05")
    status.sig <- filter_(cor.status.labeled, filter.cond)
    me.genes.status <- select(ME.genes, one_of(as.character(status.sig$Color)))
    me.genes.status$Status <- status.vector
    split.cols <- str_split(status.col, "\\.")[[1]]
    me.status.melt <- melt(me.genes.status, id.vars = "Status") %>% filter(Status == split.cols[1] | Status == split.cols[2])
    colnames(me.status.melt)[2] <- "Module"

    sum.fun <- function(data.vector){ data.frame(ymin = min(data.vector), ymax = max(data.vector), y = mean(data.vector)) }
    me.status.melt$Module %<>% as.character
    p <- ggplot(me.status.melt, aes(x = factor(Status), y = value, col = Module)) + geom_point(position = "jitter")
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("Eigengene")
    p <- p + theme(axis.title.x = element_blank()) + stat_summary(aes(group = 1), fun.y = mean, geom = "line", col = "black", position = position_dodge(width = 0.9))
    p <- p + scale_color_manual(values = sort(unique(me.status.melt$Module)))
    p <- p + facet_wrap(~ Module, scales = "free_y")
    p <- p + theme(legend.position = "none")

    filename <- paste(split.cols[1], "_", split.cols[2], "_eigengenes_05", sep = "")
    CairoPDF(filename, height = 13, width = 20)
    print(p)
    dev.off()
}

Workbook <- function(dataset, filename) {
    pval.detect <- colnames(dataset) %>% str_detect("pvalue") 
    pval.cols <- which(pval.detect)
    cor.detect <- colnames(dataset) %>% str_detect("MM.*") 
    cor.cols <- which(cor.detect & !(pval.detect))
    description.cols <- colnames(dataset) %>% str_detect("Description") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.005", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = cor.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = 20)
    setColWidths(wb, 1, cols = 2, widths = "auto")
    setColWidths(wb, 1, cols = 4:ncol(dataset), widths = "auto")
    setColWidths(wb, 1, cols = description.cols, widths = 45)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

EigengeneModel <- function(ME.vector, trait.vector, contrasts.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.contrasts <- posthoc.kruskal.dunn.test(ME ~ Trait, trait.df, p.adjust.method = "none") 
    trait.contrasts

    trait.pvals <- trait.contrasts$p.value
    pval.ccdc <- trait.pvals["colon.crypts","duodenum.crypts"]
    pval.cmdm <- trait.pvals["colon.mucosa","duodenum.mucosa"]
    pval.dcdm <- trait.pvals["duodenum.crypts","duodenum.mucosa"]
    pval.cccm <- trait.pvals["colon.crypts","colon.mucosa"]
    pvals.subset <- p.adjust(c(pval.ccdc, pval.cmdm, pval.dcdm, pval.cccm), method = "fdr")

    trait.medians <- group_by(trait.df, Trait) %>% summarise(median(ME)) %>% data.frame
    colnames(trait.medians)[2] <- "Median"
    rownames(trait.medians) <- trait.medians$Trait
    pval.ccdc <- trait.medians["duodenum.crypts","Median"] - trait.medians["colon.crypts","Median"] 
    pval.cmdm <- trait.medians["duodenum.mucosa","Median"] - trait.medians["colon.mucosa","Median"] 
    pval.dcdm <- trait.medians["duodenum.mucosa","Median"] - trait.medians["duodenum.crypts","Median"] 
    pval.cccm <- trait.medians["colon.mucosa","Median"] - trait.medians["colon.crypts","Median"] 

    diffs.subset <- c(pval.ccdc, pval.cmdm, pval.dcdm, pval.cccm)
    anova.return <- data.frame(Diff = diffs.subset, P.value = pvals.subset)
    rownames(anova.return) <- contrasts.vector

    return(anova.return)
}

EigengeneANOVA <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- kruskal.test(ME ~ Trait, trait.df) %>% tidy
    return(trait.anova$p.value)
}

PlotHubGenes <- function(module.color, module.table, adjacency.mat) {
    filter.string <- str_c("Gene.Type == 'protein_coding' & Module == '", module.color, "'")
    module.select <- filter_(module.table, filter.string) %>% slice(1:15)
    adjacency.select <- adjacency.expr[module.select$Ensembl.ID,module.select$Ensembl.ID] 
    rownames(adjacency.select) <- module.select$Symbol
    colnames(adjacency.select) <- module.select$Symbol
    adjacency.igraph <- graph_from_adjacency_matrix(adjacency.select, mode = "undirected", weighted = TRUE, diag = FALSE)
    layout.igraph <- graph_from_adjacency_matrix(matrix(rep(1, nrow(adjacency.select)^2), nrow = nrow(adjacency.select), ncol = ncol(adjacency.select)), mode = "undirected", diag = FALSE)
    CairoPDF(str_c(module.color, "_igraph"), width = 6, height = 6, bg = "transparent")
    plot.igraph(adjacency.igraph, layout = layout_with_fr(layout.igraph), vertex.size = 3, vertex.label.degree = -pi/2, vertex.label.dist = 0.5, vertex.label.font = 2, vertex.label.color = "black", edge.color = "green", edge.width = E(adjacency.igraph)$weight ^ (1/14))
    dev.off()
}

FilterEnrichr <- function(enrichr.df, size = 100) {
    enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
    enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
    if (nrow(enrichr.df) > size) {
        enrichr.filter %<>% slice(1:size)
    }
    enrichr.filter
}

promoters.rmcov <- ReadRDSgz("../baseline_methylation/save/promoters.compare.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

sft <- pickSoftThreshold(t(promoters.rmcov), powerVector = powers, verbose = 5, networkType = "signed")
sft.bicor <- pickSoftThreshold(t(promoters.rmcov), powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
sft.spearman <- pickSoftThreshold(t(promoters.rmcov), powerVector = powers, verbose = 5, corFnc = cor, corOptions = list(method = "spearman"), networkType = "signed")
sft.df <- sft.bicor$fitIndices
SaveRDSgz(sft, file = "./save/sft.rda")

#Use bicor, compare 14 and 16 for soft power
#Plot scale indendence and mean connectivity as functions of power
sft.df$multiplied <- sft.df$SFT.R.sq * -sign(sft.df$slope)
p <- ggplot(sft.df, aes(x = Power,  y = multiplied, label = Power)) + geom_point() + geom_text(vjust = -0.6, size = 4, col = "red")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(aes(yintercept = 0.9))
p <- p + xlab("Soft Threshold") + ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale Independence")
CairoPDF(file = "./scaleindependence", width = 6, height = 6)
print(p)
dev.off()

p <- ggplot(sft.df, aes(x = Power,  y = mean.k.)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Soft Threshold") + ylab("Mean Connectivity") + ggtitle("Mean Connectivity")
CairoPDF(file = "./meanconnectivity", width = 6, height = 6)
print(p)
dev.off()

softPower <- 14
adjacency.expr <- adjacency(t(promoters.rmcov), power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
#SaveRDSgz(adjacency.expr, file = "./save/adjacency.rda")

TOM <- TOMsimilarity(adjacency.expr, verbose = 5)
dissimilarity.TOM <- 1 - TOM
rm(TOM)
gc()
#SaveRDSgz(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
gc()
SaveRDSgz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 20

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
SaveRDSgz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(t(promoters.rmcov), colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
SaveRDSgz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50_tree", height = 10, width = 15, bg = "transparent")
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(t(promoters.rmcov), dynamic.colors, cutHeight = ME.dissimilarity.threshold, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 3) 
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering_min50", height = 10, width = 15, bg = "transparent")
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("Original Modules", "Merged Modules"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
SaveRDSgz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
SaveRDSgz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
SaveRDSgz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 6, width = 8, bg = "transparent")
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,5,1,2), plotPreservation = "standard")
dev.off()

EigengeneBayes <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- anovaBF(ME ~ Trait, data = trait.df) #%>% anova %>% tidy
    trait.anova
}

pheno.import <- readRDS.gz("../baseline_methylation/save/pheno.compare.rda")
pheno.import$Combined %<>% factor %>% droplevels
anova.status <- map_dbl(ME.genes, EigengeneANOVA, pheno.import$Combined) %>% p.adjust("fdr") %>% signif(3)
bayes.status <- map(ME.genes, EigengeneBayes, pheno.import$Combined) 
bf.status <- map(bayes.status, extractBF) %>% map_dbl(extract2, "bf")
posterior.status <- map(bayes.status, posterior, iterations = 100000) 
contrasts.vector <- c("duodenum_crypts_vs._colon_crypts", "duodenum_mucosa_vs._colon_mucosa", "duodenum_mucosa_vs._duodenum_crypts", "colon_mucosa_vs._colon_crypts")

color.values <- unique(module.colors)
cor.status <- map(ME.genes, EigengeneModel, pheno.import$Combined, contrasts.vector)

status.diff <- map(cor.status, select, Diff) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.diff) <- names(cor.status)

status.pval <- map(cor.status, select, P.value) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.pval) <- names(cor.status)

pval.adjust <- map(status.pval, p.adjust, method = "fdr", n = length(color.values)) %>% reduce(cbind) %>% data.frame
rownames(pval.adjust) <- rownames(status.pval)
colnames(pval.adjust) <- paste(colnames(status.pval), ".pval")

text.matrix.traits <- paste(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(pval.adjust), 1), ')', sep = '')
dim(text.matrix.traits) = dim(status.diff)

heatmap.range <- c(min(as.matrix(status.diff)) * 1.1, max(as.matrix(status.diff)) * 1.1)
width.dynamic <- 3 + (1 * ncol(text.matrix.traits))
CairoPDF("module-trait relationships", width = width.dynamic, height = 10, bg = "transparent")
par(mar = c(14, 9, 3, 3))
labeledHeatmap(Matrix = as.matrix(status.diff), xLabels = colnames(status.diff), yLabels = colnames(ME.genes), ySymbols = colnames(ME.genes), yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix.traits, setStdMargins = F, zlim = c(-1,1), main = "")
dev.off()

all.degrees <- intramodularConnectivity(adjacency.expr, module.colors)
gene.info <- data.frame(Ensembl.ID = rownames(all.degrees), Module = module.colors, all.degrees)
gene.info$kscaled <- by(gene.info, gene.info$Module, select, kWithin) %>% map(function(kWithin) kWithin / max(kWithin)) %>% reduce(c)
gene.info[,3:ncol(gene.info)] <- signif(gene.info[,3:ncol(gene.info)], 3)

#Annotate gene table
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
vega <- useMart("ENSEMBL_MART_VEGA", dataset = "hsapiens_gene_vega")

promoters.bm.table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description', 'gene_biotype'), filters = 'ensembl_gene_id', values = as.character(gene.info$Ensembl.ID), mart = ensembl)
promoters.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(promoters.bm.table) <- c("Ensembl.ID", "Symbol", "Description", "Gene.Type")

promoters.bm.table.vega <- getBM(attributes = c('ens_gene', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ens_gene', values = as.character(gene.info$Ensembl.ID), mart = vega)
colnames(promoters.bm.table.vega) <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

gene.info.annot <- left_join(gene.info, promoters.bm.table) %>% left_join(promoters.bm.table.vega) 
gene.module.membership <- data.frame(bicor(t(promoters.rmcov), ME.genes, maxPOutliers = 0.05)) %>% signif(3)
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(t(promoters.rmcov)))) %>% signif(3)
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")
gene.module.membership$Ensembl.ID <- rownames(gene.module.membership)
module.membership.pvalue$Ensembl.ID <- rownames(module.membership.pvalue)

final.genetable <- left_join(gene.info.annot, gene.module.membership) %>% 
    left_join(module.membership.pvalue, by = "Ensembl.ID") %>% 
    dplyr::select(Ensembl.ID, Symbol:Gene.Type, Module, kTotal:kscaled, matches("MM.*"), Vega.ID:Gene.Type.Vega) %>%
    arrange(Module, desc(kscaled))
Workbook(final.genetable, "./final_genetable_promoters.xlsx")

#Plots
ME.genes.plot <- mutate(ME.genes, Combined = pheno.import$Combined) %>%
    gather(Module.Color, Eigengene, -Combined) 
ME.genes.plot$Module.Color %<>% str_replace("ME", "")
ME.genes.plot$Combined %<>% str_replace("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))

p <- ggplot(ME.genes.plot, aes(x = Combined, y = Eigengene, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + scale_fill_manual(values = sort(unique(ME.genes.plot$Module.Color)))
p <- p + theme(panel.border = element_rect(size = 1, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + facet_wrap(~ Module.Color, ncol = 5, scales = "free")
CairoPDF("eigengene_plots", height = 12, width = 16, bg = "transparent")
print(p)
dev.off()

red.plot <- filter(ME.genes.plot, Module.Color == "red")
p <- ggplot(red.plot, aes(x = Combined, y = Eigengene, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
CairoPDF("eigengene_red", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

red.mu <- data.frame(posterior.status$MEred) %>% select(mu) %>% unlist
red.estimate <- data.frame(posterior.status$MEred) %>% select(matches("^Trait")) %>% sweep(1, red.mu)
colnames(red.estimate) %<>% str_replace("Trait\\.", "")
red.estimate.plot <- gather(red.estimate, Combined, Estimate) 
red.estimate.plot$Combined %<>% str_replace("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))
red.estimate.plot %<>% filter(!(Combined == "duodenum crypts" & Estimate > 0)) %>% filter(!(Combined == "duodenum mucosa" & Estimate < -0.15))

p <- ggplot(red.estimate.plot, aes(x = Combined, y = Estimate, fill = Combined)) + geom_violin(color = "black", draw_quantiles = c(0.025, 0.5, 0.975)) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank()) + ylim(c(-0.3, 0.3)) + ylab("Eigengene Estimate")
p <- p + theme(panel.border = element_rect(color = "black", size = 1)) + ggtitle("Red Module")
p <- p + theme(plot.title = element_text(hjust = 0.5))
CairoPDF("eigengene_red_estimate", height = 6, width = 7, bg = "transparent")
print(p)
dev.off()

greenyellow.mu <- data.frame(posterior.status$MEgreenyellow) %>% select(mu) %>% unlist
greenyellow.estimate <- data.frame(posterior.status$MEgreenyellow) %>% select(matches("^Trait")) %>% sweep(1, greenyellow.mu)
colnames(greenyellow.estimate) %<>% str_replace("Trait\\.", "")
greenyellow.estimate.plot <- gather(greenyellow.estimate, Combined, Estimate) 
greenyellow.estimate.plot$Combined %<>% str_replace("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))

p <- ggplot(greenyellow.estimate.plot, aes(x = Combined, y = Estimate, fill = Combined)) + geom_violin(color = "black", draw_quantiles = c(0.025, 0.5, 0.975)) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank()) + ylim(c(-0.3, 0.37)) + ylab("Eigengene Estimate")
p <- p + theme(panel.border = element_rect(color = "black", size = 1)) + ggtitle("GreenYellow Module")
p <- p + theme(plot.title = element_text(hjust = 0.5))
CairoPDF("eigengene_greenyellow_estimate", height = 6, width = 7, bg = "transparent")
print(p)
dev.off()

greenyellow.plot <- filter(ME.genes.plot, Module.Color == "greenyellow")
p <- ggplot(greenyellow.plot, aes(x = Combined, y = Eigengene, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
CairoPDF("eigengene_greenyellow", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

dnam.filter <- ReadRDSgz("../baseline_methylation/save/dnam.filter.rda")
dnam.filter$red <- ME.genes$MEred
dnam.filter$greenyellow <- ME.genes$MEgreenyellow

red.difference.cor <- bicor(dnam.filter$red, dnam.filter$Difference)
red.difference.pvalue <- corPvalueStudent(red.difference.cor, nrow(dnam.filter))
p <- ggplot(dnam.filter, aes(x = Difference, y = red, color = Combined)) + stat_smooth(color = "blue", method = "lm") + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(plot.title = element_text(hjust = 0.5), legend.background = element_blank()) + ylab("Eigenpromoter") + xlab("DNAm Age - Chronological Age")
p <- p + ggtitle(str_c("Red Module (", expression(rho), " = ", signif(red.difference.cor,3), ", p < ", signif(red.difference.pvalue, 3), ")"))
CairoPDF("red_difference", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

greenyellow.difference.cor <- bicor(dnam.filter$greenyellow, dnam.filter$Difference)
greenyellow.difference.pvalue <- corPvalueStudent(greenyellow.difference.cor, nrow(dnam.filter))
p <- ggplot(dnam.filter, aes(x = Difference, y = greenyellow, color = Combined)) + stat_smooth(color = "blue", method = "lm") + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(plot.title = element_text(hjust = 0.5), legend.background = element_blank()) + ylab("Eigenpromoter") + xlab("DNAm Age - Chronological Age")
p <- p + ggtitle(str_c("GreenYellow Module (", expression(rho), " = ", signif(greenyellow.difference.cor,3), ", p < ", signif(greenyellow.difference.pvalue, 3), ")"))
CairoPDF("greenyellow_difference", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

#Other plots
all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% llply(`[`, "y")
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$x <- as.factor(1:nrow(smooth.df))
smooth.plot <- melt(smooth.df, id.vars = "x")

PCAPlot("all_principal_components", smooth.plot, FALSE, 10, 15)
PCAPlot("facet_principal_components", smooth.plot, TRUE, 13, 25)

sample.ids <- factor(rownames(expr.collapse), levels = rownames(expr.collapse))
colnames(ME.genes) %<>% str_replace("ME", "")
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(expr.collapse), module.colors)
split(expr.data.plot, expr.data.plot$module.colors) %>% map(Heatmap, ME.genes.plot)

modules.filter <- filter(final.genetable, nchar(Symbol) > 0) 
PlotHubGenes("red", modules.filter, adjacency.expr)
PlotHubGenes("greenyellow", modules.filter, adjacency.expr)

EnrichrPlot <- function(enrichr.df, filename, plot.title, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    enrichr.df$Adj.P.value <- p.adjust(enrichr.df$P.value, method = "fdr")
    enrichr.df$Log.P.value <- -log10(enrichr.df$Adj.P.value)
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Database, ": ", enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(Log.P.value)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)

    p <- ggplot(enrichr.df, aes(Format.Name, Log.P.value)) + geom_bar(fill = "blue", stat = "identity", size = 1) 
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(plot.background = element_blank(), legend.background = element_blank(), axis.text.y = element_text(size = 12), panel.border = element_blank())
    p <- p + theme(panel.background = element_blank(), axis.line.x = element_line()) + geom_hline(color = "red", yintercept = -log10(0.05))
    p <- p + ylab(expression(paste(Log[10], ' P-Value')))
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

#Enrichr
source("../../code/GO/enrichr.R")
modules.submit <- select(modules.filter, Ensembl.ID, Symbol, Module)
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016", "Human_Gene_Atlas", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up") 
color.names <- unique(module.colors) %>% sort
trap1 <- map(color.names, EnrichrSubmit, modules.submit, enrichr.terms, FALSE)

red.only <- filter(modules.submit, Module == "red")$Symbol
red.gobiol.file <- "./enrichr/red/red_GO_Biological_Process_2015.xlsx"
red.gobiol <- read.xlsx(red.gobiol.file) 
red.gobiol.filter <- FilterEnrichr(red.gobiol)
GetKappaCluster(file_path_sans_ext(red.gobiol.file), red.gobiol.filter, red.only)
red.gobiol.final <- slice(red.gobiol.filter, c(14, 1, 2, 11))
red.gobiol.final$Database <- "GO Biological Process"

red.gomole.file <- "./enrichr/red/red_GO_Molecular_Function_2015.xlsx"
red.gomole <- read.xlsx(red.gomole.file) 
red.gomole.filter <- FilterEnrichr(red.gomole)
GetKappaCluster(file_path_sans_ext(red.gomole.file), red.gomole.filter, red.only)
red.gomole.final <- slice(red.gomole.filter, c(5))
red.gomole.final$Database <- "GO Molecular Function"

red.reactome.file <- "./enrichr/red/red_Reactome_2016.xlsx"
red.reactome <- read.xlsx(red.reactome.file) 
red.reactome.filter <- FilterEnrichr(red.reactome)
GetKappaCluster(file_path_sans_ext(red.reactome.file), red.reactome.filter, red.only)
red.reactome.final <- slice(red.reactome.filter, c(2))
red.reactome.final$Database <- "Reactome"

red.enrichr.final <- rbind(red.gobiol.final, red.gomole.final, red.reactome.final)
EnrichrPlot(red.enrichr.final, "red.enrichr", "Red Module")

greenyellow.only <- filter(modules.submit, Module == "greenyellow")$Symbol
greenyellow.gobiol.file <- "./enrichr/greenyellow/greenyellow_GO_Biological_Process_2015.xlsx"
greenyellow.gobiol <- read.xlsx(greenyellow.gobiol.file) 
greenyellow.gobiol.filter <- FilterEnrichr(greenyellow.gobiol)
GetKappaCluster(file_path_sans_ext(greenyellow.gobiol.file), greenyellow.gobiol.filter, greenyellow.only)
greenyellow.gobiol.final <- slice(greenyellow.gobiol.filter, c(6, 4, 1, 18))
greenyellow.gobiol.final$Database <- "GO Biological Process"

#SKIP
greenyellow.only <- filter(modules.submit, Module == "greenyellow")$Symbol
greenyellow.gomole.file <- "./enrichr/greenyellow/greenyellow_GO_Molecular_Function_2015.xlsx"
greenyellow.gomole <- read.xlsx(greenyellow.gomole.file) 
greenyellow.gomole.filter <- FilterEnrichr(greenyellow.gomole)
GetKappaCluster(file_path_sans_ext(greenyellow.gomole.file), greenyellow.gomole.filter, greenyellow.only)

greenyellow.reactome.file <- "./enrichr/greenyellow/greenyellow_Reactome_2016.xlsx"
greenyellow.reactome <- read.xlsx(greenyellow.reactome.file) 
greenyellow.reactome.filter <- FilterEnrichr(greenyellow.reactome)
GetKappaCluster(file_path_sans_ext(greenyellow.reactome.file), greenyellow.reactome.filter, greenyellow.only)
greenyellow.reactome.final <- slice(greenyellow.reactome.filter, c(3))
greenyellow.reactome.final$Database <- "Reactome"

greenyellow.kegg.file <- "./enrichr/greenyellow/greenyellow_KEGG_2016.xlsx"
greenyellow.kegg <- read.xlsx(greenyellow.kegg.file) 
greenyellow.kegg.filter <- FilterEnrichr(greenyellow.kegg)
GetKappaCluster(file_path_sans_ext(greenyellow.kegg.file), greenyellow.kegg.filter, greenyellow.only)
greenyellow.kegg.final <- slice(greenyellow.kegg.filter, c(1))
greenyellow.kegg.final$Database <- "KEGG"

greenyellow.enrichr.final <- rbind(greenyellow.gobiol.final, greenyellow.kegg.final, greenyellow.reactome.final)
EnrichrPlot(greenyellow.enrichr.final, "greenyellow.enrichr", "greenyellow Module")

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

#Clock genes
clock.transcripts <- read.xlsx("../baseline_methylation/clock.transcripts.xlsx")
clock.coef.table <- read_csv("../baseline_methylation/13059_2013_3156_MOESM23_ESM.csv")
colnames(clock.coef.table)[1] <- "Illumina.ID"
clock.join <- left_join(clock.transcripts, clock.coef.table)
clock.red <- filter(clock.join, nearestGeneSymbol %in% red.only) %>% arrange(desc(abs(CoefficientTraining)))
clock.greenyellow <- filter(clock.join, nearestGeneSymbol %in% greenyellow.only) %>% arrange(desc(abs(CoefficientTraining)))


