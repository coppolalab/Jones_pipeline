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

EigenCor <- function(module.traits.pval, status.col, status.vector) {
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

EigengeneBayes <- function(ME.vector, trait.vector) {
    trait.df <- data.frame(Trait = trait.vector, ME = ME.vector)
    trait.anova <- anovaBF(ME ~ Trait, data = trait.df) #%>% anova %>% tidy
    trait.anova
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

genes.rmcov <- readRDS.gz("../baseline_methylation/save/genes.compare.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

#sft for combat only = 10
sft <- pickSoftThreshold(t(genes.rmcov), powerVector = powers, verbose = 5, networkType = "signed")
sft.bicor <- pickSoftThreshold(t(genes.rmcov), powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
sft.spearman <- pickSoftThreshold(t(genes.rmcov), powerVector = powers, verbose = 5, corFnc = cor, corOptions = list(method = "spearman"), networkType = "signed")
sft.df <- sft.bicor$fitIndices
saveRDS.gz(sft.bicor, file = "./save/sft.bicor.rda")

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

softPower <- 20
adjacency.expr <- adjacency(t(genes.rmcov), power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
#saveRDS.gz(adjacency.expr, file = "./save/adjacency.rda")

TOM <- TOMsimilarity(adjacency.expr, verbose = 5)
dissimilarity.TOM <- 1 - TOM
rm(TOM)
gc()
#saveRDS.gz(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
gc()
saveRDS.gz(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 20

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
saveRDS.gz(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(t(genes.rmcov), colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
saveRDS.gz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50_tree", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(t(genes.rmcov), dynamic.colors, cutHeight = ME.dissimilarity.threshold, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 3) 
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("module_eigengene_clustering_min50", height = 10, width = 15, bg = "transparent")
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("Original Modules", "Merged Modules"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
saveRDS.gz(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
saveRDS.gz(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
saveRDS.gz(ME.genes, file = "./save/me.genes.rda")

CairoPDF("eigengenes", height = 6, width = 8, bg = "transparent")
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,5,1,2), plotPreservation = "standard")
dev.off()

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

genes.bm.table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description', 'gene_biotype'), filters = 'ensembl_gene_id', values = as.character(gene.info$Ensembl.ID), mart = ensembl)
genes.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(genes.bm.table) <- c("Ensembl.ID", "Symbol", "Description", "Gene.Type")

genes.bm.table.vega <- getBM(attributes = c('ens_gene', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ens_gene', values = as.character(gene.info$Ensembl.ID), mart = vega)
colnames(genes.bm.table.vega) <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

gene.info.annot <- left_join(gene.info, genes.bm.table) %>% left_join(genes.bm.table.vega) 
gene.module.membership <- data.frame(bicor(t(genes.rmcov), ME.genes, maxPOutliers = 0.05)) %>% signif(3)
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(t(genes.rmcov)))) %>% signif(3)
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")
gene.module.membership$Ensembl.ID <- rownames(gene.module.membership)
module.membership.pvalue$Ensembl.ID <- rownames(module.membership.pvalue)

final.genetable <- left_join(gene.info.annot, gene.module.membership) %>% 
    left_join(module.membership.pvalue, by = "Ensembl.ID") %>% 
    dplyr::select(Ensembl.ID, Symbol:Gene.Type, Module, kTotal:kscaled, matches("MM.*"), Vega.ID:Gene.Type.Vega) %>%
    arrange(Module, desc(kscaled))
Workbook(final.genetable, "./final_genetable_genes.xlsx")

#Plots
ME.genes.plot <- mutate(ME.genes, Combined = pheno.import$Combined) %>%
    gather(Module.Color, Eigengene, -Combined) 
ME.genes.plot$Module.Color %<>% str_replace("ME", "")
ME.genes.plot$Combined %<>% str_replace("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))

p <- ggplot(ME.genes.plot, aes(x = Combined, y = Eigengene, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
p <- p + scale_fill_manual(values = sort(unique(ME.genes.plot$Module.Color))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + facet_wrap(~ Module.Color, ncol = 5, scales = "free")
CairoPDF("eigengene_plots", height = 12, width = 16, bg = "transparent")
print(p)
dev.off()

black.plot <- filter(ME.genes.plot, Module.Color == "black")
p <- ggplot(black.plot, aes(x = Combined, y = Eigengene, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
CairoPDF("eigengene_black", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

black.mu <- data.frame(posterior.status$MEblack) %>% select(mu) %>% unlist
black.estimate <- data.frame(posterior.status$MEblack) %>% select(matches("^Trait")) %>% sweep(1, black.mu)
colnames(black.estimate) %<>% str_replace("Trait\\.", "")
black.estimate.plot <- gather(black.estimate, Combined, Estimate) 
black.estimate.plot$Combined %<>% str_replace("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))

p <- ggplot(black.estimate.plot, aes(x = Combined, y = Estimate, fill = Combined)) + geom_violin(color = "black", draw_quantiles = c(0.025, 0.5, 0.975)) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank()) +  ylab("Eigengene Estimate")
p <- p + theme(panel.border = element_rect(color = "black", size = 1)) + ggtitle("Black Module")
p <- p + theme(plot.title = element_text(hjust = 0.5))
CairoPDF("eigengene_black_estimate", height = 6, width = 7, bg = "transparent")
print(p)
dev.off()

salmon.plot <- filter(ME.genes.plot, Module.Color == "salmon")
p <- ggplot(salmon.plot, aes(x = Combined, y = Eigengene, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
CairoPDF("eigengene_salmon", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

salmon.mu <- data.frame(posterior.status$MEsalmon) %>% select(mu) %>% unlist
salmon.estimate <- data.frame(posterior.status$MEsalmon) %>% select(matches("^Trait")) %>% sweep(1, salmon.mu)
colnames(salmon.estimate) %<>% str_replace("Trait\\.", "")
salmon.estimate.plot <- gather(salmon.estimate, Combined, Estimate) 
salmon.estimate.plot$Combined %<>% str_replace("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))

p <- ggplot(salmon.estimate.plot, aes(x = Combined, y = Estimate, fill = Combined)) + geom_violin(color = "black", draw_quantiles = c(0.025, 0.5, 0.975)) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank()) + ylab("Eigengene Estimate")
p <- p + theme(panel.border = element_rect(color = "black", size = 1)) + ggtitle("salmon Module")
p <- p + theme(plot.title = element_text(hjust = 0.5)) + ylim(-0.3, 0.4)
CairoPDF("eigengene_salmon_estimate", height = 6, width = 7, bg = "transparent")
print(p)
dev.off()

dnam.filter <- ReadRDSgz("../baseline_methylation/save/dnam.filter.rda")
dnam.filter$black <- ME.genes$MEblack
dnam.filter$salmon <- ME.genes$MEsalmon

black.difference.cor <- bicor(dnam.filter$black, dnam.filter$Difference)
black.difference.pvalue <- corPvalueStudent(black.difference.cor, nrow(dnam.filter))
p <- ggplot(dnam.filter, aes(x = Difference, y = black, color = Combined)) + stat_smooth(color = "blue", method = "lm") + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(plot.title = element_text(hjust = 0.5), legend.background = element_blank()) + ylab("Eigenpromoter") + xlab("DNAm Age - Chronological Age")
p <- p + ggtitle(str_c("black Module (", expression(rho), " = ", signif(black.difference.cor,3), ", p < ", signif(black.difference.pvalue, 3), ")"))
CairoPDF("black_difference", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

salmon.difference.cor <- bicor(dnam.filter$salmon, dnam.filter$Difference)
salmon.difference.pvalue <- corPvalueStudent(salmon.difference.cor, nrow(dnam.filter))
p <- ggplot(dnam.filter, aes(x = Difference, y = salmon, color = Combined)) + stat_smooth(color = "blue", method = "lm") + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(plot.title = element_text(hjust = 0.5), legend.background = element_blank()) + ylab("Eigenpromoter") + xlab("DNAm Age - Chronological Age")
p <- p + ggtitle(str_c("salmon Module (", expression(rho), " = ", signif(salmon.difference.cor,3), ", p < ", signif(salmon.difference.pvalue, 3), ")"))
CairoPDF("salmon_difference", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

#Other Plots

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
cluster <- makeForkCluster(8)
split(expr.data.plot, expr.data.plot$module.colors) %>% map(Heatmap, ME.genes.plot)

modules.filter <- filter(final.genetable, nchar(Symbol) > 0) 
modules.submit <- select(modules.filter, Ensembl.ID, Symbol, Module)

PlotHubGenes("black", modules.filter, adjacency.expr)
PlotHubGenes("salmon", modules.filter, adjacency.expr)

source("../../code/GO/enrichr.R")
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016", "Human_Gene_Atlas", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up") 
color.names <- unique(module.colors) %>% sort
trap1 <- map(color.names, EnrichrSubmit, modules.submit, enrichr.terms, FALSE)

black.only <- filter(modules.submit, Module == "black")$Symbol
black.gobiol.file <- "./enrichr/black/black_GO_Biological_Process_2015.xlsx"
black.gobiol <- read.xlsx(black.gobiol.file) 
black.gobiol.filter <- FilterEnrichr(black.gobiol)
GetKappaCluster(file_path_sans_ext(black.gobiol.file), black.gobiol.filter, black.only)
black.gobiol.final <- slice(black.gobiol.filter, c(1, 46, 3, 6, 14))
black.gobiol.final$Database <- "GO Biological Process"

black.gomole.file <- "./enrichr/black/black_GO_Molecular_Function_2015.xlsx"
black.gomole <- read.xlsx(black.gomole.file) 
black.gomole.filter <- FilterEnrichr(black.gomole)
GetKappaCluster(file_path_sans_ext(black.gomole.file), black.gomole.filter, black.only)
black.gomole.final <- slice(black.gomole.filter, c(2))
black.gomole.final$Database <- "GO Molecular Function"

black.reactome.file <- "./enrichr/black/black_Reactome_2016.xlsx"
black.reactome <- read.xlsx(black.reactome.file) 
black.reactome.filter <- FilterEnrichr(black.reactome)
GetKappaCluster(file_path_sans_ext(black.reactome.file), black.reactome.filter, black.only)
black.reactome.final <- slice(black.reactome.filter, c(9, 13))
black.reactome.final$Database <- "Reactome"

black.enrichr.final <- rbind(black.gobiol.final, black.gomole.final, black.reactome.final)
EnrichrPlot(black.enrichr.final, "black.enrichr", "black Module")

salmon.only <- filter(modules.submit, Module == "salmon")$Symbol
salmon.gobiol.file <- "./enrichr/salmon/salmon_GO_Biological_Process_2015.xlsx"
salmon.gobiol <- read.xlsx(salmon.gobiol.file) 
salmon.gobiol$Num.Genes <- map(salmon.gobiol$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
salmon.gobiol.filter <- FilterEnrichr(salmon.gobiol)
GetKappaCluster(file_path_sans_ext(salmon.gobiol.file), salmon.gobiol.filter, salmon.only)
salmon.gobiol.final <- slice(salmon.gobiol, c(15, 36, 2, 78))
salmon.gobiol.final$Database <- "GO Biological Process"

salmon.gomole.file <- "./enrichr/salmon/salmon_GO_Molecular_Function_2015.xlsx"
salmon.gomole <- read.xlsx(salmon.gomole.file) 
salmon.gomole$Num.Genes <- map(salmon.gomole$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
salmon.gomole.filter <- FilterEnrichr(salmon.gomole)
GetKappaCluster(file_path_sans_ext(salmon.gomole.file), salmon.gomole.filter, salmon.only)
salmon.gomole.final <- slice(salmon.gomole, c(8))
salmon.gomole.final$Database <- "GO Molecular Function"

salmon.reactome.file <- "./enrichr/salmon/salmon_Reactome_2016.xlsx"
salmon.reactome <- read.xlsx(salmon.reactome.file) 
salmon.reactome$Num.Genes <- map(salmon.reactome$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
salmon.reactome.filter <- FilterEnrichr(salmon.reactome)
GetKappaCluster(file_path_sans_ext(salmon.reactome.file), salmon.reactome.filter, salmon.only)
salmon.reactome.final <- slice(salmon.reactome, c(17))
salmon.reactome.final$Database <- "Reactome"

salmon.enrichr.final <- rbind(salmon.gobiol.final, salmon.gomole.final, salmon.reactome.final)
EnrichrPlot(salmon.enrichr.final, "salmon.enrichr", "Salmon Module")

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

