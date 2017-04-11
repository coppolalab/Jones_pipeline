library(WGCNA)
library(flashClust)
enableWGCNAThreads()
library(biomaRt)
library(Biobase)

library(PMCMR)
library(BayesFactor)

library(Cairo)
library(igraph)
library(openxlsx)

library(stringr)
library(magrittr)
library(broom)
library(tidyverse)

source("../../code/common_functions.R")

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
    enrichr.data <- map(enrichr.terms, GetEnrichrData, dataset, FALSE)
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

EigengeneBayes <- function(ME.vector, trait.df) {
    trait.final <- mutate(trait.df, ME = ME.vector)
    trait.final$Batch %<>% factor
    full.lm <- lmBF(ME ~ Combined + Age.cont + Sex + Batch, data = trait.final) 
    noage.lm <- lmBF(ME ~ Combined + Sex + Batch, data = trait.final) 
    covar.lm <- lmBF(ME ~ Age.cont + Sex + Batch, data = trait.final) 
    c(full.lm, covar.lm, noage.lm)
}

Top5Plot <- function(rank.column, module.object, meth.object, pheno.object, pheno.col, levels.vector, maintitle, filename) {
    col.name <- str_c("desc(", rank.column, ")")

    top5.ensembl <- arrange_(module.object, col.name)$Ensembl.ID[1:10]
    top5.symbol <- arrange_(module.object, col.name)$Symbol[1:10]
    top5.expr <- t(meth.object[top5.ensembl,])
    colnames(top5.expr) <- top5.symbol
    top5.df <- data.frame(Sample = pheno.object[[pheno.col]], top5.expr) %>% gather(Gene, Expression, -Sample)

    top5.df$Gene %<>% factor(levels = make.names(top5.symbol))

    p <- ggplot(top5.df, aes(x = Sample, y = Expression, color = Sample)) + geom_boxplot() + geom_jitter() + theme_bw()
    p <- p + facet_wrap(~ Gene, ncol = 5, scales = "free") + theme(legend.position = "none")
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
    p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
    p <- p + theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank())
    p <- p + ggtitle(maintitle)
    CairoPDF(filename, height = 8, width = 16, bg = "transparent")
    print(p)
    dev.off()
}

FilterEnrichr <- function(enrichr.df, size = 200) {
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

#PlotHubGenes <- function(module.color, module.table, adjacency.mat) {
    #filter.string <- str_c("Gene.Type == 'protein_coding' & Module == '", module.color, "'")
    #module.select <- filter_(module.table, filter.string) %>% slice(1:15)
    #adjacency.select <- adjacency.expr[module.select$Ensembl.ID,module.select$Ensembl.ID] 
    #rownames(adjacency.select) <- module.select$Symbol
    #colnames(adjacency.select) <- module.select$Symbol
    #adjacency.igraph <- graph_from_adjacency_matrix(adjacency.select, mode = "undirected", weighted = TRUE, diag = FALSE)
    #layout.igraph <- graph_from_adjacency_matrix(matrix(rep(1, nrow(adjacency.select)^2), nrow = nrow(adjacency.select), ncol = ncol(adjacency.select)), mode = "undirected", diag = FALSE)
    #CairoPDF(str_c(module.color, "_igraph"), width = 6, height = 6, bg = "transparent")
    #plot.igraph(adjacency.igraph, layout = layout_with_fr(layout.igraph), vertex.size = 3, vertex.label.degree = -pi/2, vertex.label.dist = 0.5, vertex.label.font = 2, vertex.label.color = "black", edge.color = "green", edge.width = E(adjacency.igraph)$weight ^ (1/14))
    #dev.off()
#}

promoters.compare <- ReadRDSgz("../baseline_methylation/save/promoters.compare.combat.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

sft <- pickSoftThreshold(t(promoters.compare), powerVector = powers, verbose = 5, networkType = "signed")
sft.bicor <- pickSoftThreshold(t(promoters.compare), powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
sft.spearman <- pickSoftThreshold(t(promoters.compare), powerVector = powers, verbose = 5, corFnc = cor, corOptions = list(method = "spearman"), networkType = "signed")
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
adjacency.expr <- adjacency(t(promoters.compare), power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
#SaveRDSgz(adjacency.expr, file = "./save/adjacency.rda")
gc()

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
ME.list <- moduleEigengenes(t(promoters.compare), colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
SaveRDSgz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50_tree", height = 10, width = 15, bg = "transparent")
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(t(promoters.compare), dynamic.colors, cutHeight = ME.dissimilarity.threshold, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 3) 
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

pheno.import <- readRDS.gz("../baseline_methylation/save/pheno.compare.rda")
pheno.import$Combined <- str_c(pheno.import$Tissue, pheno.import$Cell.Type, sep = " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))
anova.status <- map_dbl(ME.genes, EigengeneANOVA, pheno.import$Combined) %>% p.adjust("fdr") %>% signif(3)

bayes.status <- map(ME.genes, EigengeneBayes, pheno.import) 
bf.status <- map(bayes.status, extractBF) 
combined.bf <- map(bf.status, extract2, "bf") %>% map(magrittr::extract, c(1,2)) %>% map_dbl(reduce, divide_by) %>% log10 %>% signif(3)
age.bf <- map(bf.status, extract2, "bf") %>% map(magrittr::extract, c(1,3)) %>% map_dbl(reduce, divide_by) %>% log10 %>% signif(3)
posterior.status <- map(bayes.status, magrittr::extract, 1) %>% map(posterior, iterations = 100000) 
contrasts.vector <- c("duodenum_crypts_vs._colon_crypts", "duodenum_mucosa_vs._colon_mucosa", "duodenum_mucosa_vs._duodenum_crypts", "colon_mucosa_vs._colon_crypts")

color.values <- unique(module.colors)
cor.status <- map(ME.genes, EigengeneModel, pheno.import$Combined, contrasts.vector)

cor.age <- bicor(ME.genes, pheno.import$Age) %>% as.vector
cor.age.pval <- corPvalueStudent(cor.age, nrow(pheno.import)) %>% p.adjust("fdr")
cor.df <- tibble(Module = colnames(ME.genes), Correlation = cor.age, P.Value = cor.age.pval)

status.diff <- map(cor.status, select, Diff) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.diff) <- names(cor.status)

status.pval <- map(cor.status, select, P.value) %>% map(t) %>% reduce(rbind) %>% data.frame
rownames(status.pval) <- names(cor.status)

pval.adjust <- map(status.pval, p.adjust, method = "fdr", n = length(color.values)) %>% reduce(cbind) %>% data.frame
rownames(pval.adjust) <- rownames(status.pval)
colnames(pval.adjust) <- paste(colnames(status.pval), ".pval")

text.matrix.traits <- paste(signif(as.matrix(status.diff), 2), '\n(', signif(as.matrix(status.pval), 1), ')', sep = '')
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

promoters.bm.table.vega <- getBM(attributes = c('ensg', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ensg', values = as.character(gene.info$Ensembl.ID), mart = vega)
colnames(promoters.bm.table.vega) <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

gene.info.annot <- left_join(gene.info, promoters.bm.table) %>% left_join(promoters.bm.table.vega) 
gene.module.membership <- data.frame(bicor(t(promoters.compare), ME.genes, maxPOutliers = 0.05)) %>% signif(3)
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(t(promoters.compare)))) %>% signif(3)
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
ME.genes.plot$Combined %<>% str_replace("\\_", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))

p <- ggplot(ME.genes.plot, aes(x = Combined, y = Eigengene, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + scale_fill_manual(values = sort(unique(ME.genes.plot$Module.Color))) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + theme(panel.border = element_rect(size = 1, color = "black"))
p <- p + facet_wrap(~ Module.Color, ncol = 4, scales = "free")
CairoPDF("eigengene_plots", height = 12, width = 16, bg = "transparent")
print(p)
dev.off()

ME.genes.age <- mutate(ME.genes, Age = pheno.import$Age.cont, Combined = pheno.import$Combined) %>%
    gather(Module.Color, Eigengene, -Age, -Combined) 
ME.genes.age$Module.Color %<>% str_replace("ME", "")
ME.genes.age$Combined %<>% str_replace("\\_", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))

p <- ggplot(ME.genes.age, aes(x = Age, y = Eigengene, color = Combined)) + stat_smooth(color = "blue", method = "lm") + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + scale_fill_manual(values = sort(unique(ME.genes.plot$Module.Color)))
p <- p + theme(panel.border = element_rect(size = 1, color = "black"))
p <- p + facet_wrap(~ Module.Color, ncol = 4, scales = "free")
CairoPDF("eigengene_age", height = 10.5, width = 17, bg = "transparent")
print(p)
dev.off()

red.plot <- filter(ME.genes.plot, Module.Color == "red")
p <- ggplot(red.plot, aes(x = Combined, y = Eigengene, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
CairoPDF("eigengene_red", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

red.mu <- as_tibble(posterior.status$MEred) %>% select(mu) %>% unlist
red.estimate <- as_tibble(posterior.status$MEred) %>% select(matches("^Combined")) %>% sweep(1, red.mu) %>% as_tibble
colnames(red.estimate) %<>% str_replace("Combined\\-", "")
red.estimate.plot <- gather(red.estimate, Combined, Estimate) 
red.estimate.plot$Combined %<>% str_replace("_", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))
red.estimate.plot %<>% filter(!(Combined == "duodenum crypts" & Estimate < -0.4)) %>% filter(!(Combined == "duodenum mucosa" & Estimate < -0.15))

p <- ggplot(red.estimate.plot, aes(x = Combined, y = Estimate, fill = Combined)) + geom_violin(scale = "width", trim = FALSE, color = "black", draw_quantiles = c(0.025, 0.5, 0.975)) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank()) + ylab("Eigengene Estimate")
p <- p + theme(panel.border = element_rect(color = "black", size = 1)) + ggtitle("Red Module")
p <- p + theme(plot.title = element_text(hjust = 0.5))
CairoPDF("eigengene_red_estimate", height = 6, width = 7, bg = "transparent")
print(p)
dev.off()

magenta.mu <- as_tibble(posterior.status$MEmagenta) %>% select(mu) %>% unlist
magenta.estimate <- as_tibble(posterior.status$MEmagenta) %>% select(matches("^Combined")) %>% sweep(1, magenta.mu) %>% as_tibble
colnames(magenta.estimate) %<>% str_replace("Combined\\-", "")
magenta.estimate.plot <- gather(magenta.estimate, Combined, Estimate) 
magenta.estimate.plot$Combined %<>% str_replace("\\_", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))
magenta.estimate.plot <- filter(magenta.estimate.plot, !(Combined == "colon mucosa" & Estimate > 0)) %>% 
    filter(!(Combined == "colon crypts" & (Estimate > 0.1 | Estimate < -0.25))) %>% 
    filter(!(Combined == "duodenum mucosa" & (Estimate > 0.15 | Estimate < -0.15))) %>%
    filter(!(Combined == "duodenum crypts" & (Estimate > 0.35 )))

p <- ggplot(magenta.estimate.plot, aes(x = Combined, y = Estimate, fill = Combined)) + geom_violin(scale = "width", trim = FALSE, color = "black", draw_quantiles = c(0.025, 0.5, 0.975)) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank()) + ylab("Eigengene Estimate")
p <- p + theme(panel.border = element_rect(color = "black", size = 1)) + ggtitle("magenta Module")
p <- p + theme(plot.title = element_text(hjust = 0.5))
CairoPDF("eigengene_magenta_estimate", height = 6, width = 7, bg = "transparent")
print(p)
dev.off()

dnam.filter <- ReadRDSgz("../baseline_methylation/save/dnam.filter.rda")
dnam.filter$red <- ME.genes$MEred
dnam.filter$magenta <- ME.genes$MEmagenta

red.difference.cor <- bicor(dnam.filter$red, dnam.filter$Difference)
red.difference.pvalue <- corPvalueStudent(red.difference.cor, nrow(dnam.filter))
p <- ggplot(dnam.filter, aes(x = Difference, y = red, color = Combined)) + stat_smooth(color = "blue", method = "lm") + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(plot.title = element_text(hjust = 0.5), legend.background = element_blank()) + ylab("Eigenpromoter") + xlab("DNAm Age - Chronological Age")
p <- p + ggtitle(str_c("Red Module (", expression(rho), " = ", signif(red.difference.cor,3), ", p < ", signif(red.difference.pvalue, 3), ")"))
CairoPDF("red_difference", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

magenta.difference.cor <- bicor(dnam.filter$magenta, dnam.filter$Difference)
magenta.difference.pvalue <- corPvalueStudent(magenta.difference.cor, nrow(dnam.filter))
p <- ggplot(dnam.filter, aes(x = Difference, y = magenta, color = Combined)) + stat_smooth(color = "blue", method = "lm") + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(plot.title = element_text(hjust = 0.5), legend.background = element_blank()) + ylab("Eigenpromoter") + xlab("DNAm Age - Chronological Age")
p <- p + ggtitle(str_c("Magenta Module (", expression(rho), " = ", signif(magenta.difference.cor,3), ", p < ", signif(magenta.difference.pvalue, 3), ")"))
CairoPDF("magenta_difference", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

##Other plots
#all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% llply(`[`, "y")
#smooth.df <- data.frame(all.smooth)
#colnames(smooth.df) <- names(all.smooth)
#smooth.df$x <- as.factor(1:nrow(smooth.df))
#smooth.plot <- melt(smooth.df, id.vars = "x")

#PCAPlot("all_principal_components", smooth.plot, FALSE, 10, 15)
#PCAPlot("facet_principal_components", smooth.plot, TRUE, 13, 25)

#sample.ids <- factor(rownames(expr.collapse), levels = rownames(expr.collapse))
#colnames(ME.genes) %<>% str_replace("ME", "")
#ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
#expr.data.plot <- data.frame(t(expr.collapse), module.colors)
#split(expr.data.plot, expr.data.plot$module.colors) %>% map(Heatmap, ME.genes.plot)

modules.filter <- filter(final.genetable, nchar(Symbol) > 0 & Gene.Type == "protein_coding") %>% as_tibble
promoters.table <- read.xlsx("../baseline_methylation/promoters.combined.xlsx") 
promoters.filter <- as_tibble(promoters.table) %>% 
    select(Ensembl.ID, duodenum.mucosa_vs._duodenum.crypts.log2FC, duodenum.mucosa_vs._duodenum.crypts.p.val.adj, duodenum.crypts_vs._colon.crypts.log2FC, duodenum.crypts_vs._colon.crypts.p.val.adj) %>%
    filter(duodenum.mucosa_vs._duodenum.crypts.p.val.adj < 0.05 & duodenum.crypts_vs._colon.crypts.p.val.adj < 0.05)
promoters.up <- filter(promoters.filter, duodenum.mucosa_vs._duodenum.crypts.log2FC < 0 & duodenum.crypts_vs._colon.crypts.log2FC > 0)
promoters.down <- filter(promoters.filter, duodenum.mucosa_vs._duodenum.crypts.log2FC > 0 & duodenum.crypts_vs._colon.crypts.log2FC < 0)

promoters.filter2 <- as_tibble(promoters.table) %>% 
    select(Ensembl.ID, duodenum.mucosa_vs._duodenum.crypts.log2FC, duodenum.mucosa_vs._duodenum.crypts.p.val.adj, colon.mucosa_vs._colon.crypts.log2FC, colon.mucosa_vs._colon.crypts.p.val.adj) %>%
    filter(duodenum.mucosa_vs._duodenum.crypts.p.val.adj < 0.05 & colon.mucosa_vs._colon.crypts.p.val.adj < 0.05)
promoters.up2 <- filter(promoters.filter2, duodenum.mucosa_vs._duodenum.crypts.log2FC < 0 & colon.mucosa_vs._colon.crypts.log2FC < 0)
promoters.down2 <- filter(promoters.filter2, duodenum.mucosa_vs._duodenum.crypts.log2FC > 0 & colon.mucosa_vs._colon.crypts.log2FC > 0)

red.down <- filter(modules.filter, Module == "red") %>% inner_join(promoters.down)
Top5Plot("kscaled", red.down, promoters.compare, pheno.import, "Combined", contrasts.levels, "", "top5.red.down")
magenta.down <- filter(modules.filter, Module == "magenta") %>% inner_join(promoters.up)
Top5Plot("kscaled", magenta.down, promoters.compare, pheno.import, "Combined", contrasts.levels, "", "top5.magenta.down")

turquoise.down <- filter(modules.filter, Module == "turquoise") %>% inner_join(promoters.down2)
Top5Plot("kscaled", turquoise.down, promoters.compare, pheno.import, "Combined", contrasts.levels, "", "top5.turquoise.down")
brown.up <- filter(modules.filter, Module == "brown") %>% inner_join(promoters.up2)
Top5Plot("kscaled", brown.up, promoters.compare, pheno.import, "Combined", contrasts.levels, "", "top5.brown.up")

#Enrichr
source("../../code/GO/enrichr.R")
modules.enrichr <- filter(final.genetable, nchar(Symbol) > 0)
modules.submit <- select(modules.enrichr, Ensembl.ID, Symbol, Module)
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016") 
color.names <- unique(module.colors) %>% sort
trap1 <- map(color.names, EnrichrSubmit, modules.submit, enrichr.terms, FALSE)

red.only <- filter(modules.submit, Module == "red")$Symbol
red.gobiol.file <- "./enrichr/red/red_GO_Biological_Process_2015.xlsx"
red.gobiol <- read.xlsx(red.gobiol.file) 
red.gobiol.filter <- FilterEnrichr(red.gobiol) %>% as_tibble
GetKappaCluster(file_path_sans_ext(red.gobiol.file), red.gobiol.filter, red.only)
red.gobiol.final <- slice(red.gobiol.filter, c(2, 5, 7, 13))
red.gobiol.final$Database <- "GO BP"

red.gomole.file <- "./enrichr/red/red_GO_Molecular_Function_2015.xlsx"
red.gomole <- read.xlsx(red.gomole.file) 
red.gomole.filter <- FilterEnrichr(red.gomole) %>% as_tibble
GetKappaCluster(file_path_sans_ext(red.gomole.file), red.gomole.filter, red.only)
red.gomole.final <- slice(red.gomole.filter, c(1, 4))
red.gomole.final$Database <- "GO MF"

red.reactome.file <- "./enrichr/red/red_Reactome_2016.xlsx"
red.reactome <- read.xlsx(red.reactome.file) 
red.reactome.filter <- FilterEnrichr(red.reactome) %>% as_tibble
GetKappaCluster(file_path_sans_ext(red.reactome.file), red.reactome.filter, red.only)
red.reactome.final <- slice(red.reactome.filter, c(3, 4))
red.reactome.final$Database <- "Reactome"

red.kegg.file <- "./enrichr/red/red_KEGG_2016.xlsx"
red.kegg <- read.xlsx(red.kegg.file) 
red.kegg.filter <- FilterEnrichr(red.kegg) %>% as_tibble
GetKappaCluster(file_path_sans_ext(red.kegg.file), red.kegg.filter, red.only)
red.kegg.final <- slice(red.kegg.filter, c(3))
red.kegg.final$Database <- "KEGG"

red.enrichr.final <- rbind(red.gobiol.final, red.gomole.final, red.reactome.final, red.kegg.final)
EnrichrPlot(red.enrichr.final, "red.enrichr", "Red Module")

magenta.only <- filter(modules.submit, Module == "magenta")$Symbol
magenta.gobiol.file <- "./enrichr/magenta/magenta_GO_Biological_Process_2015.xlsx"
magenta.gobiol <- read.xlsx(magenta.gobiol.file) 
magenta.gobiol.filter <- FilterEnrichr(magenta.gobiol) %>% as_tibble
GetKappaCluster(file_path_sans_ext(magenta.gobiol.file), magenta.gobiol.filter, magenta.only)
magenta.gobiol.final <- slice(magenta.gobiol.filter, c(1, 4))
magenta.gobiol.final$Database <- "GO BP"

magenta.only <- filter(modules.submit, Module == "magenta")$Symbol
magenta.gomole.file <- "./enrichr/magenta/magenta_GO_Molecular_Function_2015.xlsx"
magenta.gomole <- read.xlsx(magenta.gomole.file) 
magenta.gomole.filter <- FilterEnrichr(magenta.gomole) %>% as_tibble
GetKappaCluster(file_path_sans_ext(magenta.gomole.file), magenta.gomole.filter, magenta.only)
magenta.gomole.final <- slice(magenta.gomole.filter, c(1))
magenta.gomole.final$Database <- "GO MF"

magenta.reactome.file <- "./enrichr/magenta/magenta_Reactome_2016.xlsx"
magenta.reactome <- read.xlsx(magenta.reactome.file) 
magenta.reactome.filter <- FilterEnrichr(magenta.reactome) %>% as_tibble
GetKappaCluster(file_path_sans_ext(magenta.reactome.file), magenta.reactome.filter, magenta.only)
magenta.reactome.final <- slice(magenta.reactome.filter, c(1))
magenta.reactome.final$Database <- "Reactome"

magenta.kegg.file <- "./enrichr/magenta/magenta_KEGG_2016.xlsx"
magenta.kegg <- read.xlsx(magenta.kegg.file) 
magenta.kegg.filter <- FilterEnrichr(magenta.kegg) %>% as_tibble
GetKappaCluster(file_path_sans_ext(magenta.kegg.file), magenta.kegg.filter, magenta.only)
magenta.kegg.final <- slice(magenta.kegg.filter, c(1))
magenta.kegg.final$Database <- "KEGG"

magenta.enrichr.final <- rbind(magenta.gobiol.final, magenta.gomole.final, magenta.kegg.final, magenta.reactome.final)
EnrichrPlot(magenta.enrichr.final, "magenta.enrichr", "magenta Module")

brown.only <- filter(modules.submit, Module == "brown")$Symbol
brown.gobiol.file <- "./enrichr/brown/brown_GO_Biological_Process_2015.xlsx"
brown.gobiol <- read.xlsx(brown.gobiol.file) 
brown.gobiol.filter <- FilterEnrichr(brown.gobiol) %>% as_tibble
GetKappaCluster(file_path_sans_ext(brown.gobiol.file), brown.gobiol.filter, brown.only)
brown.gobiol.final <- slice(brown.gobiol.filter, c(2, 4, 7, 8))
brown.gobiol.final$Database <- "GO BP"

brown.only <- filter(modules.submit, Module == "brown")$Symbol
brown.gomole.file <- "./enrichr/brown/brown_GO_Molecular_Function_2015.xlsx"
brown.gomole <- read.xlsx(brown.gomole.file) 
brown.gomole.filter <- FilterEnrichr(brown.gomole) %>% as_tibble
GetKappaCluster(file_path_sans_ext(brown.gomole.file), brown.gomole.filter, brown.only)
brown.gomole.final <- slice(brown.gobiol.filter, c(1))
brown.gomole.final$Database <- "GO MF"

brown.reactome.file <- "./enrichr/brown/brown_Reactome_2016.xlsx"
brown.reactome <- read.xlsx(brown.reactome.file) 
brown.reactome.filter <- FilterEnrichr(brown.reactome) %>% as_tibble
GetKappaCluster(file_path_sans_ext(brown.reactome.file), brown.reactome.filter, brown.only)
brown.reactome.final <- slice(brown.reactome.filter, c(1, 2))
brown.reactome.final$Database <- "Reactome"

brown.kegg.file <- "./enrichr/brown/brown_KEGG_2016.xlsx"
brown.kegg <- read.xlsx(brown.kegg.file) 
brown.kegg.filter <- FilterEnrichr(brown.kegg) %>% as_tibble
GetKappaCluster(file_path_sans_ext(brown.kegg.file), brown.kegg.filter, brown.only)
brown.kegg.final <- slice(brown.kegg.filter, c(1, 2))
brown.kegg.final$Database <- "KEGG"

brown.enrichr.final <- rbind(brown.gobiol.final, brown.gomole.final, brown.kegg.final, brown.reactome.final)
EnrichrPlot(brown.enrichr.final, "brown.enrichr", "brown Module", plot.width = 10)

turquoise.only <- filter(modules.submit, Module == "turquoise")$Symbol
turquoise.gobiol.file <- "./enrichr/turquoise/turquoise_GO_Biological_Process_2015.xlsx"
turquoise.gobiol <- read.xlsx(turquoise.gobiol.file) 
turquoise.gobiol.filter <- FilterEnrichr(turquoise.gobiol) %>% as_tibble
GetKappaCluster(file_path_sans_ext(turquoise.gobiol.file), turquoise.gobiol.filter, turquoise.only)
turquoise.gobiol.final <- slice(turquoise.gobiol.filter, c(1, 6, 11, 5))
turquoise.gobiol.final$Database <- "GO BP"

turquoise.only <- filter(modules.submit, Module == "turquoise")$Symbol
turquoise.gomole.file <- "./enrichr/turquoise/turquoise_GO_Molecular_Function_2015.xlsx"
turquoise.gomole <- read.xlsx(turquoise.gomole.file) 
turquoise.gomole.filter <- FilterEnrichr(turquoise.gomole) %>% as_tibble
GetKappaCluster(file_path_sans_ext(turquoise.gomole.file), turquoise.gomole.filter, turquoise.only)
turquoise.gomole.final <- slice(turquoise.gomole.filter, c(1, 4, 5))
turquoise.gomole.final$Database <- "GO MF"

turquoise.reactome.file <- "./enrichr/turquoise/turquoise_Reactome_2016.xlsx"
turquoise.reactome <- read.xlsx(turquoise.reactome.file) 
turquoise.reactome.filter <- FilterEnrichr(turquoise.reactome) %>% as_tibble
GetKappaCluster(file_path_sans_ext(turquoise.reactome.file), turquoise.reactome.filter, turquoise.only)
turquoise.reactome.final <- slice(turquoise.reactome.filter, c(1, 2, 4))
turquoise.reactome.final$Database <- "Reactome"

turquoise.kegg.file <- "./enrichr/turquoise/turquoise_KEGG_2016.xlsx"
turquoise.kegg <- read.xlsx(turquoise.kegg.file) 
turquoise.kegg.filter <- FilterEnrichr(turquoise.kegg) %>% as_tibble
GetKappaCluster(file_path_sans_ext(turquoise.kegg.file), turquoise.kegg.filter, turquoise.only)
turquoise.kegg.final <- slice(turquoise.kegg.filter, c(3))
turquoise.kegg.final$Database <- "KEGG"

turquoise.enrichr.final <- rbind(turquoise.gobiol.final, turquoise.gomole.final, turquoise.kegg.final, turquoise.reactome.final)
EnrichrPlot(turquoise.enrichr.final, "turquoise.enrichr", "turquoise Module", plot.width = 9)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

#Clock genes
clock.transcripts <- read.xlsx("../baseline_methylation/clock.transcripts.xlsx")
clock.coef.table <- read_csv("../baseline_methylation/13059_2013_3156_MOESM23_ESM.csv")
colnames(clock.coef.table)[1] <- "Illumina.ID"
clock.join <- left_join(clock.transcripts, clock.coef.table)
clock.red <- filter(clock.join, nearestGeneSymbol %in% red.only) %>% arrange(desc(abs(CoefficientTraining)))
clock.magenta <- filter(clock.join, nearestGeneSymbol %in% magenta.only) %>% arrange(desc(abs(CoefficientTraining)))
