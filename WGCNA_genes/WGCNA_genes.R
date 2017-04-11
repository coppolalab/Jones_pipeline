library(WGCNA)
library(flashClust)
enableWGCNAThreads()
library(biomaRt)
library(limma)
library(Biobase)

library(PMCMR)
library(BayesFactor)

library(Cairo)
library(igraph)

library(openxlsx)
library(broom)
library(magrittr)
library(stringr)
library(tidyverse)

source("../../code/common_functions.R")

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
    pval.ccdc <- trait.pvals["colon crypts","duodenum crypts"]
    pval.cmdm <- trait.pvals["colon mucosa","duodenum mucosa"]
    pval.dcdm <- trait.pvals["duodenum crypts","duodenum mucosa"]
    pval.cccm <- trait.pvals["colon crypts","colon mucosa"]
    pvals.subset <- p.adjust(c(pval.ccdc, pval.cmdm, pval.dcdm, pval.cccm), method = "fdr")

    trait.medians <- group_by(trait.df, Trait) %>% summarise(median(ME)) %>% data.frame
    colnames(trait.medians)[2] <- "Median"
    rownames(trait.medians) <- trait.medians$Trait
    pval.ccdc <- trait.medians["duodenum crypts","Median"] - trait.medians["colon crypts","Median"] 
    pval.cmdm <- trait.medians["duodenum mucosa","Median"] - trait.medians["colon mucosa","Median"] 
    pval.dcdm <- trait.medians["duodenum mucosa","Median"] - trait.medians["duodenum crypts","Median"] 
    pval.cccm <- trait.medians["colon mucosa","Median"] - trait.medians["colon crypts","Median"] 

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

source("../../code/common_functions.R")
genes.compare <- ReadRDSgz("../baseline_methylation/save/genes.compare.combat.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

#sft for combat only = 10
sft.bicor <- pickSoftThreshold(t(genes.compare), powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
sft.df <- sft.bicor$fitIndices
SaveRDSgz(sft.bicor, file = "./save/sft.bicor.rda")

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

softPower <- 18
adjacency.expr <- adjacency(t(genes.compare), power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
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
ME.list <- moduleEigengenes(t(genes.compare), colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
SaveRDSgz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50_tree", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(t(genes.compare), dynamic.colors, cutHeight = ME.dissimilarity.threshold, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 3) 
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

pheno.import <- ReadRDSgz("../baseline_methylation/save/pheno.compare.rda")
pheno.import$Combined %<>% str_c(pheno.import$Tissue, pheno.import$Cell.Type, sep = " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))
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

genes.bm.table.vega <- getBM(attributes = c('ensg', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ensg', values = as.character(gene.info$Ensembl.ID), mart = vega)
colnames(genes.bm.table.vega) <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

gene.info.annot <- left_join(gene.info, genes.bm.table) %>% left_join(genes.bm.table.vega) 
gene.module.membership <- data.frame(bicor(t(genes.compare), ME.genes, maxPOutliers = 0.05)) %>% signif(3)
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(t(genes.compare)))) %>% signif(3)
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

black.mu <- data.frame(posterior.status$MEblack) %>% select(mu) %>% unlist
black.estimate <- data.frame(posterior.status$MEblack) %>% select(matches("^Combined")) %>% sweep(1, black.mu)
colnames(black.estimate) %<>% str_replace("Combined\\.", "")
black.estimate.plot <- gather(black.estimate, Combined, Estimate) %>% as_tibble
black.estimate.plot$Combined %<>% str_replace("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))
black.estimate.plot %<>% filter(!(Combined == "duodenum crypts" & Estimate < -0.4)) %>% filter(!(Combined == "colon mucosa" & Estimate < 0))

p <- ggplot(black.estimate.plot, aes(x = Combined, y = Estimate, fill = Combined)) + geom_violin(color = "black", draw_quantiles = c(0.025, 0.5, 0.975)) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank()) +  ylab("Eigengene Estimate")
p <- p + theme(panel.border = element_rect(color = "black", size = 1)) + ggtitle("Black Module")
p <- p + theme(plot.title = element_text(hjust = 0.5))
CairoPDF("eigengene_black_estimate", height = 6, width = 7, bg = "transparent")
print(p)
dev.off()

magenta.mu <- data.frame(posterior.status$MEmagenta) %>% select(mu) %>% unlist
magenta.estimate <- data.frame(posterior.status$MEmagenta) %>% select(matches("^Combined")) %>% sweep(1, magenta.mu)
colnames(magenta.estimate) %<>% str_replace("Combined\\.", "")
magenta.estimate.plot <- gather(magenta.estimate, Combined, Estimate) 
magenta.estimate.plot$Combined %<>% str_replace("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))
magenta.estimate.plot %<>% filter(!(Combined == "duodenum crypts" & Estimate < -0)) %>% 
    filter(!(Combined == "colon mucosa" & Estimate > 0)) %>%
    filter(!(Combined == "colon crypts" & Estimate > 0.1))

p <- ggplot(magenta.estimate.plot, aes(x = Combined, y = Estimate, fill = Combined)) + geom_violin(color = "black", draw_quantiles = c(0.025, 0.5, 0.975)) 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank()) + ylab("Eigengene Estimate")
p <- p + theme(panel.border = element_rect(color = "black", size = 1)) + ggtitle("magenta Module")
p <- p + theme(plot.title = element_text(hjust = 0.5)) 
CairoPDF("eigengene_magenta_estimate", height = 6, width = 7, bg = "transparent")
print(p)
dev.off()

dnam.filter <- ReadRDSgz("../baseline_methylation/save/dnam.filter.rda")
dnam.filter$black <- ME.genes$MEblack
dnam.filter$magenta <- ME.genes$MEmagenta

black.difference.cor <- bicor(dnam.filter$black, dnam.filter$Difference)
black.difference.pvalue <- corPvalueStudent(black.difference.cor, nrow(dnam.filter))
p <- ggplot(dnam.filter, aes(x = Difference, y = black, color = Combined)) + stat_smooth(color = "blue", method = "lm") + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(plot.title = element_text(hjust = 0.5), legend.background = element_blank()) + ylab("Eigenpromoter") + xlab("DNAm Age - Chronological Age")
p <- p + ggtitle(str_c("black Module (", expression(rho), " = ", signif(black.difference.cor,3), ", p < ", signif(black.difference.pvalue, 3), ")"))
CairoPDF("black_difference", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

magenta.difference.cor <- bicor(dnam.filter$magenta, dnam.filter$Difference)
magenta.difference.pvalue <- corPvalueStudent(magenta.difference.cor, nrow(dnam.filter))
p <- ggplot(dnam.filter, aes(x = Difference, y = magenta, color = Combined)) + stat_smooth(color = "blue", method = "lm") + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(plot.title = element_text(hjust = 0.5), legend.background = element_blank()) + ylab("Eigenpromoter") + xlab("DNAm Age - Chronological Age")
p <- p + ggtitle(str_c("magenta Module (", expression(rho), " = ", signif(magenta.difference.cor,3), ", p < ", signif(magenta.difference.pvalue, 3), ")"))
CairoPDF("magenta_difference", height = 6, width = 8, bg = "transparent")
print(p)
dev.off()

#Other Plots
modules.filter <- filter(final.genetable, nchar(Symbol) > 0 & Gene.Type == "protein_coding") %>% as_tibble
genes.table <- read.xlsx("../baseline_methylation/genes.combined.xlsx") 
genes.filter <- as_tibble(genes.table) %>% 
    select(Ensembl.ID, duodenum.mucosa_vs._duodenum.crypts.log2FC, duodenum.mucosa_vs._duodenum.crypts.p.val.adj, duodenum.crypts_vs._colon.crypts.log2FC, duodenum.crypts_vs._colon.crypts.p.val.adj) %>%
    filter(duodenum.mucosa_vs._duodenum.crypts.p.val.adj < 0.05 & duodenum.crypts_vs._colon.crypts.p.val.adj < 0.05)
genes.up <- filter(genes.filter, duodenum.mucosa_vs._duodenum.crypts.log2FC < 0 & duodenum.crypts_vs._colon.crypts.log2FC > 0)
genes.down <- filter(genes.filter, duodenum.mucosa_vs._duodenum.crypts.log2FC > 0 & duodenum.crypts_vs._colon.crypts.log2FC < 0)

genes.filter2 <- as_tibble(genes.table) %>% 
    select(Ensembl.ID, duodenum.mucosa_vs._duodenum.crypts.log2FC, duodenum.mucosa_vs._duodenum.crypts.p.val.adj, colon.mucosa_vs._colon.crypts.log2FC, colon.mucosa_vs._colon.crypts.p.val.adj) %>%
    filter(duodenum.mucosa_vs._duodenum.crypts.p.val.adj < 0.05 & colon.mucosa_vs._colon.crypts.p.val.adj < 0.05)
genes.up2 <- filter(genes.filter2, duodenum.mucosa_vs._duodenum.crypts.log2FC < 0 & colon.mucosa_vs._colon.crypts.log2FC < 0)
genes.down2 <- filter(genes.filter2, duodenum.mucosa_vs._duodenum.crypts.log2FC > 0 & colon.mucosa_vs._colon.crypts.log2FC > 0)

black.down <- filter(modules.filter, Module == "black") %>% inner_join(genes.down)
Top5Plot("kscaled", black.down, genes.compare, pheno.import, "Combined", contrasts.levels, "", "top5.black.down")
magenta.down <- filter(modules.filter, Module == "magenta") %>% inner_join(genes.up)
Top5Plot("kscaled", magenta.down, genes.compare, pheno.import, "Combined", contrasts.levels, "", "top5.magenta.down")

turquoise.up <- filter(modules.filter, Module == "turquoise") %>% inner_join(genes.up2)
Top5Plot("kscaled", turquoise.up, genes.compare, pheno.import, "Combined", contrasts.levels, "", "top5.turquoise.up")
blue.down <- filter(modules.filter, Module == "blue") %>% inner_join(genes.down2)
Top5Plot("kscaled", blue.down, genes.compare, pheno.import, "Combined", contrasts.levels, "", "top5.blue.down")

source("../../code/GO/enrichr.R")
modules.submit <- filter(final.genetable, nchar(Symbol) > 0) %>% as_tibble

enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016") 
color.names <- unique(module.colors) %>% sort
trap1 <- map(color.names, EnrichrSubmit, modules.submit, enrichr.terms, FALSE)

black.only <- filter(modules.submit, Module == "black")$Symbol
black.gobiol.file <- "./enrichr/black/black_GO_Biological_Process_2015.xlsx"
black.gobiol <- read.xlsx(black.gobiol.file) 
black.gobiol.filter <- FilterEnrichr(black.gobiol) %>% as_tibble
GetKappaCluster(file_path_sans_ext(black.gobiol.file), black.gobiol.filter, black.only)
black.gobiol.final <- slice(black.gobiol.filter, c(3, 19, 20, 4))
black.gobiol.final$Database <- "GO BP"

black.gomole.file <- "./enrichr/black/black_GO_Molecular_Function_2015.xlsx"
black.gomole <- read.xlsx(black.gomole.file) 
black.gomole.filter <- FilterEnrichr(black.gomole) %>% as_tibble
GetKappaCluster(file_path_sans_ext(black.gomole.file), black.gomole.filter, black.only)
black.gomole.final <- slice(black.gomole.filter, c(9, 14, 15))
black.gomole.final$Database <- "GO MF"

black.reactome.file <- "./enrichr/black/black_Reactome_2016.xlsx"
black.reactome <- read.xlsx(black.reactome.file) 
black.reactome.filter <- FilterEnrichr(black.reactome) %>% as_tibble
GetKappaCluster(file_path_sans_ext(black.reactome.file), black.reactome.filter, black.only)
black.reactome.final <- slice(black.reactome.filter, c(2, 3, 5))
black.reactome.final$Database <- "Reactome"

black.kegg.file <- "./enrichr/black/black_KEGG_2016.xlsx"
black.kegg <- read.xlsx(black.kegg.file) 
black.kegg.filter <- FilterEnrichr(black.kegg) %>% as_tibble
GetKappaCluster(file_path_sans_ext(black.kegg.file), black.kegg.filter, black.only)
black.kegg.final <- slice(black.kegg.filter, c(3))
black.kegg.final$Database <- "Reactome"

black.enrichr.final <- rbind(black.gobiol.final, black.gomole.final, black.reactome.final, black.kegg.final)
EnrichrPlot(black.enrichr.final, "black.enrichr", "Black Module")

magenta.only <- filter(modules.submit, Module == "magenta")$Symbol
magenta.gobiol.file <- "./enrichr/magenta/magenta_GO_Biological_Process_2015.xlsx"
magenta.gobiol <- read.xlsx(magenta.gobiol.file) 
magenta.gobiol.filter <- FilterEnrichr(magenta.gobiol) %>% as_tibble
GetKappaCluster(file_path_sans_ext(magenta.gobiol.file), magenta.gobiol.filter, magenta.only)
magenta.gobiol.final <- slice(magenta.gobiol, c(57, 26, 22))
magenta.gobiol.final$Database <- "GO BP"

magenta.gomole.file <- "./enrichr/magenta/magenta_GO_Molecular_Function_2015.xlsx"
magenta.gomole <- read.xlsx(magenta.gomole.file) 
magenta.gomole.filter <- FilterEnrichr(magenta.gomole) %>% as_tibble
GetKappaCluster(file_path_sans_ext(magenta.gomole.file), magenta.gomole.filter, magenta.only)
magenta.gomole.final <- slice(magenta.gomole, c(12))
magenta.gomole.final$Database <- "GO MF"

magenta.reactome.file <- "./enrichr/magenta/magenta_Reactome_2016.xlsx"
magenta.reactome <- read.xlsx(magenta.reactome.file) 
magenta.reactome.filter <- FilterEnrichr(magenta.reactome) %>% as_tibble
GetKappaCluster(file_path_sans_ext(magenta.reactome.file), magenta.reactome.filter, magenta.only)
magenta.reactome.final <- slice(magenta.reactome, c(4))
magenta.reactome.final$Database <- "Reactome"

magenta.kegg.file <- "./enrichr/magenta/magenta_KEGG_2016.xlsx"
magenta.kegg <- read.xlsx(magenta.kegg.file) 
magenta.kegg.filter <- FilterEnrichr(magenta.kegg) %>% as_tibble
GetKappaCluster(file_path_sans_ext(magenta.kegg.file), magenta.kegg.filter, magenta.only)
magenta.kegg.final <- slice(magenta.kegg, c(3))
magenta.kegg.final$Database <- "KEGG"

magenta.enrichr.final <- rbind(magenta.gobiol.final, magenta.gomole.final, magenta.reactome.final, magenta.kegg.final)
EnrichrPlot(magenta.enrichr.final, "magenta.enrichr", "Magenta Module")

blue.only <- filter(modules.submit, Module == "blue")$Symbol
blue.gobiol.file <- "./enrichr/blue/blue_GO_Biological_Process_2015.xlsx"
blue.gobiol <- read.xlsx(blue.gobiol.file) 
blue.gobiol.filter <- FilterEnrichr(blue.gobiol) %>% as_tibble
GetKappaCluster(file_path_sans_ext(blue.gobiol.file), blue.gobiol.filter, blue.only)
blue.gobiol.final <- slice(blue.gobiol.filter, c(1, 4, 6, 8))
blue.gobiol.final$Database <- "GO BP"

blue.gomole.file <- "./enrichr/blue/blue_GO_Molecular_Function_2015.xlsx"
blue.gomole <- read.xlsx(blue.gomole.file) 
blue.gomole.filter <- FilterEnrichr(blue.gomole) %>% as_tibble
GetKappaCluster(file_path_sans_ext(blue.gomole.file), blue.gomole.filter, blue.only)
blue.gomole.final <- slice(blue.gomole.filter, c(1, 3))
blue.gomole.final$Database <- "GO MF"

blue.reactome.file <- "./enrichr/blue/blue_Reactome_2016.xlsx"
blue.reactome <- read.xlsx(blue.reactome.file) 
blue.reactome.filter <- FilterEnrichr(blue.reactome) %>% as_tibble
GetKappaCluster(file_path_sans_ext(blue.reactome.file), blue.reactome.filter, blue.only)
blue.reactome.final <- slice(blue.reactome.filter, c(1))
blue.reactome.final$Database <- "Reactome"

blue.kegg.file <- "./enrichr/blue/blue_KEGG_2016.xlsx"
blue.kegg <- read.xlsx(blue.kegg.file) 
blue.kegg.filter <- FilterEnrichr(blue.kegg) %>% as_tibble
GetKappaCluster(file_path_sans_ext(blue.kegg.file), blue.kegg.filter, blue.only)
blue.kegg.final <- slice(blue.kegg.filter, c(1, 2))
blue.kegg.final$Database <- "Reactome"

blue.enrichr.final <- rbind(blue.gobiol.final, blue.gomole.final, blue.reactome.final, blue.kegg.final)
EnrichrPlot(blue.enrichr.final, "blue.enrichr", "Blue Module")

turquoise.only <- filter(modules.submit, Module == "turquoise")$Symbol
turquoise.gobiol.file <- "./enrichr/turquoise/turquoise_GO_Biological_Process_2015.xlsx"
turquoise.gobiol <- read.xlsx(turquoise.gobiol.file) 
turquoise.gobiol.filter <- FilterEnrichr(turquoise.gobiol) %>% as_tibble
GetKappaCluster(file_path_sans_ext(turquoise.gobiol.file), turquoise.gobiol.filter, turquoise.only)
turquoise.gobiol.final <- slice(turquoise.gobiol, c(1, 5, 11, 12))
turquoise.gobiol.final$Database <- "GO BP"

turquoise.gomole.file <- "./enrichr/turquoise/turquoise_GO_Molecular_Function_2015.xlsx"
turquoise.gomole <- read.xlsx(turquoise.gomole.file) 
turquoise.gomole.filter <- FilterEnrichr(turquoise.gomole) %>% as_tibble
GetKappaCluster(file_path_sans_ext(turquoise.gomole.file), turquoise.gomole.filter, turquoise.only)
turquoise.gomole.final <- slice(turquoise.gomole, c(1))
turquoise.gomole.final$Database <- "GO MF"

turquoise.reactome.file <- "./enrichr/turquoise/turquoise_Reactome_2016.xlsx"
turquoise.reactome <- read.xlsx(turquoise.reactome.file) 
turquoise.reactome.filter <- FilterEnrichr(turquoise.reactome) %>% as_tibble
GetKappaCluster(file_path_sans_ext(turquoise.reactome.file), turquoise.reactome.filter, turquoise.only)
turquoise.reactome.final <- slice(turquoise.reactome, c(1))
turquoise.reactome.final$Database <- "Reactome"

turquoise.kegg.file <- "./enrichr/turquoise/turquoise_KEGG_2016.xlsx"
turquoise.kegg <- read.xlsx(turquoise.kegg.file) 
turquoise.kegg.filter <- FilterEnrichr(turquoise.kegg) %>% as_tibble
GetKappaCluster(file_path_sans_ext(turquoise.kegg.file), turquoise.kegg.filter, turquoise.only)
turquoise.kegg.final <- slice(turquoise.kegg, c(1))
turquoise.kegg.final$Database <- "KEGG"

turquoise.enrichr.final <- rbind(turquoise.gobiol.final, turquoise.gomole.final, turquoise.reactome.final, turquoise.kegg.final)
EnrichrPlot(turquoise.enrichr.final, "turquoise.enrichr", "turquoise Module")

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

