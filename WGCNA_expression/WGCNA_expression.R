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
    filter.string <- str_c("Module == '", module.color, "'")
    module.select <- filter_(module.table, filter.string) %>% slice(1:15)
    adjacency.select <- adjacency.expr[module.select$Symbol,module.select$Symbol] 
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

genes.rmcov <- ReadRDSgz("../differential_expression/save/rmcov.lumi.rda")

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

sft <- pickSoftThreshold(t(exprs(genes.rmcov)), powerVector = powers, verbose = 5, networkType = "signed")
sft.bicor <- pickSoftThreshold(t(exprs(genes.rmcov)), powerVector = powers, verbose = 5, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), networkType = "signed")
sft.spearman <- pickSoftThreshold(t(exprs(genes.rmcov)), powerVector = powers, verbose = 5, corFnc = cor, corOptions = list(method = "spearman"), networkType = "signed")
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

softPower <- 9
adjacency.expr <- adjacency(t(exprs(genes.rmcov)), power = softPower, type = "signed", corFnc = "bicor", corOptions = "maxPOutliers = 0.05")
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
ME.list <- moduleEigengenes(t(exprs(genes.rmcov)), colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - bicor(ME.genes, maxPOutliers = 0.05)
METree <- flashClust(as.dist(MEDiss), method = "average")
SaveRDSgz(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./module_eigengene_clustering_min50_tree", height = 10, width = 15, bg = "transparent")
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.2
merge.all <- mergeCloseModules(t(exprs(genes.rmcov)), dynamic.colors, cutHeight = ME.dissimilarity.threshold, corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 3) 
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

pheno.import <- ReadRDSgz("../differential_expression/save/rmcov.pheno.rda")
pheno.import$Combined %<>% factor %>% droplevels
anova.status <- map_dbl(ME.genes, EigengeneANOVA, pheno.import$Combined) %>% p.adjust("fdr") %>% signif(3)
bayes.status <- map(ME.genes, EigengeneBayes, pheno.import$Combined) 
bf.status <- map(bayes.status, extractBF) %>% map_dbl(extract2, "bf")
posterior.status <- map(bayes.status, posterior, iterations = 100000) 
contrasts.vector <- c("duodenum_crypts_vs._colon_crypts")

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
gene.info <- data.frame(Symbol = rownames(all.degrees), Module = module.colors, all.degrees)
gene.info$kscaled <- by(gene.info, gene.info$Module, select, kWithin) %>% map(function(kWithin) kWithin / max(kWithin)) %>% reduce(c)
gene.info[,3:ncol(gene.info)] <- signif(gene.info[,3:ncol(gene.info)], 3)

#Annotate gene table
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes.bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = featureNames(genes.rmcov), mart = ensembl)
genes.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(genes.bm.table) <- c("Symbol", "Description")

gene.info.annot <- left_join(gene.info, genes.bm.table) 
gene.module.membership <- data.frame(bicor(t(exprs(genes.rmcov)), ME.genes, maxPOutliers = 0.05)) %>% signif(3)
module.membership.pvalue <- data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(t(exprs(genes.rmcov))))) %>% signif(3)
colnames(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")
gene.module.membership$Symbol <- rownames(gene.module.membership)
module.membership.pvalue$Symbol <- rownames(module.membership.pvalue)

final.genetable <- left_join(gene.info.annot, gene.module.membership) %>% 
    left_join(module.membership.pvalue, by = "Symbol") %>% 
    dplyr::select(Symbol, Description, Module, kTotal:kscaled, matches("MM.*")) %>%
    arrange(Module, desc(kscaled))
Workbook(final.genetable, "./final_genetable.xlsx")

#Plots
ME.genes.plot <- mutate(ME.genes, Combined = pheno.import$Combined) %>%
    gather(Module.Color, Eigengene, -Combined) 
ME.genes.plot$Module.Color %<>% str_replace("ME", "")
ME.genes.plot$Combined %<>% str_replace("\\_", " ") %>% factor(levels = c("duodenum crypts", "colon crypts"))

p <- ggplot(ME.genes.plot, aes(x = Combined, y = Eigengene, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
p <- p + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
p <- p + scale_fill_manual(values = sort(unique(ME.genes.plot$Module.Color)))
p <- p + theme(panel.border = element_rect(size = 1, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + facet_wrap(~ Module.Color, ncol = 5, scales = "free")
CairoPDF("eigengene_plots", height = 8, width = 16, bg = "transparent")
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

PlotHubGenes("blue", final.genetable, adjacency.expr)
PlotHubGenes("brown", final.genetable, adjacency.expr)
PlotHubGenes("turquoise", final.genetable, adjacency.expr)

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
modules.submit <- select(final.genetable, Symbol, Module)
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016", "Human_Gene_Atlas", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up") 
color.names <- unique(module.colors) %>% sort
map(color.names, EnrichrSubmit, modules.submit, enrichr.terms, FALSE)

blue.only <- filter(modules.submit, Module == "blue")$Symbol
blue.gobiol.file <- "./enrichr/blue/blue_GO_Biological_Process_2015.xlsx"
blue.gobiol <- read.xlsx(blue.gobiol.file) 
blue.gobiol.filter <- FilterEnrichr(blue.gobiol)
GetKappaCluster(file_path_sans_ext(blue.gobiol.file), blue.gobiol.filter, blue.only)
blue.gobiol.final <- slice(blue.gobiol.filter, c(21, 13, 18, 6))
blue.gobiol.final$Database <- "GO Biological Process"

blue.gomole.file <- "./enrichr/blue/blue_GO_Molecular_Function_2015.xlsx"
blue.gomole <- read.xlsx(blue.gomole.file) 
blue.gomole.filter <- FilterEnrichr(blue.gomole)
GetKappaCluster(file_path_sans_ext(blue.gomole.file), blue.gomole.filter, blue.only)
blue.gomole.final <- slice(blue.gomole.filter, c(2, 3))
blue.gomole.final$Database <- "GO Molecular Function"

blue.reactome.file <- "./enrichr/blue/blue_Reactome_2016.xlsx"
blue.reactome <- read.xlsx(blue.reactome.file) 
blue.reactome.filter <- FilterEnrichr(blue.reactome)
GetKappaCluster(file_path_sans_ext(blue.reactome.file), blue.reactome.filter, blue.only)
blue.reactome.final <- slice(blue.reactome.filter, c(4, 29))
blue.reactome.final$Database <- "Reactome"

blue.kegg.file <- "./enrichr/blue/blue_KEGG_2016.xlsx"
blue.kegg <- read.xlsx(blue.kegg.file) 
blue.kegg.filter <- FilterEnrichr(blue.kegg)
GetKappaCluster(file_path_sans_ext(blue.kegg.file), blue.kegg.filter, blue.only)
blue.kegg.final <- slice(blue.kegg.filter, c(3))
blue.kegg.final$Database <- "KEGG"

blue.enrichr.final <- rbind(blue.gobiol.final, blue.gomole.final, blue.reactome.final, blue.kegg.final)
EnrichrPlot(blue.enrichr.final, "blue.enrichr", "Blue Module", plot.width = 10)

brown.only <- filter(modules.submit, Module == "brown")$Symbol
brown.gobiol.file <- "./enrichr/brown/brown_GO_Biological_Process_2015.xlsx"
brown.gobiol <- read.xlsx(brown.gobiol.file) 
brown.gobiol.filter <- FilterEnrichr(brown.gobiol)
GetKappaCluster(file_path_sans_ext(brown.gobiol.file), brown.gobiol.filter, brown.only)
brown.gobiol.final <- slice(brown.gobiol.filter, c(4, 2, 10, 7, 11))
brown.gobiol.final$Database <- "GO Biological Process"

brown.gomole.file <- "./enrichr/brown/brown_GO_Molecular_Function_2015.xlsx"
brown.gomole <- read.xlsx(brown.gomole.file) 
brown.gomole.filter <- FilterEnrichr(brown.gomole)
GetKappaCluster(file_path_sans_ext(brown.gomole.file), brown.gomole.filter, brown.only)
brown.gomole.final <- slice(brown.gomole.filter, c(2, 9))
brown.gomole.final$Database <- "GO Molecular Function"

brown.reactome.file <- "./enrichr/brown/brown_Reactome_2016.xlsx"
brown.reactome <- read.xlsx(brown.reactome.file) 
brown.reactome.filter <- FilterEnrichr(brown.reactome)
GetKappaCluster(file_path_sans_ext(brown.reactome.file), brown.reactome.filter, brown.only)
brown.reactome.final <- slice(brown.reactome.filter, c(3, 5))
brown.reactome.final$Database <- "Reactome"

brown.kegg.file <- "./enrichr/brown/brown_KEGG_2016.xlsx"
brown.kegg <- read.xlsx(brown.kegg.file) 
brown.kegg.filter <- FilterEnrichr(brown.kegg)
GetKappaCluster(file_path_sans_ext(brown.kegg.file), brown.kegg.filter, brown.only)
brown.kegg.final <- slice(brown.kegg.filter, c(1))
brown.kegg.final$Database <- "KEGG"

brown.enrichr.final <- rbind(brown.gobiol.final, brown.kegg.final, brown.reactome.final, brown.gomole.final)
EnrichrPlot(brown.enrichr.final, "brown.enrichr", "brown Module", plot.width = 9)

turquoise.only <- filter(modules.submit, Module == "turquoise")$Symbol
turquoise.gobiol.file <- "./enrichr/turquoise/turquoise_GO_Biological_Process_2015.xlsx"
turquoise.gobiol <- read.xlsx(turquoise.gobiol.file) 
turquoise.gobiol.filter <- FilterEnrichr(turquoise.gobiol)
GetKappaCluster(file_path_sans_ext(turquoise.gobiol.file), turquoise.gobiol.filter, turquoise.only)
turquoise.gobiol.final <- slice(turquoise.gobiol.filter, c(1, 3, 7, 8, 15))
turquoise.gobiol.final$Database <- "GO Biological Process"

turquoise.gomole.file <- "./enrichr/turquoise/turquoise_GO_Molecular_Function_2015.xlsx"
turquoise.gomole <- read.xlsx(turquoise.gomole.file) 
turquoise.gomole.filter <- FilterEnrichr(turquoise.gomole)
GetKappaCluster(file_path_sans_ext(turquoise.gomole.file), turquoise.gomole.filter, turquoise.only)
turquoise.gomole.final <- slice(turquoise.gomole.filter, c(2, 3))
turquoise.gomole.final$Database <- "GO Molecular Function"

turquoise.reactome.file <- "./enrichr/turquoise/turquoise_Reactome_2016.xlsx"
turquoise.reactome <- read.xlsx(turquoise.reactome.file) 
turquoise.reactome.filter <- FilterEnrichr(turquoise.reactome)
GetKappaCluster(file_path_sans_ext(turquoise.reactome.file), turquoise.reactome.filter, turquoise.only)
turquoise.reactome.final <- slice(turquoise.reactome.filter, c(1))
turquoise.reactome.final$Database <- "Reactome"

turquoise.kegg.file <- "./enrichr/turquoise/turquoise_KEGG_2016.xlsx"
turquoise.kegg <- read.xlsx(turquoise.kegg.file) 
turquoise.kegg.filter <- FilterEnrichr(turquoise.kegg)
GetKappaCluster(file_path_sans_ext(turquoise.kegg.file), turquoise.kegg.filter, turquoise.only)
turquoise.kegg.final <- slice(turquoise.kegg.filter, c(1))
turquoise.kegg.final$Database <- "KEGG"

turquoise.enrichr.final <- rbind(turquoise.gobiol.final, turquoise.kegg.final, turquoise.reactome.final, turquoise.gomole.final)
EnrichrPlot(turquoise.enrichr.final, "turquoise.enrichr", "turquoise Module", plot.width = 8)

test <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(test) <- ls()
unlist(test) %>% sort

