#For DE
library(lumi)
library(lumiHumanIDMapping)
library(lumiHumanAll.db)
library(limma)
library(annotate)
library(biomaRt)
library(siggenes)
library(BayesFactor)

library(Cairo)
library(openxlsx)
library(WGCNA)

library(tools)
library(broom)
library(magrittr)
library(rlist)
library(stringr)
library(tidyverse)

#Boxplot function
BoxPlot <- function(filename, lumi.object, colorscheme, maintext, ylabtext) {
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Status = sampleNames(lumi.object), Batch = lumi.object$Batch)
    dataset.m <- gather(dataset.addvars, id = c("Sample.Status", "Batch"))
    p <- ggplot(dataset.m, aes(x = Sample.Status, y = value, fill = factor(Batch))) + geom_boxplot() + theme_bw()
    p <- p + scale_fill_manual(values = colorscheme)
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5))     
    p <- p + ggtitle(maintext) + ylab(ylabtext) + xlab("Sample") + theme(axis.text.x = element_text(size = 3))
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(filename = filename, plot = p, family = "Noto Sans", width = 20 , height = 8)
}

#Histogram
Histogram <- function(filename, lumi.object) {
    expr.df <- exprs(lumi.object) %>% t %>% data.frame
    dataset.addvars <- mutate(expr.df, Sample.Name = sampleNames(lumi.object), Combined = lumi.object$Combined)
    dataset.m <- gather(dataset.addvars, nuID, Expression, -Sample.Name, -Combined)

    p <- ggplot(dataset.m, aes(Expression, group = Sample.Name, col = factor(Combined))) + geom_density() + theme_bw()
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle("Histogram of VST Expression") + ylab("Density") + xlab("VST Expression") 
    CairoPDF(filename, height = 5, width = 9)
    print(p)
    dev.off()
}

#Samplewise connectivity plot
ConnectivityPlot <- function(filename, dataset, maintitle) {
    norm.adj <- (0.5 + 0.5 * bicor(exprs(dataset)))
    colnames(norm.adj) <- dataset$Slide.ID
    rownames(norm.adj) <- dataset$Slide.ID
    net.summary <- fundamentalNetworkConcepts(norm.adj)
    net.connectivity <- net.summary$Connectivity
    connectivity.zscore <- (net.connectivity - mean(net.connectivity)) / sqrt(var(net.connectivity))

    connectivity.plot <- data.frame(Slide.ID = names(connectivity.zscore), Z.score = connectivity.zscore, Sample.Num = 1:length(connectivity.zscore))
    p <- ggplot(connectivity.plot, aes(x = Sample.Num, y = Z.score, label = Slide.ID) )
    p <- p + geom_text(size = 4, colour = "red")
    p <- p + geom_hline(aes(yintercept = -2)) + geom_hline(yintercept = -3) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle(maintitle)
    CairoPDF(filename, width = 10, height = 10)
    print(p)
    dev.off()
    connectivity.zscore
}

#MDS function 
MDSPlot <- function(filename, dataset, targetset, colorscheme = "none", variablename) {
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$Slide.ID, factor(targetset[[variablename]]))
    colnames(target.data) <- c("Slide.ID", variablename)
    colnames(dataset.plot) <- c("Slide.ID", "Component.1", "Component.2")
    dataset.plot <- merge(dataset.plot, target.data)
    dataset.plot[[variablename]] %<>% str_replace_all("_", " ")
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(plot.title = element_text(hjust = 0.5), plot.background = element_blank(), legend.background = element_blank())
    if (colorscheme != "none") {
        p <- p + scale_color_manual(values = colorscheme) 
    }
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle(variablename)
    CairoPDF(file = filename, height = 6, width = 8, bg= "transparent")
    print(p)
    dev.off()
}

#Make Excel spreadsheet
DEWorkbook <- function(de.table, filename) { 
    pval.cols <- colnames(de.table) %>% str_detect("P.Value") %>% which
    adj.pval.cols <- colnames(de.table) %>% str_detect("adj.P.Val") %>% which
    coef.cols <- colnames(de.table) %>% str_detect("logFC") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = de.table)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(de.table), rule = "<0.005", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = adj.pval.cols, rows = 1:nrow(de.table), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(de.table), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = 3:12, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 11)
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#Plot number of genes at each threshold
DecidePlot <- function(decide.plot, filename, width.plot = 6, height.plot = 7) {
    decide.ggplot <- ggplot() + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = Contrast, y = Num.Genes), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    decide.ggplot <- decide.ggplot + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Contrast, y = Num.Genes, ymax = max(Num.Genes) + 110, label = Num.Genes), hjust = -0.3, position = position_dodge(width = 1))
    decide.ggplot <- decide.ggplot + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = Contrast, y = Num.Genes), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    decide.ggplot <- decide.ggplot + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Contrast, y = Num.Genes, ymax = min(Num.Genes) - 110, label = abs(Num.Genes)), hjust = 1.3, position = position_dodge(width = 1))
    decide.ggplot <- decide.ggplot + facet_grid(Test + Num ~ .) 
    decide.ggplot <- decide.ggplot + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    decide.ggplot <- decide.ggplot + theme(axis.title.y = element_blank(), axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0))
    CairoPDF(filename, width = width.plot, height = height.plot)
    print(decide.ggplot)
    dev.off()
}

#Volcano plot
VolcanoPlot <- function(top.table, filename, pval.column = "P.Value", log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Log.Pvalue") {
    top.table$Log.Pvalue <- -log10(top.table[[pval.column]])
    p <- ggplot(top.table, aes_string(x = log.column, y = "Log.Pvalue")) + geom_point(aes(color = Significant))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(legend.position = "none", plot.background = element_blank())
    p <- p + xlab(xlabel) + ylab(ylabel)
    CairoPDF(filename, bg = "transparent")
    print(p)
    dev.off()
}

BayesPlot <- function(siggene, filename, threshold, log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Posterior Probability") {
    p <- ggplot(siggene, aes_string(x = log.column, y = "Posterior")) + geom_point(aes(color = Significant))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(legend.position = "none", plot.background = element_blank())
    p <- p + xlab(xlabel) + ylab(ylabel)
    CairoPDF(filename, width = 4, height = 4, bg = "transparent")
    print(p)
    dev.off()
}

GeneBoxplot <- function(lumi.object, gene.symbol, show.label = FALSE) {
    gene.expr <- as.vector(exprs(lumi.object[gene.symbol,]))
    gene.df <- data.frame(Combined = lumi.object$Combined, Expression = gene.expr)
    gene.df$Combined %<>% factor(levels = c("duodenum_crypts", "colon_crypts"))

    p <- ggplot(gene.df, aes(x = Combined, y = Expression, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter() + theme_bw()
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p +  theme(plot.background = element_blank(), panel.border = element_rect(color = "black", size = 1))
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    if (show.label == TRUE) {
        p <- p + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(gene.symbol)
    }
    CairoPDF(str_c(gene.symbol, ".pdf"), width = 4, height = 3, bg = "transparent")
    print(p)
    dev.off()
}

#GO functions
EnrichrWorkbook <- function(database, full.df, colname) {
    dataset <- full.df[[database]]

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    setColWidths(wb, 1, cols = 2, widths = 45)
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    
    dir.create(file.path("./enrichr", colname), recursive = TRUE)
    filename = str_c(file.path("./enrichr", colname, database), ".xlsx")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

EnrichrPlot <- function(enrichr.df, enrichr.expr, filename, plot.height = 5, plot.width = 8) {
    enrichr.df$Gene.Count <- map(enrichr.df$Genes, str_split, ",") %>% map(unlist) %>% map_int(length)
    enrichr.df$Adj.P.value <- p.adjust(enrichr.df$P.value, method = "fdr")
    enrichr.df$Log.P.value <- log10(enrichr.df$Adj.P.value)
    enrichr.updown <- map(enrichr.df$Genes, UpDown, enrichr.expr) %>% reduce(rbind)
    colnames(enrichr.updown) <- c("Down", "Up")
    enrichr.df <- cbind(enrichr.df, enrichr.updown)
    enrichr.df$Log.Up <- -enrichr.df$Log.P.value * enrichr.df$Up / enrichr.df$Gene.Count
    enrichr.df$Log.Down <- -enrichr.df$Log.P.value * enrichr.df$Down / enrichr.df$Gene.Count
    enrichr.df$Term %<>% str_replace_all("\\ \\(.*$", "") %>% str_replace_all("\\_Homo.*$", "")  #Remove any thing after the left parenthesis and convert to all lower case
    enrichr.df$Format.Name <- str_c(enrichr.df$Term, " (", enrichr.df$Gene.Count, ")")
    enrichr.df %<>% arrange(-Log.P.value)
    enrichr.df$Format.Name %<>% factor(levels = enrichr.df$Format.Name)
    enrichr.df.plot <- select(enrichr.df, Format.Name, Log.Up, Log.Down) %>% gather(Direction, Length, -Format.Name) 

    p <- ggplot(enrichr.df.plot, aes(Format.Name, Length, fill = Direction)) 
    p <- p + geom_bar(stat = "identity", size = 1) 
    p <- p + scale_fill_discrete(name = "Direction", labels = c("Up", "Down")) 
    p <- p + coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab(expression(paste(Log[10], ' P-Value')))
    p <- p + theme(plot.background = element_blank(), legend.background = element_blank(), axis.text.y = element_text(size = 12), panel.border = element_blank())
    p <- p + theme(panel.background = element_blank(), axis.line.x = element_line()) + geom_hline(color = "red", yintercept = -log10(0.05)) 
    CairoPDF(filename, height = plot.height, width = plot.width, bg = "transparent")
    print(p)
    dev.off()
}

UpDown <- function(filter.vector, enrichr.df) {
    split.vector <- str_split(filter.vector, ",")[[1]]
    enrichr.filter <- filter(enrichr.df, is.element(Symbol, split.vector))
    enrichr.vector <- c("Up" = length(which(sign(enrichr.filter$logFC) == 1)), "Down" = length(which(sign(enrichr.filter$logFC) == -1)))
    enrichr.vector
}


#Code for getting size of objects in memory
objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

source("../../FRDA project/common_functions.R") #Load shared functions
targets.final <- read.xlsx("../Jones_Genomic Patient Sample Summary_2016.xlsx") %>% filter(!is.na(HT12.Array.Experiment))#Load phenotype data
colnames(targets.final)[7] <- "Age"
targets.final$Slide.ID <- str_c(targets.final$HT12.General.Array, targets.final$genexstripe.controling.stripe, sep = "_")
targets.reduce <- select(targets.final, Slide.ID, Tissue, Cell.type, Sex, Age) %>% arrange(Slide.ID)
targets.reduce$Age <- str_split_fixed(targets.reduce$Age, "-", 2) %>% apply(2, as.integer) %>% apply(1, mean)
targets.reduce$Combined <- str_c(targets.reduce$Tissue, targets.reduce$Cell.type, sep = "_")

#Why can't lumi read in .csv files correctly?
temp <- read_csv("../UNGC HT12 arrays/2016-9206 sample probe profile.csv")
write.table(temp, "../UNGC HT12 arrays/2016-9206.tsv", sep = "\t", row.names = FALSE) #Save table to disk
lumi.raw <- lumiR("../UNGC HT12 arrays/2016-9206.tsv", lib.mapping = "lumiHumanIDMapping", checkDupId = TRUE, convertNuID = TRUE, QC = FALSE) #Read in the joined batches
sampleNames(lumi.raw) <- targets.reduce$Slide.ID

rownames(targets.reduce) <- targets.reduce$Slide.ID #Set the row names to be the same as the sampleID (may be unnecessary)
pData(lumi.raw) <- targets.reduce #Add phenotype data
SaveRDSgz(lumi.raw, "./save/lumi.raw.rda")

age.pediatric <- lumi.raw$Age < 5
sex.unknown <- lumi.raw$Sex == "Unknown"
combined.key <- !age.pediatric & !sex.unknown
lumi.known <- lumi.raw[,combined.key]

lumi.vst <- lumiT(lumi.known) #Perform variance stabilized transformation
SaveRDSgz(lumi.vst, file = "./save/lumi.vst.rda")

lumi.norm <- lumiN(lumi.vst, method = "rsn") #Normalize with robust spline regression
lumi.cutoff <- detectionCall(lumi.norm) #Get the count of probes which passed the detection threshold per sample
lumi.expr <- lumi.norm[which(lumi.cutoff > 0),] #Drop any probe where none of the samples passed detection threshold
symbols.lumi <- getSYMBOL(rownames(lumi.expr), 'lumiHumanAll.db') %>% is.na #Determine which remaining probes are unannotated
lumi.expr.annot <- lumi.expr[!symbols.lumi,] #Drop any probe which is not annotated
SaveRDSgz(lumi.expr.annot, file = "./save/lumi.expr.annot.rda")

#Regenerate plots
#gen.boxplot("baseline_intensity_norm.jpg", lumi.expr.annot, batch.colors, "RSN normalized signal intensity", "Intensity") #Make box plot of normalized intensities
Histogram("histogram_norm", lumi.expr.annot)
mds.norm <- exprs(lumi.expr.annot) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_status_norm", mds.norm, pData(lumi.expr.annot), "none", "Combined") #label PCs by status

#Collapse the data by symbol
collapse.expr <- collapseRows(exprs(lumi.expr.annot), getSYMBOL(featureNames(lumi.expr.annot), 'lumiHumanAll.db'), rownames(lumi.expr.annot)) #collapseRows by symbol
lumi.collapse <- ExpressionSet(assayData = collapse.expr$datETcollapsed, phenoData = phenoData(lumi.expr.annot))
SaveRDSgz(lumi.collapse, file = "./save/lumi.collapse.rda")

#Remove effects of covariates
model.design <- model.matrix(~ Combined + Sex + Age, data = pData(lumi.collapse)) #Create model matrix of covariates to be removed

rmcov.expr <- removeBatchEffect(exprs(lumi.collapse), covariates = model.design[,3:5], design = model.design[,1:2]) #Remove the effects of covariates, with the difference in diagnoses being supplied as the design argument to preserve those group differences
rmcov.lumi <- lumi.collapse #Make a copy of lumi object
exprs(rmcov.lumi) <- rmcov.expr #Transfer cleaned expression values into new lumi object
SaveRDSgz(rmcov.lumi, file = "./save/rmcov.lumi.rda")
SaveRDSgz(pData(rmcov.lumi), file = "./save/rmcov.pheno.rda")

mds.rmcov <- exprs(rmcov.lumi) %>% t %>% dist(method = "manhattan") %>% cmdscale(eig = TRUE) #get first two principle components
MDSPlot("mds_status_rmcov", mds.rmcov, pData(rmcov.lumi), "none", "Combined") #label PCs by status

#Limma fit
fitb <- lmFit(lumi.collapse, model.design) %>% eBayes

#Create DE table
#Make top tables for each coefficient
toptable.dc <- topTable(fitb, coef = 2, n = Inf) 
toptable.dc$Symbol <- rownames(toptable.dc)
SaveRDSgz(toptable.dc, "./save/toptable.dc.rda")

#Retrieve annotation information from Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm.table <- getBM(attributes = c('hgnc_symbol', 'description'), filters = 'hgnc_symbol', values = as.character(toptable.dc$Symbol), mart = ensembl)
bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(bm.table) <- c("Symbol", "Definition")

toptable.annot <- left_join(toptable.dc, bm.table) %>% 
    select(Symbol, Definition, logFC, P.Value, adj.P.Val, AveExpr, t, B) %>%
    arrange(P.Value)
toptable.annot[,grepl("logFC|P.Val|AveExpr|^t$|B", colnames(toptable.annot))] %<>% signif(3)

DEWorkbook(toptable.annot, "toptable_annotated.xlsx")

#Siggenes - why doesn't this work right?
siggene.ebam.dc <- limma2ebam(fitb, coef = 2) 
ebam2excel(siggene.ebam.dc, 0.9, "siggene.ebam.dc.csv")

a0.all <- find.a0(rmcov.expr, as.integer(factor(rmcov.lumi$Combined)) - 1, B = 1000)
ebam.all <- ebam(a0.all, which.a0 = 3)
ebam.all.df <- data.frame(Z.score = ebam.all@z, Posterior = ebam.all@posterior)
ebam.all.df$Significant <- ebam.all.df$Posterior > 0.9
ebam.all.df$Symbol <- rownames(ebam.all.df)

#Make volcano plots
toptable.dc$Significant <- toptable.dc$adj.P.Val < 0.01
VolcanoPlot(toptable.dc, "volcano_dc")

FilterEnrichr <- function(enrichr.df, size = 100) {
    enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
    enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
    if (nrow(enrichr.df) > size) {
        enrichr.filter %<>% slice(1:size)
    }
    enrichr.filter
}

#Submit genes to Enrichr
source('../../code/GO/enrichr.R')
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016", "Human_Gene_Atlas") 

filter.df <- filter(toptable.dc, adj.P.Val < 0.01) %>% slice(1:1000)
enrichr.data <- map(enrichr.terms, GetEnrichrData, filter.df, FALSE)
enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
enrichr.data <- enrichr.data[!is.na(enrichr.data)]
names(enrichr.data) <- enrichr.names
map(names(enrichr.data), EnrichrWorkbook, enrichr.data, "duodenum_vs_colon")

gobiol.file <- "./enrichr/duodenum_vs_colon/GO_Biological_Process_2015.xlsx" 
gobiol <- read.xlsx(gobiol.file) 
gobiol.filter <- FilterEnrichr(gobiol)
GetKappaCluster(file_path_sans_ext(gobiol.file), gobiol.filter, filter.df$Symbol)
gobiol.filter$Database <- "GO Biological Process"
gobiol.final <- slice(gobiol.filter, c(1, 3, 4, 5, 8))

gomole.file <- "./enrichr/duodenum_vs_colon/GO_Molecular_Function_2015.xlsx"
gomole <- read.xlsx(gomole.file) 
gomole.filter <- FilterEnrichr(gomole)
GetKappaCluster(file_path_sans_ext(gomole.file), gomole.filter, filter.df$Symbol)
gomole.filter$Database <- "GO Molecular Function"
gomole.final <- slice(gomole.filter, c(3, 4, 8))

reactome.file <- "./enrichr/duodenum_vs_colon/Reactome_2016.xlsx"
reactome <- read.xlsx(reactome.file) 
reactome.filter <- FilterEnrichr(reactome)
GetKappaCluster(file_path_sans_ext(reactome.file), reactome.filter, filter.df$Symbol)
reactome.filter$Database <- "Reactome"
reactome.final <- slice(reactome.filter, c(4,8))

kegg.file <- "./enrichr/duodenum_vs_colon/KEGG_2016.xlsx"
kegg <- read.xlsx(kegg.file) 
kegg.filter <- FilterEnrichr(kegg)
GetKappaCluster(file_path_sans_ext(kegg.file), kegg.filter, filter.df$Symbol)
kegg.filter$Database <- "KEGG"
kegg.final <- slice(kegg.filter, c(4))

enrichr.final <- rbind(gobiol.final, gomole.final, kegg.final, reactome.final)
EnrichrPlot(enrichr.final, filter.df, "enrichr")

top5.symbol <- arrange(toptable.annot, adj.P.Val)$Symbol[1:10]
top5.expr <- t(exprs(lumi.collapse)[top5.symbol,])
colnames(top5.expr) <- top5.symbol
top5.df <- data.frame(Combined = as.character(lumi.collapse$Combined), top5.expr) %>% gather(Gene, Expression, -Combined)
top5.df$Gene %<>% factor(levels = top5.symbol)
top5.df$Combined %<>% str_replace_all("_", " ") %>% factor(levels = c("duodenum crypts", "colon crypts")) 

p <- ggplot(top5.df, aes(x = Combined, y = Expression, color = Combined)) + geom_jitter() + geom_boxplot() + theme_bw()
p <- p + facet_wrap(~ Gene, ncol = 5, scales = "free") + theme(legend.position = "none")
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank())
CairoPDF("top5", height = 8, width = 16, bg = "transparent")
print(p)
dev.off()

GeneBoxplot(lumi.collapse, "GATA4")
GeneBoxplot(lumi.collapse, "GATA5")
GeneBoxplot(lumi.collapse, "MUC17")
GeneBoxplot(lumi.collapse, "CRADD")
GeneBoxplot(lumi.collapse, "C10orf99")
GeneBoxplot(lumi.collapse, "HOXB8")
GeneBoxplot(lumi.collapse, "LGALS1")
GeneBoxplot(lumi.collapse, "ALPI", TRUE)
GeneBoxplot(lumi.collapse, "FABP2", TRUE)
GeneBoxplot(lumi.collapse, "SLC5A1", TRUE)

#Only different in crypts

de.sig.up <- filter(toptable.annot, adj.P.Val < 0.01 & logFC > 0) %>% select(Symbol, logFC, P.Value, adj.P.Val)
colnames(de.sig.up)[-1] %<>% str_c(".expr")
de.sig.down <- filter(toptable.annot, adj.P.Val < 0.01 & logFC < 0) %>% select(Symbol, logFC, P.Value, adj.P.Val)
colnames(de.sig.down)[-1] %<>% str_c(".expr")

all.promoters <- read.xlsx("../baseline_methylation/promoters.combined.xlsx") %>% as_tibble %>% filter(!is.na(Symbol) & nchar(Symbol) > 0) %>% filter(Symbol %in% toptable.annot$Symbol)
all.genes <- read.xlsx("../baseline_methylation/genes.combined.xlsx") %>% as_tibble %>% filter(!is.na(Symbol) & nchar(Symbol) > 0) %>% filter(Symbol %in% toptable.annot$Symbol)

promoters.crypts.up <- read.xlsx("../baseline_methylation/promoters.crypts.up.xlsx") %>% filter(Symbol %in% toptable.annot$Symbol)
promoters.crypts.down <- read.xlsx("../baseline_methylation/promoters.crypts.down.xlsx") %>% filter(Symbol %in% toptable.annot$Symbol)

genes.crypts.up <- read.xlsx("../baseline_methylation/genes.crypts.up.xlsx") %>% filter(Symbol %in% toptable.annot$Symbol)
genes.crypts.down <- read.xlsx("../baseline_methylation/genes.crypts.down.xlsx") %>% filter(Symbol %in% toptable.annot$Symbol)

expr.up.promoters <- filter(de.sig.up, Symbol %in% all.promoters$Symbol)
up.promoters <- inner_join(expr.up.promoters, promoters.crypts.down)
write.xlsx(up.promoters, "up.promoters.xlsx")
expr.down.promoters <- filter(de.sig.down, Symbol %in% all.promoters$Symbol)
down.promoters <- inner_join(de.sig.down, promoters.crypts.up)
write.xlsx(down.promoters, "down.promoters.xlsx")

promoters.down.table <- c(nrow(up.promoters), nrow(expr.up.promoters) - nrow(up.promoters), nrow(promoters.crypts.down) - nrow(up.promoters), nrow(all.promoters) - nrow(expr.up.promoters) - nrow(promoters.crypts.down) + nrow(up.promoters)) 
dim(promoters.down.table) <- c(2,2)
colnames(promoters.down.table) <- c("DE", "Not DE")
rownames(promoters.down.table) <- c("DM", "Not DM")
write.csv(promoters.down.table, "promoters.down.table.csv", quote = FALSE)

promoters.up.table <- c(nrow(down.promoters), nrow(expr.down.promoters) - nrow(down.promoters), nrow(promoters.crypts.up) - nrow(down.promoters), nrow(all.promoters) - nrow(expr.down.promoters) - nrow(promoters.crypts.up) + nrow(down.promoters)) 
dim(promoters.up.table) <- c(2,2)
colnames(promoters.up.table) <- c("DE", "Not DE")
rownames(promoters.up.table) <- c("DM", "Not DM")
write.csv(promoters.up.table, "promoters.up.table.csv", quote = FALSE)

expr.up.genes <- filter(de.sig.up, Symbol %in% all.genes$Symbol)
up.genes <- inner_join(de.sig.up, genes.crypts.down)
write.xlsx(up.genes, "up.genes.xlsx")
expr.down.genes <- filter(de.sig.down, Symbol %in% all.genes$Symbol)
down.genes <- inner_join(de.sig.down, genes.crypts.up)
write.xlsx(down.genes, "down.genes.xlsx")

genes.down.table <- c(nrow(up.genes), nrow(expr.up.genes) - nrow(up.genes), nrow(genes.crypts.down) - nrow(up.genes), nrow(all.genes) - nrow(expr.up.genes) - nrow(genes.crypts.down) + nrow(up.genes)) 
dim(genes.down.table) <- c(2,2)
colnames(genes.down.table) <- c("DE", "Not DE")
rownames(genes.down.table) <- c("DM", "Not DM")
write.csv(genes.down.table, "genes.down.table.csv", quote = FALSE)

genes.up.table <- c(nrow(down.genes), nrow(expr.down.genes) - nrow(down.genes), nrow(genes.crypts.up) - nrow(down.genes), nrow(all.genes) - nrow(expr.down.genes) - nrow(genes.crypts.up) + nrow(down.genes)) 
dim(genes.up.table) <- c(2,2)
colnames(genes.up.table) <- c("DE", "Not DE")
rownames(genes.up.table) <- c("DM", "Not DM")
write.csv(genes.up.table, "genes.up.table.csv", quote = FALSE)

#Also different in mucosa

promoters.dccc.up <- read.xlsx("../baseline_methylation/promoters.dccc.up.xlsx")
promoters.dccc.down <- read.xlsx("../baseline_methylation/promoters.dccc.down.xlsx")

genes.dccc.up <- read.xlsx("../baseline_methylation/genes.dccc.up.xlsx")
genes.dccc.down <- read.xlsx("../baseline_methylation/genes.dccc.down.xlsx")

up.promoters.dccc <- inner_join(de.sig.up, promoters.dccc.down) %>% as_tibble
write.xlsx(up.promoters.dccc, "up.promoters.dccc.xlsx")
down.promoters.dccc <- inner_join(de.sig.down, promoters.dccc.up) %>% as_tibble
write.xlsx(down.promoters.dccc, "down.promoters.dccc.xlsx")

up.genes.dccc <- inner_join(de.sig.up, genes.dccc.down) %>% as_tibble
write.xlsx(up.genes.dccc, "up.genes.dccc.xlsx")
down.genes.dccc <- inner_join(de.sig.down, genes.dccc.up) %>% as_tibble
write.xlsx(down.genes.dccc, "down.genes.dccc.xlsx")

clock.transcripts <- read.xlsx("../baseline_methylation/clock.transcripts.xlsx") %>% as_tibble

genes.bicor <- bicor(t(exprs(lumi.collapse)), lumi.collapse$Age) %>% 
    signif(3) %>%
    as_tibble %>% 
    set_colnames("Correlation") %>%
    mutate(Symbol = rownames(lumi.collapse)) 
genes.bicor$P.Value <- corPvalueStudent(genes.bicor$Correlation, nSamples = nrow(lumi.collapse)) %>% p.adjust("fdr") %>% signif(3)
genes.bicor %<>% select(Symbol, Correlation, P.Value) %>% arrange(P.Value)
write.xlsx(genes.bicor, "genes.age.cor.xlsx")

candidate.list <- c("GATA4","GATA5","MUC17","CRADD","C10orf99","HOXB8","LGALS1")

genes.candidates <- exprs(lumi.collapse[candidate.list,]) %>% 
    t %>% 
    as_tibble %>%
    mutate(Age = lumi.collapse$Age) %>%
    mutate(Cell.Type = lumi.collapse$Combined) %>%
    gather(Gene, Expression, -Age, -Cell.Type)

p <- ggplot(genes.candidates, aes(x = Age, y = Expression, color = Cell.Type)) + stat_smooth(color = "blue", method = "lm") + geom_point() + theme_bw() 
p <- p + theme(legend.background = element_blank()) + xlab("Age")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_blank())
p <- p + facet_wrap(~ Gene, ncol = 4, scales = "free")
CairoPDF("age_plot", height = 7, width = 16, bg = "transparent")
print(p)
dev.off()
