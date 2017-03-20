library(RnBeads)
library(sva)
library(biomaRt)
library(Metrics)
library(PMCMR)
#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(heatmap.plus)
library(gplots) #for heatmap.2
library(RColorBrewer)
library(UpSetR)

#Reading and writing tables
library(readr)
library(openxlsx)

#library(matrixStats)
library(R.utils)

#Data arrangement
library(reshape2)
library(dplyr)
library(parallel)
library(data.table)

#Functional programming
library(magrittr)
library(purrr)
library(tidyr)
library(broom)
library(WGCNA)

#THIS MUST BE FIXED IN PACKAGE
rnb.execute.na.removal.fix <- function (rnb.set, threshold = rnb.getOption("filtering.missing.value.quantile")) {
    if (!inherits(rnb.set, "RnBSet")) {
        stop("invalid value for rnb.set")
    }
    if (!(is.double(threshold) && length(threshold) == 1 && (!is.na(threshold)))) {
        stop("invalid value for threshold")
    }
    if (!(0 <= threshold && threshold <= 1)) {
        stop("invalid value for threshold; expected a value between 0 and 1")
    }
    filterRes <- RnBeads:::rnb.execute.na.removal.internal(rnb.set, NULL,
        threshold)
    list(dataset.before = rnb.set, dataset = remove.sites(rnb.set,
        filterRes$filtered), filtered = filterRes$filtered, threshold = threshold,
        naCounts = filterRes$naCounts)
}

GetSizes <- function(dataset, pval.column, log.column, p.val) {
    dataset.sig <- filter_(dataset, str_c(pval.column, " < ", p.val))
    dataset.up <- filter_(dataset.sig, str_c(log.column, " > 0"))
    dataset.down <- filter_(dataset.sig, str_c(log.column, " < 0"))
    return(list(positive = dim(dataset.up)[1], negative = -(dim(dataset.down)[1])))
}

DecidePlot <- function(file.name, decide.plot, y.lab, bar.padding = 300, pos.hjust = -0.4, neg.hjust = 1.3) {
    y.max <- max(decide.plot$Num.Sites) + nchar(max(decide.plot$Num.Sites)) * bar.padding
    y.min = min(decide.plot$Num.Sites) - nchar(abs(min(decide.plot$Num.Sites))) * bar.padding

    p <- ggplot()
    p <- p + geom_bar(data = filter(decide.plot, Direction == "positive"),  aes(x = Comparison, y = Num.Sites), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = filter(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Comparison, y = Num.Sites, hjust = pos.hjust, label = Num.Sites), position = position_dodge(width = 1))
    p <- p + geom_bar(data = filter(decide.plot, Direction == "negative"),  aes(x = Comparison, y = Num.Sites), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = filter(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Comparison, y = Num.Sites, hjust = neg.hjust, label = abs(Num.Sites)), position = position_dodge(width = 1))
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab(y.lab)
    p <- p + theme(plot.background = element_blank()) + ylim(y.min, y.max) + facet_grid(Test + Total.Sites ~ .) 
    CairoPDF(file.name, width = 7, height = 8, bg = "transparent")
    print(p)
    dev.off()
}

Heatmap <- function(meths, ngenes, file.name, cluster.object) {
    sites.ordered <-  apply(meths, 1, mad) %>% order(decreasing = TRUE)
    sites.plot <- meths[sites.ordered[1:ngenes],]
    cluster.tree <- cluster.object[[7]]@result
    attr(cluster.tree, "class") <- "hclust"
    cluster.dendro <- as.dendrogram(cluster.tree)

    CairoPDF(file.name, width = 10, height = 10)
    heatmap.2(sites.plot, Rowv = TRUE, Colv = cluster.dendro, dendrogram = "both", scale = "none", trace = "none", labRow = FALSE)
    dev.off()
}

UniqueNames <- function(data.list) {
    comparison.name <- names(data.list)
    df.extract <- data.list[[1]]
    colnames(df.extract)[5] <- str_c(comparison.name, ".log2FC")
    colnames(df.extract)[6] <- str_c(comparison.name, ".p.val")
    colnames(df.extract)[7] <- str_c(comparison.name, ".p.val.adj")
    colnames(df.extract)[8] <- str_c(comparison.name, ".combinedRank")
    list(df.extract)
}

DMWorkbook <- function(dataset, filename) {
    pval.cols <- colnames(dataset) %>% str_detect("p.val$") %>% which
    adj.pval.cols <- colnames(dataset) %>% str_detect("p.val.adj") %>% which
    logfc.cols <- colnames(dataset) %>% str_detect("log2FC") %>% which
    description.cols <- colnames(dataset) %>% str_detect("Description") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.005", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = adj.pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = logfc.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = 20)
    setColWidths(wb, 1, cols = 2, widths = "auto")
    setColWidths(wb, 1, cols = 4:ncol(dataset), widths = "auto")
    setColWidths(wb, 1, cols = description.cols, widths = 45)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

GetEnrichr <- function(comparison, submit.df, cutoff, enrichr.terms, site.type) {
    comparison.pval <- str_c(comparison, "p.val.adj", sep = ".")
    filter.df <- filter_(submit.df, str_c(comparison.pval, "<", cutoff)) %>% slice(1:1000)
    print(dim(filter.df))
    enrichr.data <- map(enrichr.terms, GetEnrichrData, filter.df, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), EnrichrWorkbook, enrichr.data, comparison, site.type)
}

FilterEnrichr <- function(enrichr.df, size = 100) {
    enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
    enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
    if (nrow(enrichr.df) > size) {
        enrichr.filter %<>% slice(1:size)
    }
    enrichr.filter
}

EnrichrWorkbook <- function(subindex, full.df, comparison, site.type) {
    dataset <- full.df[[subindex]]
    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    setColWidths(wb, 1, cols = c(1, 3:ncol(dataset)), widths = "auto")
    setColWidths(wb, 1, cols = 2, widths = 45)

    sub.dir <- file.path("./enrichr", site.type, comparison)
    dir.create(sub.dir, showWarnings = FALSE, recursive = TRUE)
    filename = str_c(sub.dir, "/", subindex, ".xlsx", sep = "")
    saveWorkbook(wb, filename, overwrite = TRUE) 
}

VolcanoPlot <- function(top.table, filename, maintitle, cutoff = 0.05, cutoff.column = "P.Value", log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Log.Pvalue") {
    top.table$Significant <- factor(top.table[[cutoff.column]] < cutoff)
    top.table$Log.Pvalue <- -log10(top.table[[cutoff.column]])
    p <- ggplot(top.table, aes_string(x = log.column, y = "Log.Pvalue")) + geom_point(aes(color = Significant))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(maintitle)
    p <- p + theme(legend.position = "none", plot.background = element_blank())
    p <- p + theme(panel.border = element_rect(color = "black", size = 1))
    p <- p + xlab(xlabel) + ylab(ylabel)
    CairoPDF(filename, width = 6, height = 6, bg = "transparent")
    print(p)
    dev.off()
}

MapVolcanoPlots <- function(comparison, prefix, tab.list) {
    tab.reduce <- select(as.data.frame(tab.list[[comparison]]), mean.mean.quot.log2, comb.p.adj.fdr)
    comparison.title <- str_replace_all(comparison, "_", " ") %>% str_replace_all("m\\.", "m ") %>% str_replace_all("n\\.", "n ")
    colnames(tab.reduce) <- c("logFC", "P.Value")
    VolcanoPlot(tab.reduce, str_c(prefix, comparison, sep = "_"), comparison.title)
}

Top5Plot <- function(contrast, diffmeth.table, meth.matrix, pheno.df, prefix) {
    rank.column <- str_c(contrast, "combinedRank", sep = ".")
    top5.ensembl <- arrange_(diffmeth.table, rank.column)$Ensembl.ID[1:5]
    top5.symbol <- arrange_(diffmeth.table, rank.column)$Symbol[1:5]
    top5.beta <- t(meth.matrix[top5.ensembl,])
    colnames(top5.beta) <- top5.symbol
    top5.df <- data.frame(Combined = as.character(pheno.df$Combined.Factor), top5.beta) %>% gather(Gene, Methylation, -Combined)
    top5.df$Gene %<>% factor(levels = top5.symbol)
    top5.df$Combined %<>% str_replace_all("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts")) 

    contrast.title <- str_replace_all(contrast, "_", " ") %>% str_replace_all("m\\.", "m ") %>% str_replace_all("n\\.", "n ")
    p <- ggplot(top5.df, aes(x = Combined, y = Methylation, color = Combined)) + geom_jitter() + geom_boxplot() + theme_bw()
    p <- p + facet_wrap(~ Gene, ncol = 5, scales = "free") + theme(legend.position = "none")
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
    p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + ggtitle(contrast.title) + ylim(0,1)
    CairoPDF(str_c("top5", prefix, contrast, sep = "_"), height = 4, width = 16, bg = "transparent")
    print(p)
    dev.off()
}

EnrichrPlot <- function(enrichr.df, enrichr.expr, filename, plot.title, plot.height = 5, plot.width = 8) {
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

objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

source("../../FRDA project/common_functions.R")

sample.data <- read.xlsx("../Joneslab_Methylation_Annotated_new.xlsx")
colnames(sample.data)[5:6] <- c("barcode", "Slide")
sample.data$barcode %<>% str_replace("9016", "90016") %>% str_replace("9006", "90006")
sample.age <- sample.data$Chronological.Age %>% as.character
sample.age[str_detect(sample.age, "months")] <- "0-5"
sample.age[str_detect(sample.age, "fetus")] <- "0-5"
sample.age[str_detect(sample.age, "years")] <- "70-75"
sample.age %<>% factor
sample.data$Chronological.Age <- sample.age 
sample.data$Cell.Type %<>% str_replace(" ", "") %>% tolower
sample.data$Batch <- factor(sample.data$Slide) %>% as.integer
sample.reduce <- select(sample.data, Tissue, Cell.Type, barcode, Chronological.Age, Sex)
sample.reduce$Tissue %<>% str_replace_all("small intestine", "duodenum")
sample.reduce$Batch <- 1
sample.reduce[grepl("2005", sample.reduce$barcode),]$Batch <- 2

new.pheno <- read.xlsx("./pheno_gender.xlsx") %>% select(Tissue, Cells, barcode, Age, Predicted.Gender)
colnames(new.pheno) <- c("Tissue", "Cell.Type", "barcode", "Chronological.Age", "Sex")
new.pheno$Tissue %<>% str_replace_all(" ", "")
new.pheno$Cell.Type %<>% str_replace_all("spheroids", "spheroid") 
new.pheno$Batch <- 3
new.pheno$Chronological.Age %<>% str_replace_all("0-1", "0-5") %>% str_replace_all("17", "15-20")
sample.data.final <- rbind(sample.reduce, new.pheno)
write_csv(sample.data.final, "../sample_annotation.csv")

#new.targets <- read.csv("../new_methylation/2016-9177 Sample Sheet.csv")
#new.targets$External.Sample.ID %<>% str_replace_all("-", "/") %>% str_replace_all("_", " ")
#pheno.join <- left_join(new.pheno, new.targets)
#pheno.join$barcode <- str_c(pheno.join$chip.ID, pheno.join$stripe, sep = "_")
#write_csv(pheno.join, "../new_sample_annotations.csv")
#write.xlsx(pheno(rnb.set), "pheno_gender.xlsx")

idat.dir <- "../Raw_Data"
sample.annotation <- "../sample_annotation.csv"
report.dir <- "./reports"

rnb.initialize.reports(report.dir)
logger.start(fname = NA)
parallel.setup(7)
options(fftempdir = "~/tmp/Rtmp")
#set.tempdir("~/tmp/Rtmp")
rnb.options(logging.disk = TRUE)

data.source <- c(idat.dir, sample.annotation)
rnb.options(import.gender.prediction = FALSE)
result <- rnb.run.import(data.source = data.source, data.type = "infinium.idat.dir", dir.reports = report.dir)
rnb.set <- result$rnb.set
SaveRDSgz(rnb.set, "./save/rnb.set.rda")
remove.tissue <- grepl("unknown", pheno(rnb.set)$Tissue)
remove.celltype <- grepl("wholebowel", pheno(rnb.set)$Cell.Type)
remove.age <- grepl("^0-5$", pheno(rnb.set)$Chronological.Age)

#Code for fixing various issues with age
remove.all <- remove.tissue | remove.celltype | remove.age
rnb.known <- remove.samples(rnb.set, remove.all)
SaveRDSgz(rnb.known, "./save/rnb.known.rda")

pheno.known <- pheno(rnb.known)
pheno.known$Cell.Type %<>% droplevels
lm.celltype.age <- lm(as.integer(Cell.Type) ~ Chronological.Age, pheno.known) %>% anova %>% tidy
lm.celltype.sex <- lm(as.integer(Cell.Type) ~ Sex, pheno.known) %>% anova %>% tidy
lm.tissue.age <- lm(as.integer(Tissue) ~ Chronological.Age, pheno.known) %>% anova %>% tidy
lm.tissue.sex <- lm(as.integer(Tissue) ~ Sex, pheno.known) %>% anova %>% tidy

rnb.run.qc(rnb.known, report.dir)
raw.beta <- meth(rnb.known)
rownames(raw.beta) <- rownames(annotation(rnb.known, type = "sites"))
cpgset <- read_csv("./datMiniAnnotation.csv")
cpg.match <- match(cpgset$Name, rownames(raw.beta))
reduce.beta <- raw.beta[cpgset$Name,]
write.csv(raw.beta, "./allsites_beta.csv")

pheno.export <- mutate(pheno(rnb.known), ID = 1:ncol(raw.beta))
colnames(pheno.export)[4:5] <- c("Age", "Female")
pheno.export$Age <- str_split_fixed(pheno.export$Age, "-", 2) %>% apply(2, as.integer) %>% apply(1, mean)
pheno.export$Female %<>% equals("female") %>% as.integer
write.csv(pheno.export, "./pheno_export.csv", row.names = FALSE)

rnb.filter <- rnb.execute.context.removal(rnb.known)$dataset

rnb.filter <- rnb.execute.snp.removal(rnb.filter, snp = "any")$dataset

rnb.filter <- rnb.execute.sex.removal(rnb.filter)$dataset

rnb.greedy <- rnb.execute.greedycut(rnb.filter)
filter.sites <- rnb.greedy$sites

rnb.filter <- remove.sites(rnb.filter, filter.sites)


rnb.filter <- rnb.execute.na.removal.fix(rnb.filter, 0)$dataset

rnb.filter <- rnb.execute.variability.removal(rnb.filter, 0.005)$dataset

rnb.norm <- rnb.execute.normalization(rnb.filter, method = "bmiq", bgcorr.method = "methylumi.noob")
SaveRDSgz(rnb.norm, "./save/rnb.norm.rda")

dred.sites <- rnb.execute.dreduction(rnb.norm)
dred.promoters <- rnb.execute.dreduction(rnb.norm, target = "promoters")
dred <- list(sites = dred.sites, promoters = dred.promoters)

dred.plot <- data.frame(dred.promoters$mds$manhattan[,1:2])
colnames(dred.plot) <- c("PC1", "PC2")
dred.plot %<>% mutate(Cell.Type = pheno(rnb.norm)$Cell.Type)
dred.plot %<>% mutate(Age = pheno(rnb.norm)$Chronological.Age)

#Convert this plot ggplot
p <- ggplot(data = dred.plot, aes(x = PC1, y = PC2, col = Cell.Type)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
CairoPDF("promoters_celltype", width = 7, height = 6)
print(p)
dev.off()

p <- ggplot(data = dred.plot, aes(x = PC1, y = PC2, col = Age)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
CairoPDF("promoters_age", width = 7, height = 6)
print(p)
dev.off()

dred.allsites <- data.frame(dred.sites$mds$manhattan[,1:2])
colnames(dred.allsites) <- c("PCA1", "PCA2")
dred.allsites %<>% mutate(Cell.Type = pheno(rnb.norm)$Cell.Type)
dred.allsites %<>% mutate(Age = pheno(rnb.norm)$Chronological.Age)

CairoPDF("allsites_celltype", width = 7, height = 6)
p <- ggplot(data = dred.allsites, aes(x = PCA1, y = PCA2, col = Cell.Type)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()

CairoPDF("allsites_age", width = 7, height = 6)
p <- ggplot(data = dred.allsites, aes(x = PCA1, y = PCA2, col = Age)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p)
dev.off()

promoters.beta <- meth(rnb.norm, type = "promoters")
rownames(promoters.beta) <- rownames(annotation(rnb.norm, type = "promoters"))
SaveRDSgz(promoters.beta, "./save/promoters.beta.rda")
genes.beta <- meth(rnb.norm, type = "genes")
rownames(genes.beta) <- rownames(annotation(rnb.norm, type = "genes"))
sites.beta <- meth(rnb.norm)
pheno.export <- pheno(rnb.norm)
SaveRDSgz(pheno.export, "./save/pheno.export.rda")

model.design <- model.matrix( ~ droplevels(Tissue) + droplevels(Cell.Type) + droplevels(Chronological.Age) + Sex, pheno(rnb.norm))
promoters.combat <- ComBat(dat = promoters.beta, batch = pheno(rnb.norm)$Batch, mod = model.design)
rownames(promoters.combat) <- rownames(annotation(rnb.norm, type = "promoters"))
SaveRDSgz(promoters.combat, "./save/promoters.combat.rda")

genes.combat <- ComBat(dat = genes.beta, batch = pheno(rnb.norm)$Batch, mod = model.design)
site.combat <- ComBat(dat = sites.beta, batch = pheno(rnb.norm)$Batch, mod = model.design)

model.cov <- model.matrix( ~ droplevels(Chronological.Age) + Sex, pheno.export)
SaveRDSgz(model.design, "./save/model.design.rda")
SaveRDSgz(model.cov, "./save/model.cov.rda")

promoters.rmcov <- removeBatchEffect(promoters.combat, covariates = model.cov[,-1], design = model.design[,1:5])
SaveRDSgz(promoters.rmcov, "./save/promoters.rmcov.rda")

sites.rmcov <- removeBatchEffect(site.combat, covariates = model.cov[,-1], design = model.design[,1:5])
SaveRDSgz(sites.rmcov, "./save/sites.rmcov.rda")

genes.rmcov <- removeBatchEffect(genes.combat, covariates = model.cov[,-1], design = model.design[,1:5])
SaveRDSgz(genes.rmcov, "./save/genes.rmcov.rda")

promoters.mds <- isoMDS(dist(t(promoters.rmcov), method = "manhattan"))
dred.promoters.rmcov <- data.frame(promoters.mds$points)
colnames(dred.promoters.rmcov) <- c("PC1", "PC2")
dred.promoters.rmcov %<>% mutate(Cell.Type = pheno.export$Cell.Type)
dred.promoters.rmcov %<>% mutate(Age = pheno.export$Chronological.Age)
dred.promoters.rmcov %<>% mutate(Sex = pheno.export$Sex)
dred.promoters.rmcov %<>% mutate(Tissue = pheno.export$Tissue)
dred.promoters.rmcov %<>% mutate(Tissue.CellType = factor(str_c(pheno.export$Tissue, pheno.export$Cell.Type, sep = " ")))

sites.mds <- isoMDS(dist(t(sites.rmcov), method = "manhattan"))
dred.sites.rmcov <- data.frame(sites.mds$points)
colnames(dred.sites.rmcov) <- c("PC1", "PC2")
dred.sites.rmcov %<>% mutate(Cell.Type = pheno.export$Cell.Type)
dred.sites.rmcov %<>% mutate(Age = pheno.export$Chronological.Age)
dred.sites.rmcov %<>% mutate(Sex = pheno.export$Sex)
dred.sites.rmcov %<>% mutate(Tissue = pheno.export$Tissue)
dred.sites.rmcov %<>% mutate(Tissue.CellType = factor(str_c(pheno.export$Tissue, pheno.export$Cell.Type, sep = " ")))

genes.mds <- isoMDS(dist(t(genes.rmcov), method = "manhattan"))
dred.genes.rmcov <- data.frame(genes.mds$points)
colnames(dred.genes.rmcov) <- c("PC1", "PC2")
dred.genes.rmcov %<>% mutate(Cell.Type = pheno.export$Cell.Type)
dred.genes.rmcov %<>% mutate(Age = pheno.export$Chronological.Age)
dred.genes.rmcov %<>% mutate(Sex = pheno.export$Sex)
dred.genes.rmcov %<>% mutate(Tissue = pheno.export$Tissue)
dred.genes.rmcov %<>% mutate(Tissue.CellType = factor(str_c(pheno.export$Tissue, pheno.export$Cell.Type, sep = " ")))

#dred.sites.dist <- dist(sites.mds$points, method = "euclidean", diag = TRUE, upper = TRUE)
#dred.sites.cor <- bicor(c(dist(t(sites.rmcov), "euclidean", diag = TRUE, upper = TRUE)), c(dred.sites.dist))
#dred.sites.rsquared <- dred.sites.cor ^ 2
#dred.freedom <- length(c(dred.sites.dist)) - 2
#dred.fvalue <- dred.sites.rsquared / ((1 - dred.sites.rsquared) / dred.freedom)
#dred.pvalue <- pf(dred.fvalue, 1, dred.freedom, lower.tail = FALSE)

#dred.promoters.dist <- dist(promoters.mds$points, method = "euclidean", diag = TRUE, upper = TRUE)
#dred.promoters.cor <- bicor(c(dist(t(promoters.rmcov), "euclidean", diag = TRUE, upper = TRUE)), c(dred.promoters.dist))
#dred.promoters.rsquared <- dred.promoters.cor ^ 2
#dred.promoters.freedom <- length(c(dred.promoters.dist)) - 2
#dred.promoters.fvalue <- dred.promoters.rsquared / ((1 - dred.promoters.rsquared) / dred.promoters.freedom)
#dred.promoters.pvalue <- pf(dred.promoters.fvalue, 1, dred.promoters.freedom, lower.tail = FALSE)

PCAPlot <- function(plot.df, color.column, scale.name, filename){
    p <- ggplot(data = plot.df, aes_string(x = "PC1", y = "PC2", col = color.column)) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank(), legend.background = element_blank())
    p <- p + theme(panel.border = element_rect(color = "black", size = 1))
    p <- p + scale_color_discrete(name = scale.name, labels = capitalize(levels(plot.df[[color.column]])))
    CairoPDF(filename, width = 8, height = 6, bg = "transparent")
    print(p)
    dev.off()
}

PCAPlot(dred.promoters.rmcov, "Cell.Type", "Cell Type", "promoters_celltype_rmcov")
PCAPlot(dred.promoters.rmcov, "Age", "Age", "promoters_age_rmcov")
PCAPlot(dred.promoters.rmcov, "Sex", "Sex", "promoters_sex_rmcov")
PCAPlot(dred.promoters.rmcov, "Tissue.CellType", "Tissue + Cell Type", "promoters_combined_rmcov")

PCAPlot(dred.genes.rmcov, "Cell.Type", "Cell Type", "genes_celltype_rmcov")
PCAPlot(dred.genes.rmcov, "Age", "Age", "genes_age_rmcov")
PCAPlot(dred.genes.rmcov, "Sex", "Sex", "genes_sex_rmcov")
PCAPlot(dred.genes.rmcov, "Tissue.CellType", "Tissue + Cell Type", "genes_combined_rmcov")

PCAPlot(dred.sites.rmcov, "Cell.Type", "Cell Type", "sites_celltype_rmcov")
PCAPlot(dred.sites.rmcov, "Age", "Age", "sites_age_rmcov")
PCAPlot(dred.sites.rmcov, "Sex", "Sex", "sites_sex_rmcov")
PCAPlot(dred.sites.rmcov, "Tissue.CellType", "Tissue + Cell Type", "sites_combined_rmcov")

age.fix <- pheno(rnb.norm)$Chronological.Age %>% droplevels %>% str_replace_all(" ", "")
rnb.norm <- addPheno(rnb.norm, age.fix, "Age.fixed")
main.comparisons <- !grepl("crypt|mucosa", pheno(rnb.norm)$Combined)
rnb.compare <- remove.samples(rnb.norm, main.comparisons)
compare.design <- model.matrix( ~ droplevels(Tissue) + droplevels(Cell.Type) + droplevels(Chronological.Age) + Sex, pheno(rnb.compare))
SaveRDSgz(rnb.compare, "./save/rnb.compare.rda")

promoters.compare <- meth(rnb.compare, type = "promoters")
rownames(promoters.compare) <- rownames(annotation(rnb.compare, type = "promoters"))
promoters.compare.combat <- ComBat(dat = promoters.compare, batch = pheno(rnb.compare)$Batch, mod = compare.design)

genes.compare <- meth(rnb.compare, type = "genes")
rownames(genes.compare) <- rownames(annotation(rnb.compare, type = "genes"))
genes.compare.combat <- ComBat(dat = genes.compare, batch = pheno(rnb.compare)$Batch, mod = compare.design)

promoters.compare.rmcov <- removeBatchEffect(promoters.compare.combat, covariates = compare.design[,4:ncol(compare.design)], design = compare.design[,1:3])
genes.compare.rmcov <- removeBatchEffect(genes.compare.combat, covariates = compare.design[,4:ncol(compare.design)], design = compare.design[,1:3])
SaveRDSgz(promoters.compare.rmcov, "./save/promoters.compare.rda")
SaveRDSgz(genes.compare.rmcov, "./save/genes.compare.rda")

rnb.options("covariate.adjustment.columns" = c("Age.fixed", "Sex", "Batch"))
combined <- str_c(pheno(rnb.compare)$Tissue, pheno(rnb.compare)$Cell.Type, sep = ".") %>% factor(levels = c("duodenum.mucosa", "duodenum.crypts", "colon.mucosa", "colon.crypts")) %>% droplevels
rnb.compare <- addPheno(rnb.compare, combined, "Combined.Factor")
SaveRDSgz(rnb.compare, "./save/rnb.compare.rda")
SaveRDSgz(pheno(rnb.compare), "./save/pheno.compare.rda")
comp.cols <- "Combined.Factor"
reg.types <- c("genes", "promoters")

rnb.options(disk.dump.big.matrices = FALSE)
rnb.options(enforce.memory.management = TRUE)
diffmeth.adj <- rnb.execute.computeDiffMeth(rnb.compare, pheno.cols = comp.cols, pheno.cols.all.pairwise = comp.cols, region.types = reg.types)
SaveRDSgz(diffmeth.adj, "./save/diffmeth.adj.rda")

comparisons <- get.comparisons(diffmeth.adj)
comparisons.keep <- comparisons[c(1:2,5:6)]
comparisons.format <- str_replace(comparisons.keep, " \\(based on Combined.Factor\\)", "") %>% str_replace_all(" ", "_")

tab.sites <- map(comparisons.keep, get.table, object = diffmeth.adj, region.type = "sites", return.data.frame = TRUE) 
names(tab.sites) <- comparisons.format
tab.promoters <- map(comparisons.keep, get.table, object = diffmeth.adj, region.type = "promoters", return.data.frame = TRUE) 
names(tab.promoters) <- comparisons.format
tab.genes <- map(comparisons.keep, get.table, object = diffmeth.adj, region.type = "genes", return.data.frame = TRUE) 
names(tab.genes) <- comparisons.format

SaveRDSgz(tab.sites, "./save/tab.sites.rda")

map(comparisons.format, MapVolcanoPlots, "promoters", tab.promoters)
map(comparisons.format, MapVolcanoPlots, "genes", tab.genes)
map(comparisons.format, MapVolcanoPlots, "sites", tab.sites)

#Get Thresholds
thresholds <- c(0.01, 0.005, 0.001)
thresholds.fdr <- c(0.05, 0.01, 0.005)

promoters.threshold <- map(thresholds, . %>% map(.x = tab.promoters, .f = GetSizes, pval.column = "comb.p.val", log.column = "mean.mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons.keep, length(thresholds)))
promoters.threshold$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(promoters.threshold)[3] <- "Test"
promoters.threshold$Test <- str_c("p < ", promoters.threshold$Test)

promoters.threshold.fdr <- map(thresholds.fdr, . %>% map(.x = tab.promoters, .f = GetSizes, pval.column = "comb.p.adj.fdr", log.column = "mean.mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds.fdr, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons.keep, length(thresholds.fdr)))
promoters.threshold.fdr$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(promoters.threshold.fdr)[3] <- "Test"
promoters.threshold.fdr$Test <- str_c("FDR p < ", promoters.threshold.fdr$Test) 

promoters.threshold.combine <- rbind(promoters.threshold, promoters.threshold.fdr) 
promoters.threshold.combine$positive %<>% unlist
promoters.threshold.combine$negative %<>% unlist
promoters.threshold.combine$Test %<>% factor
promoters.threshold.genetotals <- select(promoters.threshold.combine, positive, negative, Test) %>% group_by(Test) %>% summarise(sum(sum(positive), sum(abs(negative))))
colnames(promoters.threshold.genetotals)[2] <- "Total.Sites"
promoters.threshold.genetotals$Total.Sites %<>% str_c(" Sites")
promoters.threshold.combine %<>% join(promoters.threshold.genetotals)
promoters.threshold.plot <- gather(promoters.threshold.combine, Direction, Num.Sites, -Test, -Comparison, -Total.Sites)

genes.threshold <- map(thresholds, . %>% map(.x = tab.genes, .f = GetSizes, pval.column = "comb.p.val", log.column = "mean.mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons.keep, length(thresholds)))
genes.threshold$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(genes.threshold)[3] <- "Test"
genes.threshold$Test <- str_c("p < ", genes.threshold$Test)

genes.threshold.fdr <- map(thresholds.fdr, . %>% map(.x = tab.genes, .f = GetSizes, pval.column = "comb.p.adj.fdr", log.column = "mean.mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds.fdr, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons.keep, length(thresholds.fdr)))
genes.threshold.fdr$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(genes.threshold.fdr)[3] <- "Test"
genes.threshold.fdr$Test <- str_c("FDR p < ", genes.threshold.fdr$Test)

genes.threshold.combine <- rbind(genes.threshold, genes.threshold.fdr) 
genes.threshold.combine$positive %<>% unlist
genes.threshold.combine$negative %<>% unlist
genes.threshold.combine$Test %<>% factor
genes.threshold.genetotals <- select(genes.threshold.combine, positive, negative, Test) %>% group_by(Test) %>% summarise(sum(sum(positive), sum(abs(negative))))
colnames(genes.threshold.genetotals)[2] <- "Total.Sites"
genes.threshold.genetotals$Total.Sites %<>% str_c(" Sites")
genes.threshold.combine %<>% join(genes.threshold.genetotals)
genes.threshold.plot <- gather(genes.threshold.combine, Direction, Num.Sites, -Test, -Comparison, -Total.Sites)

sites.threshold <- map(thresholds, . %>% map(.x = tab.sites, .f = GetSizes, pval.column = "diffmeth.p.val", log.column = "mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons.keep, length(thresholds)))
sites.threshold$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(sites.threshold)[3] <- "Test"
sites.threshold$Test <- str_c("p < ", sites.threshold$Test)

sites.threshold.fdr <- map(thresholds.fdr, . %>% map(.x = tab.sites, .f = GetSizes, pval.column = "diffmeth.p.adj.fdr", log.column = "mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds.fdr, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons.keep, length(thresholds.fdr)))
sites.threshold.fdr$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(sites.threshold.fdr)[3] <- "Test"
sites.threshold.fdr$Test <- str_c("FDR p < ", sites.threshold.fdr$Test)

sites.threshold.combine <- rbind(sites.threshold, sites.threshold.fdr) 
sites.threshold.combine$positive %<>% unlist
sites.threshold.combine$negative %<>% unlist
sites.threshold.combine$Test %<>% factor
sites.threshold.genetotals <- select(sites.threshold.combine, positive, negative, Test) %>% group_by(Test) %>% summarise(sum(sum(positive), sum(abs(negative))))
colnames(sites.threshold.genetotals)[2] <- "Total.Sites"
sites.threshold.genetotals$Total.Sites %<>% str_c(" Sites")
sites.threshold.combine %<>% join(sites.threshold.genetotals)
sites.threshold.plot <- gather(sites.threshold.combine, Direction, Num.Sites, -Test, -Comparison, -Total.Sites)

DecidePlot("promoter_threshold_selection", promoters.threshold.plot, "Differentially Methylated Promoters", bar.padding = 500)
DecidePlot("gene_threshold_selection", genes.threshold.plot, "Differentially Methylated Genes", bar.padding = 900)
DecidePlot("site_threshold_selection", sites.threshold.plot, "Differentially Methylated Sites", bar.padding = 10000)

#Gene Annotation
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
vega <- useMart("ENSEMBL_MART_VEGA", dataset = "hsapiens_gene_vega")
ensembl.attributes <- listAttributes(ensembl)
write.xlsx(ensembl.attributes, "./ensembl.attributes.xlsx")

promoters.annotation <- annotation(rnb.norm, type = "promoters")
promoters.annotation$Ensembl.ID <- rownames(promoters.annotation)
SaveRDSgz(promoters.annotation, "./save/promoters_annotation.rda")

promoters.bm.table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description', 'gene_biotype'), filters = 'ensembl_gene_id', values = as.character(promoters.annotation$Ensembl.ID), mart = ensembl)
promoters.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(promoters.bm.table) <- c("Ensembl.ID", "Symbol", "Description", "Gene.Type")

promoters.bm.table.vega <- getBM(attributes = c('ens_gene', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ens_gene', values = as.character(promoters.annotation$Ensembl.ID), mart = vega)
colnames(promoters.bm.table.vega) <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

SaveRDSgz(promoters.bm.table, "./save/promoters.bm.table.rda")
SaveRDSgz(promoters.bm.table.vega, "./save/promoters.bm.table.vega.rda")

genes.annotation <- annotation(rnb.norm, type = "genes")
genes.annotation$Ensembl.ID <- rownames(genes.annotation)
SaveRDSgz(genes.annotation, "./save/genes_annotation.rda")

genes.bm.table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol','description', 'gene_biotype'), filters = 'ensembl_gene_id', values = as.character(genes.annotation$Ensembl.ID), mart = ensembl)
genes.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(genes.bm.table) <- c("Ensembl.ID", "Symbol", "Description", "Gene.Type")

genes.bm.table.vega <- getBM(attributes = c('ens_gene', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ens_gene', values = as.character(genes.annotation$Ensembl.ID), mart = vega)
colnames(genes.bm.table.vega) <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

SaveRDSgz(genes.bm.table, "./save/genes.bm.table.rda")
SaveRDSgz(genes.bm.table.vega, "./save/genes.bm.table.vega.rda")

sites.annotation <- annotation(rnb.norm, type = "sites")
sites.annotation$Illumina.ID <- rownames(sites.annotation)
colnames(sites.annotation) %<>% make.names
SaveRDSgz(sites.annotation, "./save/sites_annotation.rda")

sites.bm.table <- getBM(attributes = c('illumina_human_methylation_450', 'hgnc_symbol', 'description', 'gene_biotype'), filters = 'illumina_human_methylation_450', values = as.character(sites.annotation$Illumina.ID), mart = ensembl)
sites.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(sites.bm.table) <- c("Illumina.ID", "Symbol", "Definition", "Gene.Type")
SaveRDSgz(sites.bm.table, "./save/sites.bm.table.rda")

tab.promoters.annot <- map(tab.promoters, mutate, Ensembl.ID = rownames(promoters.annotation)) %>% 
    map(join, promoters.annotation) %>%
    map(join, promoters.bm.table) %>%
    map(join, promoters.bm.table.vega) %>%
    map(dplyr::select, Ensembl.ID, Symbol:Gene.Type, mean.mean.quot.log2, comb.p.val, comb.p.adj.fdr, combinedRank, Chromosome:Strand, num.sites, CpG:G, Vega.ID:Gene.Type.Vega, entrezID, mean.mean.g1:mean.mean.diff, mean.num.na.g1:mean.nsamples.covg.thresh.g2) 
SaveRDSgz(tab.promoters.annot, "./save/tab.promoters.annot.rda")

tab.genes.annot <- map(tab.genes, mutate, Ensembl.ID = rownames(genes.annotation)) %>% 
    map(join, genes.annotation) %>%
    map(join, genes.bm.table) %>%
    map(join, genes.bm.table.vega) %>%
    map(dplyr::select, Ensembl.ID, Symbol:Gene.Type, mean.mean.quot.log2, comb.p.val, comb.p.adj.fdr, combinedRank, Chromosome:Strand, num.sites, CpG:G, Vega.ID:Gene.Type.Vega, entrezID, mean.mean.g1:mean.mean.diff, mean.num.na.g1:mean.nsamples.covg.thresh.g2) 
SaveRDSgz(tab.genes.annot, "./save/tab.genes.annot.rda")

tab.sites.annot <- map(tab.sites, mutate, Illumina.ID = rownames(sites.annotation)) %>% 
    map(join, sites.annotation) %>%
    map(join, sites.bm.table) %>%
    map(dplyr::select, Illumina.ID, Symbol:Gene.Type, mean.quot.log2, diffmeth.p.val, diffmeth.p.adj.fdr, combinedRank, Chromosome:Strand, CGI.Relation:GC, AddressA:Mismatches.B, SNPs.3:Cross.reactive, mean.g1:mean.diff, max.g1:min.diff, num.na.g1:covg.thresh.nsamples.g2) 
SaveRDSgz(tab.sites.annot, "./save/tab.sites.annot.rda")

names(tab.promoters.annot) <- comparisons.format
names(tab.genes.annot) <- comparisons.format
names(tab.sites.annot) <- comparisons.format

promoters.table.reduce <- map(tab.promoters.annot, dplyr::select, Ensembl.ID:Gene.Type.Vega) %>% 
    lmap(UniqueNames) %>% 
    reduce(join) %>%
    dplyr::select(Ensembl.ID:Gene.Type, dplyr::contains("log2FC"), dplyr::matches("p.val$"), dplyr::contains("p.val.adj"), dplyr::contains("combinedRank"), Chromosome:Gene.Type.Vega) %>%
    arrange(duodenum.crypts_vs._colon.crypts.combinedRank)
promoters.table.reduce[,grepl("p.val|log2FC", colnames(promoters.table.reduce))] %<>% signif(3)
DMWorkbook(promoters.table.reduce, "./promoters.combined.xlsx")

genes.table.reduce <- map(tab.genes.annot, dplyr::select, Ensembl.ID:Gene.Type.Vega) %>% 
    lmap(UniqueNames) %>% 
    reduce(join) %>%
    dplyr::select(Ensembl.ID:Gene.Type, dplyr::contains("log2FC"), dplyr::matches("p.val$"), dplyr::contains("p.val.adj"), dplyr::contains("combinedRank"), Chromosome:Gene.Type.Vega) %>%
    arrange(duodenum.crypts_vs._colon.crypts.combinedRank)
genes.table.reduce[,grepl("p.val|log2FC", colnames(genes.table.reduce))] %<>% signif(3)
DMWorkbook(genes.table.reduce, "./genes.combined.xlsx")

sites.table.reduce <- map(tab.sites.annot, dplyr::select, Illumina.ID:Context) %>% 
    lmap(UniqueNames) %>% 
    reduce(join) %>%
    dplyr::select(Illumina.ID:Definition, dplyr::contains("log2FC"), dplyr::matches("p.val$"), dplyr::contains("p.val.adj"), dplyr::contains("combinedRank"), Chromosome:Context) 
    arrange(crypts_vs._mucosa.combinedRank)
SaveRDSgz(sites.table.reduce, "./save/sites.table.reduce.rda")
sites.table.reduce[,grepl("p.val|log2FC", colnames(sites.table.reduce))] %<>% signif(3)
write.csv(sites.table.reduce, "./sites.combined.csv", row.names = FALSE)

promoters.filter <- filter(promoters.table.reduce, nchar(Symbol) > 0) 
genes.filter <- filter(genes.table.reduce, nchar(Symbol) > 0)

#Upset
format.pvals <- str_c(comparisons.format, ".p.val.adj < 0.01")

top.upset.promoters <- select(promoters.filter, Symbol, dplyr::contains("p.val.adj")) 
top.list.promoters <- map(format.pvals, filter_, .data = top.upset.promoters) %>% map(extract2, "Symbol") 
names(top.list.promoters) <- comparisons.format

top.upset.genes <- select(genes.filter, Symbol, dplyr::contains("p.val.adj")) 
top.list.genes <- map(format.pvals, filter_, .data = top.upset.genes) %>% map(extract2, "Symbol") 
names(top.list.genes) <- comparisons.format

Cairo("upset.promoters", onefile = TRUE, type = "svg", width = 10, height = 6, units = "in", bg = "transparent")
upset(fromList(top.list.promoters), nsets = 4, nintersects = NA, order.by = "freq")
dev.off()

Cairo("upset.genes", onefile = TRUE, type = "svg", width = 10, height = 6, units = "in", bg = "transparent")
upset(fromList(top.list.genes), nsets = 4, nintersects = NA, order.by = "freq")
dev.off()

promoters.submit <-  select(promoters.filter, Ensembl.ID, Symbol, matches("p.val.adj"))
genes.submit <-  select(genes.filter, Ensembl.ID, Symbol, matches("p.val.adj"))

promoters.coding <- filter(promoters.filter, Gene.Type == "protein_coding")
genes.coding <- filter(genes.filter, Gene.Type == "protein_coding")

#Get exact members
promoters.union <- reduce(top.list.promoters, union)
promoters.overlap <- map(top.list.promoters, is.element, el = promoters.union) %>% reduce(cbind) %>% data.frame
colnames(promoters.overlap) <- names(top.list.promoters)
promoters.overlap$Gene <- promoters.union
promoters.dccc.only <- filter(promoters.overlap, duodenum.crypts_vs._colon.crypts == TRUE & duodenum.mucosa_vs._duodenum.crypts == TRUE)
write.xlsx(promoters.dccc.only, "promoters.dccc.only.xlsx")

promoters.dccc.up <- filter(promoters.filter, duodenum.mucosa_vs._duodenum.crypts.p.val.adj < 0.01 & duodenum.crypts_vs._colon.crypts.p.val.adj < 0.01 & duodenum.mucosa_vs._duodenum.crypts.log2FC < 0 & duodenum.crypts_vs._colon.crypts.log2FC > 0)
write.xlsx(promoters.dccc.up, "promoters.dccc.up.xlsx")
promoters.dccc.down <- filter(promoters.filter, duodenum.mucosa_vs._duodenum.crypts.p.val.adj < 0.01 & duodenum.crypts_vs._colon.crypts.p.val.adj < 0.01 & duodenum.mucosa_vs._duodenum.crypts.log2FC > 0 & duodenum.crypts_vs._colon.crypts.log2FC < 0)
write.xlsx(promoters.dccc.down, "promoters.dccc.down.xlsx")

#Top 5 differentially methylated plots
map(comparisons.format, Top5Plot, promoters.coding, promoters.compare, pheno.compare, "promoters")
map(comparisons.format, Top5Plot, genes.coding, genes.compare, pheno.compare, "genes")

#DNAm age
dnam.import <- read_csv("./allsites_beta.output.csv")
dnam.plot <- select(dnam.import, Age, DNAmAge, Cell.Type, Tissue)
dnam.filter <- filter(dnam.plot, Cell.Type == "crypts" | Cell.Type == "mucosa")
dnam.filter$Combined <- str_c(dnam.filter$Tissue, dnam.filter$Cell.Type, sep = " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))
dnam.filter$Difference <- dnam.filter$DNAmAge - dnam.filter$Age 
#age.split <- str_split_fixed(dnam.plot$Chronological.Age, "\\-", 2) %>% data.frame %>% apply(2, as.numeric) %>% apply(1, mean)
#dnam.plot$Mean.Age <- age.split
SaveRDSgz(dnam.filter, "./save/dnam.filter.rda")
write.xlsx(dnam.filter, "dnam.age.xlsx")

p <- ggplot(dnam.filter, aes(Age, DNAmAge, col = Combined)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p <- p + xlab("Chronological Age") + ylab("DNA Methylation Age") + scale_color_discrete(name = "Cell Type")
p <- p + ylim(0,100) + xlim(0,100)
CairoPDF("dnam.plot", width = 7, height = 6)
print(p)
dev.off()

p <- ggplot(dnam.filter, aes(Combined, Difference, color = Combined)) + geom_boxplot(width = 0.5) + geom_jitter(width = 0.5)
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank()) 
p <- p + theme(panel.border = element_rect(color = "black", size = 1))
p <- p + theme(axis.title.x = element_blank()) + ylab("Age Difference") + theme(legend.position = "none")
CairoPDF("diff.plot", width = 7, height = 6, bg = "transparent")
print(p)
dev.off()

kw.statistic <- kruskal.test(Difference ~ Combined, dnam.filter) %>% tidy
dunn.statistic <- posthoc.kruskal.dunn.test(Difference ~ Combined, dnam.filter, p.adjust.method = "fdr")

dnam.bf <- anovaBF(Difference ~ Combined, data.frame(dnam.filter))
dnam.posterior <- posterior(dnam.bf, iterations = 100000) %>% data.frame 
dnam.quantiles <- map(dnam.posterior, quantile, c(0.025, 0.5, 0.975))

dnam.coefs <- select(dnam.posterior, matches("Combined\\.")) 
colnames(dnam.coefs) %<>% str_replace("Combined\\.", "")
estimate.plot <- gather(dnam.coefs, Group, Estimate) 
estimate.plot$Group %<>% str_replace("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))
estimate.plot$Estimate.final <- estimate.plot$Estimate + dnam.quantiles$mu[2]

p <- ggplot(estimate.plot, aes(x = Group, y = Estimate.final, fill = Group)) + geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) + theme_bw()
p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(color = "black", size = 1))
p <- p + theme(plot.background = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) + ylab("Age Difference Estimate") + ggtitle("Age Difference")
p <- p + theme(plot.title = element_text(hjust = 0.5))
CairoPDF("estimate", width = 7, height = 6, bg = "transparent")
print(p)
dev.off()

source("../../code/GO/enrichr.R")
enrichr.terms <- c("GO_Biological_Process_2015", "GO_Molecular_Function_2015", "KEGG_2016", "Reactome_2016", "Human_Gene_Atlas", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up") 
promoters.enrichr <- map(comparisons.format, GetEnrichr, promoters.submit, "0.01", enrichr.terms, "promoters")
genes.enrichr <- map(comparisons.format, GetEnrichr, genes.submit, "0.01", enrichr.terms, "genes")

promoters.up.enrichr <- GetEnrichr("duodenum.crypts_vs._colon.crypts", promoters.dccc.up, "0.01", enrichr.terms, "promoters_dccc_up") 
promoters.down.enrichr <- GetEnrichr("duodenum.crypts_vs._colon.crypts", promoters.dccc.down, "0.01", enrichr.terms, "promoters_dccc_down") 

#Enrichr Plots
dccc.promoters.gobiol.file <- "./enrichr/promoters/duodenum.crypts_vs._colon.crypts/GO_Biological_Process_2015.xlsx"
dccc.promoters.gobiol <- read.xlsx(dccc.promoters.gobiol.file) 
dccc.promoters.gobiol.filter <- FilterEnrichr(dccc.promoters.gobiol)
GetKappaCluster(file_path_sans_ext(dccc.promoters.gobiol.file), dccc.promoters.gobiol.filter, promoters.submit$Symbol)
dccc.promoters.gobiol.final <- slice(dccc.promoters.gobiol.filter, c(2, 4, 6, 17, 5))
dccc.promoters.gobiol.final$Database <- "GO Biological Process"

dccc.promoters.gomole.file <- "./enrichr/promoters/duodenum.crypts_vs._colon.crypts/GO_Molecular_Function_2015.xlsx"
dccc.promoters.gomole <- read.xlsx(dccc.promoters.gomole.file) 
dccc.promoters.gomole.filter <- FilterEnrichr(dccc.promoters.gomole)
GetKappaCluster(file_path_sans_ext(dccc.promoters.gomole.file), dccc.promoters.gomole.filter, promoters.submit$Symbol)
dccc.promoters.gomole.final <- slice(dccc.promoters.gomole.filter, c(1, 23, 6))
dccc.promoters.gomole.final$Database <- "GO Molecular Function"

dccc.promoters.reactome.file <- "./enrichr/promoters/duodenum.crypts_vs._colon.crypts/Reactome_2016.xlsx"
dccc.promoters.reactome <- read.xlsx(dccc.promoters.reactome.file) 
dccc.promoters.reactome.filter <- FilterEnrichr(dccc.promoters.reactome)
GetKappaCluster(file_path_sans_ext(dccc.promoters.reactome.file), dccc.promoters.reactome.filter, promoters.submit$Symbol)
dccc.promoters.reactome.final <- slice(dccc.promoters.reactome.filter, c(4))
dccc.promoters.reactome.final$Database <- "Reactome"

dccc.promoters.enrichr.final <- rbind(dccc.promoters.gobiol.final, dccc.promoters.gomole.final, dccc.promoters.reactome.final)

toptable.dccc.promoters.enrichr <- filter(promoters.filter, !is.na(Symbol) & nchar(Symbol) > 0) %>% select(Symbol, duodenum.crypts_vs._colon.crypts.log2FC)
colnames(toptable.dccc.promoters.enrichr)[2] <- "logFC"
EnrichrPlot(dccc.promoters.enrichr.final, toptable.dccc.promoters.enrichr, "dccc.promoters.enrichr")

dccc.genes.gobiol.file <- "./enrichr/genes/duodenum.crypts_vs._colon.crypts/GO_Biological_Process_2015.xlsx"
dccc.genes.gobiol <- read.xlsx(dccc.genes.gobiol.file) 
dccc.genes.gobiol.filter <- FilterEnrichr(dccc.genes.gobiol)
GetKappaCluster(file_path_sans_ext(dccc.genes.gobiol.file), dccc.genes.gobiol.filter, genes.submit$Symbol)
dccc.genes.gobiol.final <- slice(dccc.genes.gobiol.filter, c(1, 19, 21, 14))
dccc.genes.gobiol.final$Database <- "GO Biological Process"

dccc.genes.gomole.file <- "./enrichr/genes/duodenum.crypts_vs._colon.crypts/GO_Molecular_Function_2015.xlsx"
dccc.genes.gomole <- read.xlsx(dccc.genes.gomole.file) 
dccc.genes.gomole.filter <- FilterEnrichr(dccc.genes.gomole)
GetKappaCluster(file_path_sans_ext(dccc.genes.gomole.file), dccc.genes.gomole.filter, genes.submit$Symbol)
dccc.genes.gomole.final <- slice(dccc.genes.gomole.filter, c(1, 9, 14))
dccc.genes.gomole.final$Database <- "GO Molecular Function"

dccc.genes.reactome.file <- "./enrichr/genes/duodenum.crypts_vs._colon.crypts/Reactome_2016.xlsx"
dccc.genes.reactome <- read.xlsx(dccc.genes.reactome.file) 
dccc.genes.reactome.filter <- FilterEnrichr(dccc.genes.reactome)
GetKappaCluster(file_path_sans_ext(dccc.genes.reactome.file), dccc.genes.reactome.filter, genes.submit$Symbol)
dccc.genes.reactome.final <- slice(dccc.genes.reactome.filter, c(4, 8))
dccc.genes.reactome.final$Database <- "Reactome"

dccc.genes.kegg.file <- "./enrichr/genes/duodenum.crypts_vs._colon.crypts/KEGG_2016.xlsx"
dccc.genes.kegg <- read.xlsx(dccc.genes.kegg.file) 
dccc.genes.kegg.filter <- FilterEnrichr(dccc.genes.kegg)
GetKappaCluster(file_path_sans_ext(dccc.genes.kegg.file), dccc.genes.kegg.filter, genes.submit$Symbol)
dccc.genes.kegg.final <- slice(dccc.genes.kegg.filter, c(5))
dccc.genes.kegg.final$Database <- "KEGG"

dccc.genes.enrichr.final <- rbind(dccc.genes.gobiol.final, dccc.genes.gomole.final, dccc.genes.reactome.final)
toptable.dccc.genes.enrichr <- filter(genes.filter, !is.na(Symbol) & nchar(Symbol) > 0) %>% select(Symbol, duodenum.crypts_vs._colon.crypts.log2FC)
colnames(toptable.dccc.genes.enrichr)[2] <- "logFC"
EnrichrPlot(dccc.genes.enrichr.final, toptable.dccc.genes.enrichr, "dccc.genes.enrichr", plot.width = 10)

#DCCC down
dccc.down.gobiol.file <- "./enrichr/promoters_dccc_down/duodenum.crypts_vs._colon.crypts/GO_Biological_Process_2015.xlsx"
dccc.down.gobiol <- read.xlsx(dccc.down.gobiol.file) 
dccc.down.gobiol.filter <- FilterEnrichr(dccc.down.gobiol, size = 200)
GetKappaCluster(file_path_sans_ext(dccc.down.gobiol.file), dccc.down.gobiol.filter, promoters.dccc.down$Symbol)
dccc.down.gobiol.final <- slice(dccc.down.gobiol.filter, c(2, 4, 6, 17, 5))
dccc.down.gobiol.final$Database <- "GO Biological Process"

dccc.down.gomole.file <- "./enrichr/promoters_dccc_down/duodenum.crypts_vs._colon.crypts/GO_Molecular_Function_2015.xlsx"
dccc.down.gomole <- read.xlsx(dccc.down.gomole.file) 
dccc.down.gomole.filter <- FilterEnrichr(dccc.down.gomole)
GetKappaCluster(file_path_sans_ext(dccc.down.gomole.file), dccc.down.gomole.filter, promoters.dccc.down$Symbol)
dccc.down.gomole.final <- slice(dccc.down.gomole.filter, c(1, 23, 6))
dccc.down.gomole.final$Database <- "GO Molecular Function"

dccc.down.reactome.file <- "./enrichr/promoters_dccc_down/duodenum.crypts_vs._colon.crypts/Reactome_2016.xlsx"
dccc.down.reactome <- read.xlsx(dccc.down.reactome.file) 
dccc.down.reactome.filter <- FilterEnrichr(dccc.down.reactome)
GetKappaCluster(file_path_sans_ext(dccc.down.reactome.file), dccc.down.reactome.filter, promoters.dccc.down$Symbol)
dccc.down.reactome.final <- slice(dccc.down.reactome.filter, c(4))
dccc.down.reactome.final$Database <- "Reactome"

dccc.down.kegg.file <- "./enrichr/promoters_dccc_down/duodenum.crypts_vs._colon.crypts/KEGG_2016.xlsx"
dccc.down.kegg <- read.xlsx(dccc.down.kegg.file) 
dccc.down.kegg.filter <- FilterEnrichr(dccc.down.kegg)
GetKappaCluster(file_path_sans_ext(dccc.down.kegg.file), dccc.down.kegg.filter, promoters.dccc.down$Symbol)
dccc.down.kegg.final <- slice(dccc.down.kegg.filter, c(4))
dccc.down.kegg.final$Database <- "KEGG"

dccc.down.enrichr.final <- rbind(dccc.down.gobiol.final, dccc.down.gomole.final, dccc.down.reactome.final)

#Phenotype table
pheno.table <- table(droplevels(pheno.export$Tissue), droplevels(pheno.export$Cell.Type)) %>% as.data.frame.matrix %>% select(crypts, mucosa, enteroid, spheroid)
colnames(pheno.table) %<>% capitalize
rownames(pheno.table) %<>% capitalize
pheno.table$Total <- rowSums(pheno.table)
pheno.total <- colSums(pheno.table)
pheno.table %<>% rbind(Total = pheno.total)
write.csv(pheno.table, "pheno_table.csv")

#Annotate CpGs from methylation clock
clock.bm.table <- getBM(attributes = c('illumina_human_methylation_450', 'hgnc_symbol', 'description', 'gene_biotype'), filters = 'illumina_human_methylation_450', values = as.character(clock.annotation$CpGmarker), mart = ensembl)
clock.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(clock.bm.table) <- c("Illumina.ID", "Symbol", "Definition", "Gene.Type")
write.xlsx(clock.bm.table, "clock.bm.table.xlsx")

clock.annotation <- read_csv("./13059_2013_3156_MOESM23_ESM.csv") %>% slice(-1) #Drop intercept row from table
InfiniumMethylation <- features(FDb.InfiniumMethylation.hg19) #Load InfiniumMethylation database
hm450 <- get450k() #get annotations for 450k
clock.transcripts <- getNearestTranscript(hm450[clock.annotation$CpGmarker]) #get nearest transcript for the CpGs of interest
clock.transcripts$Illumina.ID <- rownames(clock.transcripts) #Add Illumina ID column
write.xlsx(clock.transcripts, "clock.transcripts.xlsx") #Save to .xlsx

gata4.meth <- promoters.compare["ENSG00000136574",]
gata4.df <- data.frame(Tissue = pheno(rnb.compare)$Tissue, Cell.Type = pheno(rnb.compare)$Cell.Type, Methylation = gata4.meth) 
gata4.df$Combined <- str_c(gata4.df$Tissue, gata4.df$Cell.Type, sep = " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))

p <- ggplot(gata4.df, aes(x = Combined, y = Methylation, color = Combined)) + geom_boxplot() + geom_jitter() + theme_bw()
p <- p + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1)
CairoPDF("GATA4_methylation", height = 4, width = 4, bg = "transparent")
print(p)
dev.off()

gata5.meth <- promoters.compare["ENSG00000130700",]
gata5.df <- data.frame(Tissue = pheno(rnb.compare)$Tissue, Cell.Type = pheno(rnb.compare)$Cell.Type, Methylation = gata5.meth) 
gata5.df$Combined <- str_c(gata5.df$Tissue, gata5.df$Cell.Type, sep = " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts"))

p <- ggplot(gata5.df, aes(x = Combined, y = Methylation, color = Combined)) + geom_boxplot() + geom_jitter() + theme_bw()
p <- p + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1)
CairoPDF("GATA5_methylation", height = 4, width = 4, bg = "transparent")
print(p)
dev.off()


