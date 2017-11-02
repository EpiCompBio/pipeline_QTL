#############################
# To be run after 03 differential gene expression analysis of array data

#############################


#############################
# Illumina gene expression microarray analysis steps (Ritchie et al. 2014 PLoS Comp Bio; 
# Ritchie et al. 2015 Nucleic Acids Res.):

# 1) Data acquisition
# 2) Preprocessing and quality assessment

# 3) Background correction, normalisation, transformation and filtering

# 4) Experimental design matrix specification
# 5) Descriptive statistics
# 6) Differential gene expression analysis

# 7) Higher level analysis: Gene ontology, co-expression, gene set analyses
# 8) Integration with other data

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")
# setwd('/Users/antoniob/Desktop/BEST_D_12_FEB.DIR/results_1/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_higher_level",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
#load('R_session_saved_image_read_and_QC.RData', verbose=T)
#load('R_session_saved_image_diff_expression_full.RData', verbose=T)
load('R_session_saved_image_diff_expression_3.RData', verbose=T)

# To load multiple .RData files:
#rdata_filenames <- c('.RData')
#lapply(rdata_filenames, load, .GlobalEnv)


# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_higher_level', '.RData', sep='')
R_session_saved_image_full <- paste('R_session_saved_image_higher_level_full', '.RData', sep='')
#############################


#############################
## Update packages if necessary and load them:
#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#install.packages('GRDAVIDWebService')
# biocLite("RDAVIDWebService")

library(limma)
library(illuminaHumanv4.db)
library(reshape2)
library(gridExtra)
library(plyr)
library(ReactomePA)
library(clusterProfiler)
library(GOstats)
library(ggplot2)
library(AnnotationHub)
library(RDAVIDWebService)

# Get additional functions needed:
source('/Users/antoniob/Desktop/BEST_D_12_FEB.DIR/scripts/gene_expression_functions.R')
#source('/ifs/devel/antoniob/projects/BEST-D/gene_expression_functions.R')
#############################


#############################

# 7) Integration with other data
# a) eQTL analysis
# b) Disease variant overlap
# c) Pathway analysis

# Check results from differential expression analysis:
# TO DO: pass toptable name or read in as file provided:
topTable_results <- topTable_pairing_joint_treated
# Enforce row names as these may have been lost:
row.names(topTable_results) <- topTable_results$probe_ID
dim(topTable_results)
head(topTable_results)
summary(topTable_results)

# Get top genes and ENTREZ IDs if needed (may already be if topTable was annotated):
# TO DO: pass adj p value as parameter:
top_genes <- get_EntrezIDs_top_genes(topTable_results, 0.05)
head(top_genes)
length(top_genes)

# Get ENTREZ IDs universe for background comparisons:
filtered_illumina_ids  <- as.character(row.names(topTable_results))
head(filtered_illumina_ids)
length(filtered_illumina_ids)

ENTREZID_universe  <- unlist(mget(filtered_illumina_ids, illuminaHumanv4ENTREZID, ifnotfound = NA))
head(ENTREZID_universe)
summary(ENTREZID_universe)
str(ENTREZID_universe)
count(!is.na(ENTREZID_universe))

# Gene ontology analysis: goana, geneSetTest, camera, roast, romer
goana_top_genes <- goana(top_genes, universe=ENTREZID_universe, species='Hs')
head(goana_top_genes)
topG0_top_genes <- topGO(goana_top_genes, number=Inf)
head(topG0_top_genes)
dim(topG0_top_genes)

# Get different processes from GO:
topG0_top_genes_BP <- topGO(goana_top_genes, number=50, ontology='BP')
topG0_top_genes_CC <- topGO(goana_top_genes, number=50, ontology='CC')
topG0_top_genes_MF <- topGO(goana_top_genes, number=50, ontology='MF')
head(topG0_top_genes_BP)
head(topG0_top_genes_MF)
head(topG0_top_genes_CC)

barplot(topG0_top_genes_BP$P.DE, names.arg=topG0_top_genes_BP$Term, horiz=T)
#barcodeplot(topG0_top_genes_BP$P.DE, names.arg=topG0_top_genes_BP$Term, horiz=T)

png('barplot_topGO_BP.png', width = 12, height = 12, units = 'in', res = 300)
ggplot(topG0_top_genes_BP, aes(y=DE, x=Term)) + geom_point(stat='identity', aes(colour = P.DE, size=N)) + coord_flip()
dev.off()

png('barplot_topGO_MF.png', width = 12, height = 12, units = 'in', res = 300)
ggplot(topG0_top_genes_MF, aes(y=DE, x=Term)) + geom_point(stat='identity', aes(colour = P.DE)) + coord_flip()
dev.off()

png('barplot_topGO_CC.png', width = 12, height = 12, units = 'in', res = 300)
ggplot(topG0_top_genes_CC, aes(y=DE, x=Term)) + geom_point(stat='identity', aes(colour = P.DE)) + coord_flip()
dev.off()

# TO DO: update limma to new version that includes kegga()

# GOstats package:
# ?GOHyperGParams
# goHyperG <- hyperGTest(p=unique(topLocusLinkIDs), lib="hgu95av2", what="MF")
# ?hyperGTest


# ReactomePA analysis:
# TO DO: PA_reactome <- enrichPathway(gene=top_genes_5_ENTREZID, minGSSize=25, pvalueCutoff=0.05, qvalueCutoff=0.05, readable=T, organism='human') #universe=ENTREZID_universe,
#head(summary(PA_reactome))
#barplot(PA_reactome, showCategory=15)

# clusterProfiler analysis:
# http://www.bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.pdf

clusterP_enrich_MF <- enrichGO(gene=top_genes_ENTREZID, universe=ENTREZID_universe, pAdjustMethod='BH', ont='MF', minGSSize=10,
                     pvalueCutoff=0.01, qvalueCutoff=0.01, readable=T, organism='human')
head(summary(clusterP_enrich_MF), 200)
tail(summary(clusterP_enrich_MF))

png('barplot_clusterP_enrich_MF.png', width = 12, height = 12, units = 'in', res = 300)
barplot(clusterP_enrich_MF, drop=TRUE, showCategory=20)
dev.off()

png('enrichMap_clusterP_enrich_MF.png', width = 12, height = 12, units = 'in', res = 300)
enrichMap(clusterP_enrich_MF)#, fixed = FALSE)
dev.off()


clusterP_enrich_BP <- enrichGO(gene=top_genes_ENTREZID, universe=ENTREZID_universe, pAdjustMethod='BH', ont='BP', minGSSize=10,
                               pvalueCutoff=0.01, qvalueCutoff=0.01, readable=T, organism='human')
head(summary(clusterP_enrich_BP), 100)

png('barplot_clusterP_enrich_BP.png', width = 12, height = 12, units = 'in', res = 300)
barplot(clusterP_enrich_BP, drop=TRUE, showCategory=20)
dev.off()

png('enrichMap_clusterP_enrich_BP.png', width = 12, height = 12, units = 'in', res = 300)
enrichMap(clusterP_enrich_BP)#, fixed=F)
dev.off()



clusterP_enrich_CC <- enrichGO(gene=top_genes_ENTREZID, universe=ENTREZID_universe, pAdjustMethod='BH', ont='CC', minGSSize=10,
                               pvalueCutoff=0.01, qvalueCutoff=0.01, readable=T, organism='human')
head(summary(clusterP_enrich_CC), 100)

png('barplot_clusterP_enrich_CC.png', width = 12, height = 12, units = 'in', res = 300)
barplot(clusterP_enrich_CC, drop=TRUE, showCategory=20)
dev.off()

png('enrichMap_clusterP_enrich_CC.png', width = 12, height = 12, units = 'in', res = 300)
enrichMap(clusterP_enrich_CC)#, fixed = FALSE)
dev.off()

# clusterProfiler::cnetplot(clusterP_enrich, categorySize="pvalue", foldChange='geneList')
# ?clusterProfiler::cnetplot

# TO DO: errors:
#clusterP_gse <- gseGO(ENTREZID_universe, pAdjustMethod='BH', pvalueCutoff=0.01, nPerm=100, organism='human')
#head(summary(clusterP_gse))

clusterP_kegg <- enrichKEGG(gene=top_genes_ENTREZID, universe=ENTREZID_universe, pvalueCutoff=0.05, qvalueCutoff=0.2, pAdjustMethod='BH', organism='human')
head(summary(clusterP_kegg))
barplot(clusterP_kegg)

# TO DO: these error:
#clusterP_gse_kegg <- gseKEGG(gene=top_genes_1_ENTREZID, pvalueCutoff=0.05, nPerm=100, organism='human')
#head(summary(clusterP_kegg))

# clusterP_DAVID <- enrichDAVID(gene=top_genes_ENTREZID, idType="ENTREZ_GENE_ID", listType='Gene', 
#                               species = 'human', pvalueCutoff=0.05, qvalueCutoff=0.2, pAdjustMethod='BH')
# annotation="KEGG_PATHWAY")
# head(summary(clusterP_DAVID))


#############################


#############################
# Check AnnotationHub, download required sources, process tables, etc., then cross with experiment data.
# See:
# https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html
# https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub-HOWTO.html

#############################


#############################
# Check Hai Feng's packages:
# XOR:
# https://github.com/hfang-bristol/XOR/blob/master/inst/INSTALLATION.md
# d-net:
# http://supfam.org/dnet/


#############################


#############################
#The end:
# TO DO: remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes)))

#rm()

# To save R workspace with all objects to use at a later time:
# Also consider dump(), dput() and dget(), see http://thomasleeper.com/Rcourse/Tutorials/savingdata.html
save.image(file=R_session_saved_image_full, compress='gzip')

# Or save specific objects:
# objects_to_save <- (c('normalised_expressed', 'membership_file_cleaned', 'FAILED_QC'))
# save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: collect data, plots and tables for a report
#############################
