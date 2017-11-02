####################
# XGR Downstream analysis
# Hai's package
# 19 May 2016
####################

####################
# See:
# https://github.com/hfang-bristol/XGR
# https://cran.r-project.org/web/packages/XGR/XGR.pdf
# http://galahad.well.ox.ac.uk:3020/XGR_vignettes.html
# http://galahad.well.ox.ac.uk:3020/enricher/genes
####################

####################
options(echo = TRUE)
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D_03_MAR.DIR/results_to_share/BEST-D_22_Apr_2016/')
# setwd('/Users/antoniob/Desktop/BEST_D_03_MAR.DIR/GAT_backgrounds/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_",Sys.Date(),".txt", sep=""))
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))
getwd()

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
# load('R_session_saved_image_XGR.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_XGR','.RData', sep='')
R_session_saved_image
####################

####################
# Libraries:
library(XGR)

# Get additional functions needed:
source('gene_expression_functions.R')
# source('/ifs/devel/antoniob/projects/BEST-D/gene_expression_functions.R')
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Gene or SNP files:
hits_file <- as.character(args[1])
# hits_file <- 'eGenes_cis_tx_fdr5_reQTLs_annot_all_Tx_joint_cis.txt'
# hits_file <- 'Diff_exp_BESTD_FDR10.txt'
# hits_file <- 'me_25_covariates.eQTL'
  
background_file <- as.character(args[2])
# background_file <- 'background_genes_BESTD_expressed.txt'
# background_file <- 'cut_chr22_Airwave_CPMG_Plasma_clean_SNPs_autosome.geno_matched.tsv'

# To test enrichment against a particular set
# annotation_file <- as.character(args[3])
# annotation_file <- 'VD_exp_genes_Ramagopalan_2010.txt'

# Ontology term to test against:
# ontology_term <- as.character(args[3])
####################

####################
# Read data:
# hits_data <- fread(hits_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
hits_data <- read.csv(hits_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
hits_data <- as.vector(hits_data$V1)
class(hits_data)
str(hits_data)
hits_data

background_data <- read.csv(background_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
background_data <- as.vector(background_data$V1)
str(background_data)

# TO DO/clean up
# annotation_data <- read.csv(annotation_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
# str(annotation_data)
# head(annotation_data)
####################

####################
ontology_terms <- c("GOBP", "GOMF", "GOCC",
             "PS", "PS2", "SF", "DO", "HPPA", "HPMI", "HPCM", "HPMA", "MP",
             "MsigdbH", "MsigdbC1", "MsigdbC2CGP", "MsigdbC2CPall", "MsigdbC2CP",
             "MsigdbC2KEGG", "MsigdbC2REACTOME", "MsigdbC2BIOCARTA", "MsigdbC3TFT", "MsigdbC3MIR",
             "MsigdbC4CGN", "MsigdbC4CM", "MsigdbC5BP", "MsigdbC5MF", "MsigdbC5CC",
             "MsigdbC6", "MsigdbC7", "DGIdb")

# Run enrichment analysis and save to file as a loop for all ontology terms:
for (i in ontology_terms){
  ontology_term <- i
  print(ontology_term)
  run_XGR_xEnricherGenes(hits_data, background_data, ontology_term, hits_file)
}

# view enrichment results for the top significant terms
# View(xEnrichViewer(enrichment_result))

# Visualise the top 10 significant terms in the ontology hierarchy
# Set an ontology term:
ontology_term <- 'MsigdbC2CPall'

# Run enrichment analysis for a particular ontology:
run_XGR_xEnricherGenes(hits_data, background_data, ontology_term, hits_file)

# Load ig.XXX (an object of class "igraph" storing as a directed graph):
g <- xRDataLoader(RData=sprintf('ig.%s', ontology_term))
g
nodes_query <- names(sort(enrichment_result$adjp)[1:10])
nodes.highlight <- rep("red", length(nodes_query))
names(nodes.highlight) <- nodes_query
subg <- dnet::dDAGinduce(g, nodes_query)

# Colour-code terms according to adjusted p-values (taking the form of 10-based negative logarithm):
png(sprintf("%s_%s_ontology_hierarchy.png", hits_file, ontology_term), width = 12, height = 12, units = 'in', res = 300)
dnet::visDAG(g=subg, data=-1*log10(enrichment_result$adjp[V(subg)$name]),
             node.info="both", zlim=c(0,2), node.attrs=list(color=nodes.highlight))
# color-code terms according to the z-scores
dnet::visDAG(g=subg, data=enrichment_result$zscore[V(subg)$name], node.info="both",
             colormap="darkblue-white-darkorange",
             node.attrs=list(color=nodes.highlight))
dev.off()
####################

####################
# # TO DO: run SNPs as separate script
# # Test SNPs:
# ontology_term <- 'EF'
# enrichment_result <- xEnricherSNPs(data = hits_data, ontology = ontology_term)
# # view enrichment results for the top significant terms
# xEnrichViewer(enrichment_result)
# # save enrichment results to file:
# save_enrichment <- xEnrichViewer(enrichment_result, top_num=length(enrichment_result$adjp), sortBy="adjp", details=TRUE)
# output_enrichment <- data.frame(term=rownames(save_enrichment), save_enrichment)
# utils::write.table(output_enrichment, file=sprintf("%s_%s_enrichments.txt", hits_file, ontology_term), sep="\t", row.names=FALSE)
####################

####################
# TO DO/clean up
# # Test overlap between gene lists using XGR:
# # xEnricherYours requires dataframes as files, simply provide file names:
# enrichment_result <- xEnricherYours(data.file = hits_file, annotation.file = annotation_file, background.file = background_file)
# xEnrichViewer(enrichment_result)
# save_enrichment <- xEnrichViewer(enrichment_result, top_num=length(enrichment_result$adjp), sortBy="adjp", details=TRUE)
# output_enrichment <- data.frame(term=rownames(save_enrichment), save_enrichment)
# utils::write.table(output_enrichment, file=sprintf("%s_%s_overlap.txt", hits_file, background_file), sep="\t", row.names=FALSE)
####################

####################

####################

####################
# Run similarity analysis:


####################

####################
# Run network analysis:
####################

####################
# Visualise results:
# xCircos()
# xVisNet()

####################

####################
# Interpretation:
# 
# 
####################


####################
# Write results to disk:

####################

####################
# The end:
# Remove objects that are not necessary to save:

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
# save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next: run the script for xxx.
####################
