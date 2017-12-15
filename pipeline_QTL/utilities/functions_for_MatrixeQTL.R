###################
# General functions for processing genotype, gene expression (or other quantitative molecular measure, 
# eg metabolomes) and covariates (phenoytpes) files.
# Antonio Berlanga-Taylor
# 06 March 2016

###################


###################
# Functions for 01_eQTL script which is project specific but processing and cleaning files can use these 'generic' functions.
# Aim to run scripts 'per file', such as 00, 02, 04 eQTL scripts. 01 is project specific but write more general use functions
# that could be applied to 018 and Airwave.
# Get master IDs:

# Function to transpose file (usually phenotype file) and name column headers to arbitrary column.
transpose_file <- function(infile, column_number_for_rownames) {
  if (class(infile)[1] == 'data.frame') {
    infile <- as.data.table(infile)
  }
  # Get row names to use as column names for file to transpose:
  row_names <- infile[[colnames(infile)[column_number_for_rownames]]] # Double brackets to return vector
  # Get column names to use as row names:
  col_names <- as.data.frame(colnames(infile)[-1])
  # Delete column as will turn numerics to factors when transposing:
  infile <- infile[, -column_number_for_rownames, with = F]
  # Tranpose and convert to data frame:
  infile_transposed <- transpose(infile) # data.table transpose() turns factors into character
  setnames(infile_transposed, colnames(infile_transposed), as.character(row_names))
  infile_transposed[, ':='(rownames=col_names[, 1])]
  return(infile_transposed)
  }


# Function to change plink column names when output as FID_IID:
split_plink_IDs <- function(geno_infile){
  # Plink outputs FID and IID as FID_IID, so new headers are needed to match 
  # the gene expression (or other molecular data) headers:
  geno_new_colnames <- strsplit(colnames(geno_infile), split='_', fixed=TRUE)
  geno_new_colnames <- matrix(unique(unlist(geno_new_colnames), ncol=1, byrow=T))
  names(geno_infile) <- geno_new_colnames
  return(geno_infile)
}

# Function to subset files based on one variable (column header) and extract second variable values as list:
get_subset_IDs <- function(master_file, column_name_1, value_for_subset, column_name_2){
  subset_file <- master_file[which(master_file[, as.character(column_name_1)] == value_for_subset), ]
  subset_IDs <- na.omit(subset_file[, as.character(column_name_2)])
  return(subset_IDs)
}

# Function to subset files based on two variables (column headers) and extract third variable values as list:
get_subset_IDs_two_vars <- function(master_file, column_name_1, value_for_subset_1, column_name_2, value_for_subset_2, column_name_3){
  subset_file <- master_file[which(master_file[, as.character(column_name_1)] == value_for_subset_1), ]
  subset_file <- subset_file[which(master_file[, as.character(column_name_2)] != value_for_subset_2), ]
  subset_IDs <- na.omit(subset_file[, as.character(column_name_3)])
  return(subset_IDs)
}


# Function to subset using data.table:
get_subset_dt <- function(master_data_table, subset_IDs, extra_cols_to_keep){
  subset_file <- master_data_table[, which(names(master_data_table) %in% as.character(subset_IDs))]
  subset_file <- master_data_table[, .SD, .SDcols = c(extra_cols_to_keep, which(names(master_data_table) %in% as.character(subset_IDs)))]
  return(subset_file)
}

# Function to order files:


# Function to match files:

# Function to save plots to disk as pdf:
dev_to_pdf <- function(filename, pdf.width = 7, pdf.height = 5) {
  dev.copy2pdf(file=filename, width = pdf.width, height = pdf.height)
}

# Function to get the name of a variable as a string:
get_var_name <- function(var1){
  deparse(substitute(var1))
}

## Add annotations:
get_illumina_annot <- function(infile, column_probes){
  # illumina_ids  <-  as.character(infile[, column_probes, with = F]) # data.table object
  if (dim(infile)[1] != 0) {
  infile <- as.data.frame(infile)
  illumina_ids  <-  as.character(infile[, column_probes])
  probes_by_ENTREZID_RE  <- unlist(mget(illumina_ids, illuminaHumanv4ENTREZREANNOTATED, ifnotfound = NA))
  # TO DO: this function takes a long time and in some cases generates many duplicates, see Venn_two_sets.R
  probes_by_ENTREZID_RE_df <- do.call(rbind, lapply(probes_by_ENTREZID_RE, data.frame, stringsAsFactors=FALSE))
  probes_by_ENTREZID_RE_df$Probe_ID <- row.names(probes_by_ENTREZID_RE_df)
  colnames(probes_by_ENTREZID_RE_df)[1] <- 'ENTREZ_ID'
  
  probes_by_symbol <- unlist(mget(illumina_ids, illuminaHumanv4SYMBOLREANNOTATED, ifnotfound = NA))
  probes_by_symbol_df <- do.call(rbind, lapply(probes_by_symbol, data.frame, stringsAsFactors=FALSE))
  probes_by_symbol_df$Probe_ID <- row.names(probes_by_symbol_df)
  colnames(probes_by_symbol_df)[1] <- 'gene'
  
  df_annot <-  merge(infile, probes_by_symbol_df)
  df_annot <-  merge(df_annot, probes_by_ENTREZID_RE_df)
  # df_annot <- data.table(df_annot)
  return(df_annot)
  } else {
    print('Data frame supplied is empty')
    return(infile)
    }
}


# Functions to compare two files and obtain a data frame with shared and unique values
# Shared returns SNP-probe pairs where shared 'column' appear, so elements may be duplicated in column:
getShared <- function(data1, data2, column, order_col) {
  # Shared:
  data1 <- as.data.frame(data1)
  data2 <- as.data.frame(data2)
  shared <- data1[which(as.character(data1[, column]) %in% as.character(data2[, column])), ]
  shared <- shared[order(shared[, order_col]), ]
  return(shared)
}

# Private returns SNP-probe pairs where elements of 'column' only appear in data1:
getPrivate <- function(data1, data2, column, order_col) {
  # Unique
  data1 <- as.data.frame(data1)
  data2 <- as.data.frame(data2)
  private_to_data1 <- data1[which(!as.character(data1[, column]) %in% as.character(data2[, column])), ]
  private_to_data1 <- private_to_data1[order(private_to_data1[, order_col]), ]
  return(private_to_data1)
}

# ## Function to write tsv files to disk:
# write_table <- function(table_to_write){
#   table_to_write <- write.table(table_to_write, sprintf('%s.tsv', deparse(substitute(table_to_write))), 
#             sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)
#   return(table_to_write)
# }
# 
# write_DT_tsv <- function(var2){
#   var_name <- get_var_name(var2)
#   # This generates an empty first header with row names as numbers, use fread with drop = 1 afterwards or cut:
#   write.table(var2, paste(var_name, '.txt', sep = ''), sep='\t', quote = FALSE, col.names = NA)
# }

###################
