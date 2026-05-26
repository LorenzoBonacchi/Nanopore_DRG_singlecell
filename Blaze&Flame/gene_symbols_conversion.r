
convert_ENSGID_to_geneSymbol <- function(gene_count_matrix_path, 
                                         id_symbol_df = isoform_gene_dict, 
                                         output_file,
                                         return_df = FALSE) {
  
  # Load the reference dictionary we made earlier - select gene-level cols
  id_symbol_df <- as_tibble(id_symbol_df) %>%
    dplyr::select(gene_id, gene_symbol)
  
  # Load the data object with ENSGID row names
  gene_count_matrix <- fread(gene_count_matrix_path, header = TRUE)
  colnames(gene_count_matrix)[1] <- "gene_id"
  
  # Replace ENSGIDs with gene symbols in original flames gene-level count matrix
  formatted_gene_count_matrix <- gene_count_matrix %>%
    merge(id_symbol_df, by.x = 'gene_id', by.y = 'gene_id') %>%   # Add gene symbol information
    distinct(gene_symbol, .keep_all = TRUE) %>%   # Remove duplicates based on gene symbol
    dplyr::select(-gene_id) %>%   # Remove the ENSGID column
    column_to_rownames(var = "gene_symbol")   # use the gene symbols we added as rownames
  
  # Write out the processed data frame
  fwrite(formatted_gene_count_matrix, 
            output_file, 
            row.names = TRUE)
  
  # Return the processed count matrix for further use if needed
  if(return_df){
    return(formatted_gene_count_matrix)
  }
}



# convert Gene_id to gene symbol for gene counts
convert_ENSGID_to_geneSymbol(
  gene_count_matrix_path = "./data/FLAMES_out/gene_count.csv",
  output_file = "./output_files/counts/geneSymbol_gene_count.csv"
)

# convert Gene_id to gene symbol for background counts
convert_ENSGID_to_geneSymbol(
  gene_count_matrix_path = "./data/background/gene_count.csv",
  output_file = "./output_files/counts/background_geneSymbol_gene_count.csv"
)