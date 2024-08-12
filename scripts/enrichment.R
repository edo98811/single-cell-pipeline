# Author: Edoardo Filippi
# mail: efilippi@uni-mainz.de

#' @details  Enrichment analysis analysis using one of the following methods: gsea, ora, panther, enrichr
#
# Arguments:
#' @param wgcna_folder
#' @param module_genes
#' @param module_threshold
#' @param wgcna_module
#' @param wgcna_exclude
#' @param cluster
#' @param count_threshold
#' @param which 
#' @param all_genes
#' @param subset
#' @param minGSSize
#' @param maxGSSize
#' @param organism
#' @param num_tries
#' @param scoring
#' @return result (list of results objects organized)

#' @examples 



enrichment_analysis <- function(name, markers_path, 
                                count_threshold = 0, wgcna_folder = FALSE, 
                                module_threshold = 500, wgcna_exclude = FALSE,
                                wgcna_module = FALSE, modules_significance_table = FALSE, 
                                cluster = FALSE, ...) {
  # source followed: (to add)
  # go options: BP, CC, MF, ALL

  # Check types of arguments
  if (!is.character(name) || length(name) != 1) stop("name argument must be a single character string")
  if (!is.character(markers_path) || length(markers_path) != 1) stop("markers_path argument must be a single character string")
  if (!is.integer(count_threshold) || length(count_threshold) != 1) stop("count_threshold argument must be a single numeric value")
  if (!is.integer(module_threshold) || length(module_threshold) != 1) stop("module_threshold argument must be a single numeric value")
  if (!((is.logical(wgcna_exclude) && isFALSE(wgcna_exclude)) || is.character(wgcna_exclude))) stop("wgcna_exclude argument must be either FALSE or a character vector")
  if (!((is.logical(wgcna_module) && isFALSE(wgcna_module)) || is.character(wgcna_module))) stop("wgcna_module argument must be either FALSE or a character vector")
  if (!((is.logical(cluster) && isFALSE(cluster)) || is.numeric(cluster))) stop("cluster argument must be either FALSE or a numerical vector or single numerical value")

  
  # Set up   output dir
  output_dir <- set_up_output(paste0(output_folder, "GSEA_", name, "/"), message)
  
  # Messages
  cat("Running Enrichment analysis...\n")
  message(paste0("Parameters: name: ", name,
                 " - count_threshold: ", count_threshold,
                 " - markers_path: ", markers_path, 
                 " - module_threshold: ", module_threshold,
                 " - wgcna_folder: ", wgcna_folder,
                 " - wgcna_exclude: ", paste(wgcna_exclude, collapse = ", "),
                 " - wgcna_modules: ", paste(wgcna_module, collapse = ", "),
                 " - cluster: ", paste(cluster, collapse = ", ")))
  
  # Check useful packages
  check_packages(c("Seurat", "tools"))
  library("Seurat")
  library("rlang")

  if (!isFALSE(modules_significance_table) || !isFALSE(wgcna_folder)) wgcna_module <- .select_mdoules_to_enrich(filename = paste0(output_folder, wgcna_folder, modules_significance_table))
  
  # List excel files
  excel_files <- list.files(paste0(output_folder, markers_path), pattern = "\\.xlsx$", full.names = TRUE)
  if (length(excel_files) == 0) stop(paste0("No excel files found in directory ", paste0(output_folder, markers_path)))
  
  # Loop through each Excel file, define needed objects
  for (file in excel_files) {
    
    # Run only for desired file
    if (is.numeric(cluster)) {
      if (!grepl(paste0("expressed_markers_", cluster), file)) {
        message("Skipping file: ", file)
        next
      }
      else message("Running for cluster ", cluster)
    }
    
    # Preparing data
    message("Loading gene data...")
    gene_rankings <- prepare_genes(file, count_threshold, ...)
    
    # Enrichment by wgcna module (enter in it if wgcna exclude is empty, if it is not default to it)
    if (is.character(wgcna_folder) && nchar(wgcna_folder) > 1 && !is.character(wgcna_exclude))  {
      
      # Messages
      message("Running enrichment by wgcna module")
      
      # Load module names
      sheet_names <- excel_sheets(paste0(output_folder, wgcna_folder, "module_genes.xlsx" ))
      
      # If only a subset was given
      if (length(wgcna_module) > 0) sheet_names <- sheet_names[sheet_names %in% wgcna_module]
      
      # Iterate though modules
      for (sheet_name in sheet_names) {
        
        # Load list of genes for module to explore
        module_genes <- c(read_excel(paste0(output_folder, wgcna_folder, "module_genes.xlsx" ), sheet = sheet_name))[[1]]
        
        # Name of this analysis
        analysis_name <- paste0(file_path_sans_ext(basename(file)), "_", sheet_name)
        
        # If the module contains enough genes use GSEA otherwise panther
        if (length(module_genes) < module_threshold) {
          message("Running for module: ", sheet_name)
          
          # Filter the gene ranking list with the module genes
          gene_rankings_subset <- gene_rankings[names(gene_rankings) %in% module_genes]

          # Run GSEA and enrichr on genes 
          result <- run_enrichment(which = c("ora" ,"panther", "enrichr"), output_dir, analysis_name, 
                                   all_genes=gene_rankings, # universe
                                   subset=select_genes_for_enrich(gene_rankings_subset, ...), ...)
          
        }
        else {
          message("Running enrichment for module: ", sheet_name)
          
          # Filter the gene ranking list with the module genes
          gene_rankings_subset <- gene_rankings[names(gene_rankings) %in% module_genes]

          # Run GSEA and enrichr on genes 
          result <- run_enrichment(which = c("gsea" ,"panther", "enrichr"), output_dir, analysis_name, 
                                   all_genes=gene_rankings_subset, 
                                   subset=select_genes_for_enrich(gene_rankings_subset, ...), ...)
        }

      }
      
    } 
    
    # Enrichment excluding wgcna module
    else if(is.character(wgcna_folder) && nchar(wgcna_folder) > 1 && is.character(wgcna_exclude))  {
      
      # Message
      message("Running enrichment excluding wgcna modules: ", unlist(wgcna_exclude))
      
      # Load module names
      sheet_names <- excel_sheets(paste0(output_folder, wgcna_folder, "module_genes.xlsx" ))
      
      # Iterate throug modules
      for (sheet_name in sheet_names) {
        if (sheet_name %in% wgcna_exclude) {
          
          # Load list of genes for module to explore
          module_genes <- c(read_excel(paste0(output_folder, wgcna_folder, "module_genes.xlsx" ), sheet = sheet_name))[[1]]
          
          # Filter the gene ranking list with the module genes
          gene_rankings_subset <- gene_rankings[!names(gene_rankings) %in% module_genes]
          
          # Name of this analyis
          analysis_name <- paste0(file_path_sans_ext(basename(file)),  "_no_", sheet_name)
          
          # Run GSEA and enrichr on genes 
          result <- run_enrichment(which = c("gsea" ,"panther", "enrichr"), output_dir, analysis_name, 
                                   all_genes=gene_rankings_subset, 
                                   subset=select_genes_for_enrich(gene_rankings_subset, ...), ...)
        }
      }
    }
    
    # Enrichment normal 
    else {
      
      # Messages
      message("Running enrichment on all data given")
      
      # Name of this analysis (will be used to create folder name)
      analysis_name <- file_path_sans_ext(basename(file))
      
      # Run GSEA and enrichr on genes 
      result <- run_enrichment(which = c("panther"), output_dir, analysis_name, 
                               all_genes = gene_rankings, 
                               subset = select_genes_for_enrich(gene_rankings, ...), ...)
    }
  }
  
}

# Main to be run for enrichment analysis
run_enrichment <- function(which = c(""), output_dir, analysis_name, all_genes = c(""), subset = c(""), ...) {
  
  # all_genes -> all the genes (a gene ranking output, so a numeric vector with the scores in which the names are the gene names)
  # subset -> used with ora and enrichr (enrichment analysis in which a threshold needs to be set and only a small number of genes used)
  
  if (!is.character(which)) {
    stop("which argument must be a character vector")
  }
  
  check_packages(c("DOSE", "enrichplot"))

  results <- list()
  # Run GSEA
  if ("gsea" %in% which) {
    results[["gsea"]][[analysis_name]] <- cluster_profiler(all_genes, "BP", "gsea",
                                                           output_dir, analysis_name, ...)
    
    results[["gsea"]][[analysis_name]] <- cluster_profiler(all_genes, "MF", "gsea",
                                                           output_dir, analysis_name, ...)
  }
  # Run ORA
  if ("ora" %in% which) {
    results[["ora"]][[analysis_name]] <- cluster_profiler(subset, "BP", "ora",
                                                          output_dir, analysis_name, universe=names(all_genes), ...)
    
    results[["ora"]][[analysis_name]] <- cluster_profiler(subset, "MF", "ora",
                                                          output_dir, analysis_name, universe=names(all_genes), ...)
  }
  # Run enrichr
  if ("enrichr" %in% which) {
    results[["enrichr"]][[analysis_name]] <- over_expression(subset, "BP", "enrichr",
                                                             output_dir, analysis_name, ...)
    
    results[["enrichr"]][[analysis_name]] <- over_expression(subset, "MF", "enrichr",
                                                             output_dir, analysis_name, ...)
  }
  # Run panther
  if ("panther" %in% which) {
    results[["panther"]][[analysis_name]] <- over_expression(subset, "BP", "panther",
                                                             output_dir, analysis_name, ...)
    
    results[["panther"]][[analysis_name]] <- over_expression(subset, "MF", "panther",
                                                             output_dir, analysis_name, ...)
  }
  return(results)
}

# Function to perfom gsea with clusterProfiler: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
cluster_profiler <- function(gene_rankings, go, type, output_dir, analysis_name, 
                             universe = c(""), minGSSize = 3, maxGSSize = 800, ...) {

  # https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
  # from this guide it seems like the gene sets should be between 15 and 500
  
  # Check types of arguments
  if (!is.numeric(gene_rankings) && !is.character(gene_rankings)) {
    stop("gene_rankings must be a numeric or character vector")
  }
  
  if (!is.character(type) || length(type) != 1 || !(type %in% c("gsea", "ora"))) {
    stop("type argument must be a single character string and one of 'gsea' or 'ora'")
  }
  
  if (!is.character(go) || length(go) != 1 || !(go %in% c("BP", "CC", "MF", "ALL"))) {
    stop("go argument must be a single character string and one of 'BP', 'CC', 'MF', or 'ALL'")
  }
  
  if (!is.character(universe)) {
    stop("universe argument must be a character vector")
  }
  
  # Check useful packages
  check_packages(c("clusterProfiler", "org.Hs.eg.db"))

  # Run GSEA
  if (type == "gsea")  result <- gseGO(geneList=gene_rankings, 
                                   ont = go, 
                                   keyType = "SYMBOL", 
                                   minGSSize = minGSSize, 
                                   maxGSSize = maxGSSize, 
                                   pvalueCutoff = 0.05, 
                                   verbose = TRUE, 
                                   OrgDb = "org.Hs.eg.db")
  # Run over representation analysis
  else if (type == "ora") result <- enrichGO(gene = gene_rankings, 
                                             ont = go, 
                                             universe = universe,
                                             keyType = "SYMBOL", 
                                             minGSSize = minGSSize, 
                                             maxGSSize = maxGSSize, 
                                             pvalueCutoff = 0.05, 
                                             OrgDb = "org.Hs.eg.db")

  else stop("cluster profile enrichment: not valid type parameter")

  # Call function to save the results, in this case it is an enrichResult object
  enrichment_save_results(output_dir, type, paste0(go, "_", analysis_name), result, ...)
  return(result)
}

# Over expression analysis using rbioapi: https://cran.r-project.org/web/packages/rbioapi/index.html
over_expression <- function(genes, go, type, output_dir, analysis_name, 
                            organism = 9606, num_tries = 3, ...) {
  
  # Argument check
  if (!is.vector(genes) || length(genes) == 0) {
    stop("Genes input should be a non-empty vector.")
  }
  
  if (!is.character(type) || length(type) == 0) {
    stop("method should be a string.")
  }
  
  # Messages
  message("Running over-expression analysis with ", type)
  
  # Dependencies
  check_packages(c("rbioapi"))
  
  # Method
  if (type == "enrichr") {
    
    # Enrichr analysis (sometimes it gives problems with connection, like this it tries again multiple times)
    # https://rdrr.io/cran/rbioapi/man/rba_enrichr.html
    if (go == "BP") go_library <- "go_Biological_Process_2017"
    else if (go == "MF") go_library <- "go_Molecular_Function_2017"
    else stop()
    
    for (i in 1:num_tries) {
      result <- try({
        result <- rba_enrichr(
          genes,
          description = NULL,
          gene_set_library = go_library,
          regex_library_name = TRUE,
          organism = "human",
          progress_bar = FALSE
        )[[1]]
        
      })
      if (class(result) != "try-error") {
        message("Success!")
        break
      }
    }
    # Conversion to enrichResult format
    result <- result %>% 
      dplyr::rename(GeneRatio = Overlap, 
                    p.adjust = Adjusted.P.value, 
                    pvalue = P.value,
                    geneID = Genes)  %>%
      tidyr::separate(Term, into = c("Description", "ID"), sep = " \\(")
    
    result$Count <- as.numeric(sub("/.*", "", result$GeneRatio))
    result$ID <- gsub("\\)", "", result$ID)
    result$geneID <- gsub(";", "/", result$geneID)
    result$BgRatio <- result$GeneRatio
    result$qvalue <- result$p.adjust
    result <- result[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")]
    # result <- result[, c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")]
    
    # Conversion to enrichResult format
    result <- new("enrichResult",
                  result = result,  # Your result data.frame
                  pvalueCutoff = 0.05,    # P-value cutoff
                  pAdjustMethod = "BH",  # Method for multiple testing adjustment
                  qvalueCutoff = 0.05,    # Q-value cutoff
                  organism = as.character(organism),     # Organism name
                  ontology = "go",        # Ontology type
                  gene = genes,          # Gene identifier type
                  keytype = "EntrezID",   # Type of gene identifiers
                  universe = "AllGenes",  # Universe of genes
                  gene2Symbol = character(), # Gene symbol mapping
                  geneSets = list(),      # List of gene sets
                  readable = FALSE,       # Readability flag
                  termsim = matrix(),     # Term similarity matrix
                  method = "enrichment", # Method used
                  dr = list()             # Data range
    )
  } else if (type == "mieaa") {
    message("Running over expression analysis with enrichr")
    stop("mieaa to implement")
    # miEAA analysis
    result <- rba_mieaa_enrich(
      genes = genes,
      species = organism,
      dataset = "" 
      # rba_mieaa_cats()
    )
    
    # Conversion to enrichResult format

  } else if (type == "panther") {
    # go = "ANNOT_TYPE_ID_PANTHER_go_SLIM_BP"
    # PANTHER analysis
    # this does not rethurn a dataframe but a results object
    # https://rdrr.io/cran/rbioapi/man/rba_panther_enrich.html
    # https://cran.r-project.org/web/packages/rbioapi/vignettes/rbioapi_panther.html
    # go BP all: "go:000815"
    
    if (go == "BP") go_library <- "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP"
    else if (go == "MF") go_library <- "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF"
    
    for (i in 1:num_tries) {
      result <- try({
        
        result <- rba_panther_enrich(
          genes = genes,
          organism = as.numeric(organism),
          annot_dataset = go_library
        )
        
      })
      if (class(result) != "try-error") {
        message("Success!")
        break
      }
    }
    
    tryCatch({
      result$result$signed_fold_enrichment <- with(result$result, 
                                                   ifelse(plus_minus == "+", 
                                                          fold_enrichment, # Yes
                                                          -fold_enrichment)) # No
    }, error = function(e) {
      warning("An error occurred during ", go, type, analysis_name, ": \n", e)
    })

    
    
  }
  
  else stop("over_expression, not valid method: ", type, " - it can be enrichr, mieaa, panther")
  
  # Call the function to save the results
  enrichment_save_results(output_dir, type, paste0(go, "_", analysis_name), result, ...)
  return(result)
}

# Save the results (excel and calls function for plots)
enrichment_save_results <- function(output_dir, type, analysis_name, result, raw = FALSE, ...) {
  
  # Check types of arguments
  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }
  
  if (!is.character(analysis_name) || length(analysis_name) != 1) {
    stop("analysis_name must be a single character string")
  }
  
  if (!is.logical(raw) || length(raw) != 1) {
    stop("raw must be a single logical value")
  }
  
  # Check useful packages
  check_packages(c("openxlsx", "readxl"))
  
  # Define folder names
  output_subdir <- paste0(output_dir, type, "_", analysis_name, "/") # (to analysis_name the go or enrichment libray name is added in the funtion where it is run)
  filename <-  paste0(type, "_", analysis_name)
  
  # If needed create directory name (unique for each analysis)
  if (!dir.exists(output_subdir)) dir.create(output_subdir)
  
  # Save data
  message("saving results:")
  if (raw) saveRDS(result, file = paste0(output_subdir, "raw_", hash(filename), ".rds"))
  if (type == "panther") write.xlsx(result$result, file = paste0(output_subdir, "results", hash(filename), ".xlsx"))
  else write.xlsx(result@result, file = paste0(output_subdir, "results", hash(filename), ".xlsx"))
  # message
  message(paste0("results saved in: ", output_subdir))
  
  # Plots
  enrichment_plotting(output_subdir, type, filename, result, ...)
  

}

# Plotting fucntion (result can be either of class enrichResult or gseaResult, posso usare anche un generic method)
enrichment_plotting <- function(output_subdir, type, filename, gse, 
                                   ridge_n = 20, qthreshold = 0.05, fc_threshold = 0.5, 
                                   extension_plot = ".png", ...) {
  
  # filename and directory are built in the enrichment_save_results function 
  
  # Check useful packages
  check_packages(c("enrichplot", "DOSE", "svglite"))
  
  ## Function definitions ----
  ridgeplot <- function(x, show_category=ridge_n, fill="NES",
                        core_enrichment = TRUE, label_format = 30, qthreshold_in=force(qthreshold)) {
    
    if (!is(x, "gseaResult"))
      stop("currently only support gseaResult")
    
    ## fill <- match.arg(fill, c("pvalue", "p.adjust", "qvalue"))
    if (fill == "qvalue") {
      fill <- "qvalues"
    }
    
    if (!fill %in% colnames(x@result)) {
      stop("'fill' variable not available ...")
    }
    
    x@result <- x@result[order(x@result$p.adjust),]
    
    
    # Prende i geni contenuti nei pathways (geneSets non e gene list ma sono i pathways con i geni)
    if (!is.na(qthreshold)) show_category <- nrow(x@result[x@result$qvalue < qthreshold_in,])
    n_plots <- show_category
    if (n_plots > 25) n_plots <- 25
    else if (n_plots < 10) n_plots <- 10
    gs2id <- geneInCategory(x)[seq_len(n_plots)]
    if (core_enrichment) {
      gs2id <- geneInCategory(x)[seq_len(n_plots)]
    } else {
      gs2id <- x@geneSets[x$ID[seq_len(n_plots)]]
    }
    
    # Assegna lo score ad ogni gene nei geneSet selezionati, prendendolo da gene list
    gs2val <- lapply(gs2id, function(id) {
      res <- x@geneList[id]
      res <- res[!is.na(res)]
    })
    
    # sostituisce gli ID dei pathways con i nomi leggibili (non ID ma nomi)
    nn <- names(gs2val)
    i <- match(nn, x$ID)
    nn <- x$Description[i]
    
    # Create ordered indexes for plotting
    j <- order(x$p.adjust[i], decreasing=TRUE)
    
    # Save the lengths of the gene sets
    len <- sapply(gs2val, length)
    
    # Creates the correct dataframe to send as input to ridge (2 cols, pathway, score)
    gs2val_df <- data.frame(category = rep(nn, times=len),
                            color = rep(x[i, fill], times=len),
                            value = unlist(gs2val))
    # Parameters of the plot
    colnames(gs2val_df)[2] <- fill
    gs2val_df$category <- factor(gs2val_df$category, levels=nn[j]) 
    
    # Ridge plot con ggplot
    
    ggplot(gs2val_df, aes(x=.data[["value"]], y=.data[["category"]], fill=.data[[fill]])) +
      ggridges::geom_density_ridges() + 
      #scale_fill_continuous(low="blue", high="red", name = fill, limits = c(-4, 4),
      #                      guide=guide_colorbar(reverse=TRUE)) + 
      scale_fill_gradient2(low="blue", mid="white", high="red", name = fill, midpoint=0, limits = c(-2.7, 4),
                           guide = guide_colorbar(reverse = FALSE)) +
      # breaks = seq(-4, 4, by = 1)) 
      scale_y_discrete(labels = function(x) str_wrap(x, width = label_format)) +
      xlab(NULL) + ylab(NULL) +  theme_dose()
  }
  
  # Volcano plot
  # all'interno c'e una variabile dipendente dall'environment di questa funzione (filename) (non posso usare questa funzione altrove)
  volcano_plot <- function(source, filename_in = force(filename), x_name = "NES", y_name = "qvalue", labels = "Description", fc_threshold = 0.5) {
    
    # Check packages
    check_packages(c("EnhancedVolcano"))

    # Check that there are no duplicate columns
    source <- source %>%
      dplyr::group_by(across(all_of(labels))) %>%
      dplyr::filter(n() == 1) %>%
      dplyr::ungroup()

    # Create dataframe to plot
    if (any(class(source) == "gseaResult")) {x_name = "qvalue" ; y_name = "NES"}
    
    source_df <- source[, c(x_name, y_name)]
    rownames(source_df) <- source[[labels]] 
    if (!nrow(source[source[[y_name]] < qthreshold, ])) qthreshold <- source[[y_name]][[ridge_n]]

    # Volcano plot
    EnhancedVolcano(source_df,
                    lab = rownames(source_df),
                    x = x_name,
                    y = y_name,     
                    pCutoff = qthreshold,
                    FCcutoff = fc_threshold) + 
      labs(subtitle = filename_in) + xlab(x_name) + ylab(paste0("log10 ", y_name)) + theme(legend.position = "none")
  }
  
  # Heatmap
  heatplot <- function(x, show_category = 30,
                       label_format = 30, genes = c("")) {
    
    gs2id <- x@geneSets[x$ID]
    
    gs2val <- lapply(gs2id, function(id) {
      res <- x@geneList[id]
      res <- res[!is.na(res)]
    })
    
    # sostituisce gli ID dei pathways con i nomi leggibili (non ID ma nomi)
    nn <- names(gs2val)
    i <- match(nn, x$ID)
    nn <- x$Description[i]
    
    # Save the lengths of the gene sets
    len <- sapply(gs2val, length)
    
    # Creates the correct dataframe to send as input to ridge (2 cols, pathway, score)
    gs2val_df <- data.frame(geneset = rep(nn, times=len),
                            foldChange = unlist(gs2val),
                            gene = sub(".*\\.", "",names(unlist(gs2val))))
    
    gs2val_df <- gs2val_df[gs2val_df$gene %in% genes,]
    
    # gs2val_df$geneset <- factor(gs2val_df$geneset, levels=nn[j])
    
    p <- ggplot(gs2val_df, aes_string(x="gene", y="geneset")) +
      geom_tile(aes_string(fill = "foldChange"), color = "white") +
      scale_fill_continuous(low="blue", high="red", name = "fold change") +
      ## scale_fill_gradientn(name = "fold change", colors = palette)
      xlab(NULL) + ylab(NULL) + theme_minimal() +
      scale_y_discrete(labels = function(x) str_wrap(x, width = label_format))  +
      theme(panel.grid.major = element_blank(),
            axis.text.x=element_text(angle = 60, hjust = 1))
  }
  
  # Upset plot
  upsetplot <- function(x, show_category = 30,
                        label_format = 30, genes = c("")) {
    check_packages(c("ggupset"))
    
    gs2id <- x@geneSets[x$ID]
    
    gs2val <- lapply(gs2id, function(id) {
      res <- x@geneList[id]
      res <- res[!is.na(res)]
    })
    
    # sostituisce gli ID dei pathways con i nomi leggibili (non ID ma nomi)
    nn <- names(gs2val)
    i <- match(nn, x$ID)
    nn <- x$Description[i]
    
    # Save the lengths of the gene sets
    len <- sapply(gs2val, length)
    
    # Creates the correct dataframe to send as input to ridge (2 cols, pathway, score)
    gs2val_df <- data.frame(geneset = rep(nn, times=len),
                            foldChange = unlist(gs2val),
                            gene = sub(".*\\.", "",names(unlist(gs2val))))
    
    gs2val_df <- gs2val_df[gs2val_df$gene %in% genes,]
    # gs2val_df$geneset <- factor(gs2val_df$geneset)
    # gs2val_df$geneset <- factor(gs2val_df$geneset, levels=nn[j])
    gs2val_df$gene <- NULL
    gs2val_df$geneset <- factor(gs2val_df$geneset)
    gs2val_df$geneset <- as.list(gs2val_df$geneset)
    
    rownames(gs2val_df) <- NULL
    
    ggplot(gs2val_df, aes_string(x = "geneset", y = "foldChange")) +
      geom_boxplot() +
      geom_jitter(width = .2, alpha = .6) +
      theme_dose(font.size = 12) +
      xlab(NULL) + ylab(NULL) + scale_x_upset(order_by = "degree")
  }
  
  # Horizontal barplot
  horizontal_barplot <- function(df, columns, pvalue = "pValue", color_limits = FALSE) {
    # Check if the columns argument is valid
    if (length(columns) != 3) {
      stop("Please provide exactly two column names.")
    }

    col_length <- columns[1]
    col_color <- columns[2]
    names <- columns[3]
    
    # Check if the specified columns exist in the dataframe
    if (!(col_length %in% colnames(df)) || !(col_color %in% colnames(df)) || !(names %in% colnames(df)) || !(pvalue %in% colnames(df))) {
      stop("One or both of the specified columns do not exist in the dataframe.")
    }
    
    # select only the significant rows, delete all na rows
    df <- na.omit(df[df[[pvalue]] < 0.05,][order(df[[pvalue]]),])
    
    # Set color limits to correct values if needed
    if (color_limits == FALSE) {
      max_n <- max(abs(df[col_color]))
      
      color_limits <- c(-max_n, max_n)
    }
    # Create a ggplot object for the horizontal barplot

    p <- ggplot(df, aes(y = .data[[names]], x = .data[[col_length]], fill = .data[[col_color]])) +
      geom_bar(stat = "identity") +
      scale_fill_gradient2(low = "blue", mid ="#99e599", high = "red", limits = color_limits) +
      labs(x = col_length, y = "") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10))
  }
  

  ## plotting ----

  if (class(gse) == "gseaResult") {
    tryCatch({
      save_plot(dotplot(gse, show_category = 10, split = ".sign") + facet_grid(.~.sign),
                paste0(output_subdir, "dotplot_", hash(filename), extension_plot), x = 10, y = 12)
      #save_plot(heatplot(gse) + ggtitle("heatplot for GSEA"),
      #          paste0(output_dir, "heatplot_", result_name, extension_plot), x = 5, y = 3)
      
      save_plot(ridgeplot(gse) + labs(x = "log2FC distribution"),
                paste0(output_subdir, "ridgeplot_", hash(filename), extension_plot), x = 10, y = 12)
      save_plot(volcano_plot(as.data.frame(gse@result)),
                paste0(output_subdir, "volcano_", hash(filename), extension_plot), x = 10, y = 12)
    }, error = function(e) {
      message("The plots could not be created, probably the enrichment didn't find anything significant, error: \n", e)
    })
  }
  
  if (type == "ora" || class(gse) == "gseaResult") {
    tryCatch({
      
      gse <- pairwise_termsim(gse) # https://rdrr.io/bioc/enrichplot/man/pairwise_termsim.html
      message("Following warning message ok to ignore: https://github.com/YuLab-SMU/ggtree/issues/577")
      save_plot(treeplot(gse),
                paste0(output_subdir,  "tree_plot_", hash(filename), extension_plot), x = 10, y = 10)
      save_plot(emapplot(gse),
                paste0(output_subdir,  "enrichment_map_", hash(filename), extension_plot), x = 10, y = 10)
    }, error = function(e) {
      message("The plots could not be created, probably the enrichment didn't find anything significant, error: \n", e)
    })
  }
  
  if (type == "ora") {
    tryCatch({
      save_plot(horizontal_barplot(as.data.frame(gse@result), c("pvalue", "Count", "Description")),
                paste0(output_subdir, "horizontal_barplot_", hash(filename), extension_plot), x = 10, y = 12)
    }, error = function(e) {
      message("An error occurred (horizontal_barplot):  \n", e)
    })
    tryCatch({
      save_plot(volcano_plot(as.data.frame(gse@result), y = "pvalue", x = "Count", labels = "Description", fc_threshold=0),
                paste0(output_subdir, "volcano_", hash(filename), extension_plot), x = 10, y = 12)
      
    }, error = function(e) {
      message("An error occurred (volcano):  \n", e)
    })
  }
  
  if (type == "panther") {
    tryCatch({
      save_plot(horizontal_barplot(as.data.frame(gse$result), c("pValue", "signed_fold_enrichment", "term.label")),
                paste0(output_subdir, "horizontal_barplot_", hash(filename), extension_plot), x = 10, y = 12)
    }, error = function(e) {
      message("An error occurred (horizontal_barplot):  \n", e)
    })
    tryCatch({
      save_plot(volcano_plot(as.data.frame(gse$result), y = "pValue", x = "signed_fold_enrichment", labels = "term.label", fc_threshold=0),
                paste0(output_subdir, "volcano_", hash(filename), extension_plot), x = 10, y = 12)
      
    }, error = function(e) {
      message("An error occurred (volcano):  \n", e)
    })
    tryCatch({
      save_plot(volcano_plot(as.data.frame(gse$result), y = "pValue", x = "signed_fold_enrichment", labels = "term.label", fc_threshold=0),
                paste0(output_subdir, "volcano_", hash(filename), extension_plot), x = 10, y = 12)
      
    }, error = function(e) {
      message("An error occurred (volcano):  \n", e)
    })
  }
  
  if (type == "enrichr") {
    tryCatch({
      save_plot(horizontal_barplot(as.data.frame(gse@result), c("Count", "pvalue", "Description"), pvalue="pvalue"),
                paste0(output_subdir, "horizontal_barplot_", hash(filename), extension_plot), x = 10, y = 12)
    }, error = function(e) {
      message("An error occurred (horizontal_barplot):  \n", e)
    })
    tryCatch({
      save_plot(horizontal_barplot(as.data.frame(gse@result), c("pvalue", "Count", "Description"), pvalue="pvalue"),
                paste0(output_subdir, "horizontal_barplot_count", hash(filename), extension_plot), x = 10, y = 12)
    }, error = function(e) {
      message("An error occurred (horizontal_barplot):  \n", e)
    })
    tryCatch({
      save_plot(volcano_plot(as.data.frame(gse@result), y = "pvalue", x = "Count", labels = "Description", fc_threshold=0),
                paste0(output_subdir, "volcano_", hash(filename), extension_plot), x = 10, y = 12)
    }, error = function(e) {
      message("An error occurred (volcano):  \n", e)
    })
  }

}

# Helper function: Prepare the seurat object for the enrichment analysis
prepare_data <- function(seurat_object) {
  
  # Delete genes that have a total expression lower than 15 (from where is the source)
  total_expression <- rowSums(seurat_object@assays[[DefaultAssay(seurat_object)]]$counts)
  genes_to_remove <- names(total_expression[total_expression < 15])
  seurat_object <- subset(seurat_object, features = rownames(seurat_object)[!(rownames(seurat_object) %in% genes_to_remove)])
  
  return(seurat_object)
}

# Helper function: to prepare the ordered gene list
prepare_genes <- function(excel_file, count_threshold, scoring = "log2FC", ...) {
  
  if (!is.character(scoring) || length(scoring) != 1 || !(scoring %in% c("log2FC", "spvalue", "paper"))) 
    stop("scoring argument must be a single character string and one of 'log2FC', 'spvalue', or 'paper'")

  
  library("readxl")

  # Read the Excel file and filter out low quality data
  data <- read_excel(excel_file)
  data <- data[data$pct.1 > count_threshold & data$pct.2 > count_threshold, ]
  
  # Extract the data
  genes <- data$gene
  if (scoring == "log2FC") scores <- data$avg_log2FC
  else if (scoring == "spvalue") scores <- sign(data$avg_log2FC) * (-log10(data$p_val))
  else if (scoring == "paper") scores <- data$avg_log2FC * (-log10(data$p_val))
  else stop("non valid scoring method selection")
  
  names(scores) <- genes
  
  # Order list and set Inf and -Inf to 1000 and -1000
  scores <- sort(scores, decreasing = TRUE)
  scores[scores == Inf] <- NA
  scores[scores == -Inf] <- NA
  scores[is.na(scores)] <- max(scores, na.rm = TRUE)
  # scores[scores == -Inf] <- min(scores)*10
  
  return(scores)
}

# Helper function: to prepare the vector of genes to enrich
select_genes_for_enrich <- function(gene_rankings, n_gene_enrich = 400, ...) {
  
  if (!is.numeric(n_gene_enrich) || length(n_gene_enrich) != 1) {
    stop("n_gene_enrich argument must be a single numeric value")
  }
  
  return(names(gene_rankings[order(abs(gene_rankings), decreasing = TRUE)])
         [seq_len(ifelse(n_gene_enrich < length(gene_rankings), # is shorter than desired lenght? 
                   n_gene_enrich, # yes -> all genes in list
                   length(gene_rankings)))]) # no -> take first n
  
}

# Helper function: to select which modules to enrich based on the previously fitted linear model pvalue
.select_mdoules_to_enrich <- function(filename) {
  data <- openxlsx::read.xlsx(filename, rowNames = TRUE)
  filtered_data <- rownames(data[data$Pr...t.. < 0.05, ])
}
