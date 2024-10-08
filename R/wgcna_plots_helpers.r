# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


write_on_excel <- function(sheet_name, data, wb = NULL, mode = "", small_cells = FALSE) {

    library(openxlsx)
    
    # Create a new workbook and add a worksheet
    if (is.null(wb)) wb <- openxlsx::createWorkbook()
    
    openxlsx::addWorksheet(wb, sheet_name)

    # Write the data to the worksheet
    openxlsx::writeDataTable(wb, sheet_name, data, rowNames = TRUE)


    # Define the style for conditional formatting
    if (mode == "colorscale")
        openxlsx::conditionalFormatting(
            wb, 
            sheet = sheet_name, 
            cols = seq_len(ncol(data) + 1), 
            rows = seq_len(nrow(data) + 1), 
            type = "colorScale", 
            rule = c(1, 0, -1),
            style = c("#78b7ff", "#ffffff", "#F8696B"),
            gradient = TRUE
        )
    else if (mode == "pvalue") 
        openxlsx::conditionalFormatting(
            wb, 
            sheet = sheet_name, 
            cols = seq_len(ncol(data) + 1), 
            rows = seq_len(nrow(data) + 1), 
            type = "expression", 
            rule = "<=0.05",
            style = createStyle(bgFill = "yellow"),
            gradient = TRUE
        )
    if (small_cells) {
        openxlsx::setColWidths(
            wb,
            sheet = sheet_name,
            cols = seq(2, ncol(data) + 1),
            widths = 4.78,
            ignoreMergedCells = FALSE
        )
        setRowHeights(
            wb,
            sheet = sheet_name,
            rows = seq(2, nrow(data) + 1), 
            heights = 30,
            fontsize = NULL,
            factor = 1,
            base_height = 15,
            extra_height = 12,
            wrap = TRUE
        )
        grid <- expand.grid(rows = seq(2, nrow(data) + 1), cols = seq(2, ncol(data) + 1))
        center_style <- createStyle(halign = "center", valign = "center", numFmt = "NUMBER")
        addStyle(
            wb, 
            sheet = sheet_name, 
            style = center_style,             
            cols = grid$cols, 
            rows = grid$rows
        )
    }
    return(wb)
}

# Per caricare la lista di markers dalla cartella, dovrebbe essere ok caricare
# anche senza fare in modo che funzioni per i clusters, magari pensare ad un
# mod di aggiungerlo, ma non priorita
load_log2fc <- function(
    markers_analysis = "",
    cluster = FALSE, ...) {
    
    tryCatch({
        if (isFALSE(cluster)) {
            message(paste0("loading... ", output_folder, "markers_", markers_analysis,
                "/expressed_markers_all_", markers_analysis, ".xlxs"))
            markers_table <- openxlsx::read.xlsx(paste0(output_folder, "markers_", markers_analysis,
                "/expressed_markers_all_", markers_analysis, ".xlsx"))
        } else {
            message(paste0("loading... ", output_folder, "markers_", markers_analysis,
                "/expressed_markers_", cluster, "_", markers_analysis, ".xlxs"))
            markers_table <- openxlsx::read.xlsx(paste0(output_folder, "markers_", markers_analysis,
                "/expressed_markers_", cluster, "_", markers_analysis, ".xlsx"))
        }
    }, error = function(e) {
      stop("Could not load the deg results to compute histogram in wgcna, probably the table does not exist?: \n", e)
    })
    return(tibble::column_to_rownames(as.data.frame(markers_table), var = "gene"))
}

merge_dataframes <- function(data1, data2, data3) {

    plot_df <- tibble::column_to_rownames(
        merge(
            tibble::column_to_rownames(merge(data1, data2, by = "row.names"), var = "Row.names"), 
            data3, 
            by = "row.names"
        ), 
        var = "Row.names")

}

load_mri_traits <- function(wgcna_subjects = list(), regions = c(), type = "zscore", data_source_mri = "aparc", ...) {

    # Check if either of the variables is empty and raise a warning if true
    if (length(wgcna_subjects) == 0) {
        stop("The wgcna_subjects list is empty, cannot load mri traits without a mapping to sigle cels subjects.")
    }

    if (length(regions) == 0) {
        stop("The regions vector is empty, need to be defined to load mri data.")
    }
    
    # Compute membershio measure
    source("R/mri_wgcna.r", local = TRUE)

    if (type == "zscore") {
        traits <- zscore_for_wgcna(data_source_mri, regions, wgcna_subjects)
    } else if (type == "mri") {
        traits <- mri_for_wgcna(data_source_mri, regions, wgcna_subjects)
    } else stop("in wgcna plot load significance wrong load type")

    if (nrow(traits) == 0) stop ("mri features were not correctly loaded")

    return(traits)
}

load_significance <- function(norm_counts, type = "zscore", ...) {
    traits <- load_mri_traits(type = type, ...)
    tryCatch({
        gene_significance <- cor(norm_counts, traits, use = "p")
        return(gene_significance)
    }, 
    error = function(e) {
        if (grepl("incompatible dimensions", e$message)) {
            stop("Error incompatible dimensions: Maybe your Seurat object contains different samples as the ones in the traits dataset?", "\n",
                "subjects in seurat object: ", rownames("norm_counts"), "\n",
                "subjects in seurat object: ", rownames("norm_counts"))
        } else {
            stop(e)
        }
    })
}

load_log2fc_df <- function(cluster, ...) {

    plot_df <- load_log2fc(cluster, ...)
    plot_df$abs_avg_log2FC <- abs(plot_df$avg_log2FC)

    return(plot_df)
}

load_module_membership <- function(bwnet, norm_counts) {

    module_membership_measure <- cor(norm_counts, bwnet$MEs, use = "p") 
}

make_heatmap <- function(bwnet, excel_filename, traits) {

    heatmap_data <- bwnet$MEs %>%
        merge(traits, by = "row.names") %>%
        tibble::column_to_rownames(var = "Row.names")

    # save correlation results
    res <- WGCNA::corAndPvalue(
        y = heatmap_data[, (length(heatmap_data) - ncol(traits) + 1):length(heatmap_data)],
        x = heatmap_data[, 1:(length(heatmap_data) - ncol(traits))])
        
    wb <- write_on_excel("correlation", as.data.frame(res$cor), mode = "colorscale", small_cells = TRUE)
    wb <- write_on_excel("p.value", as.data.frame(res$p), wb = wb, mode = "pvalue")
    wb <- write_on_excel("t.statistic", as.data.frame(res$t), wb = wb)
    openxlsx::saveWorkbook(wb, excel_filename, overwrite = TRUE)

    message("heatmap saved in: ", excel_filename)
}

# In case it is needed
# all'interno c'e una variabile dipendente dall'environment di questa funzione (filename) (non posso usare questa funzione altrove)
volcano_plot <- function(source, filename_in = force(filename), x_name = "NES", y_name = "qvalue", 
                        labels = "Description", fc_threshold = 0.5, qthreshold = 0.05) {
    
    # Check packages
    check_packages(c("EnhancedVolcano"))
    
    # Create dataframe to plot
    if (class(source) == "gseaResult") {
        x_name <- "qvalue" 
        y_name <- "NES"
        }
    
    source_df <- source[, c(x_name, y_name)]
    rownames(source_df) <- source[[labels]]
    if (!nrow(source[source[[y_name]] < qthreshold, ])) qthreshold <- source[[y_name]][[ridge_n]]

    # Volcano plot
    EnhancedVolcano(source_df,
                    lab = rownames(source_df),
                    x = x_name,
                    y = y_name,     
                    pCutoff = qthreshold,
                    fc_cutoff = fc_threshold) + 
    labs(subtitle = filename_in) + xlab(x_name) + ylab(paste0("log10 ", y_name)) + theme(legend.position = "none")
}
