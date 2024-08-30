# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


write_on_excel <- function(sheet_name, data, wb = NULL, mode = "pvalue") {

    library(openxlsx)
    
    # Create a new workbook and add a worksheet
    if (is.null(wb)) wb <- openxlsx::createWorkbook()
    # detach("package:openxlsx", unload=TRUE)
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
    else if (mode == "pvalue") {
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
    }
        
    # Save the workbook
    return(wb)
}

# Per caricare la lista di markers dalla cartella, dovrebbe essere ok caricare
# anche senza fare in modo che funzioni per i clusters, magari pensare ad un
# mod di aggiungerlo, ma non priorita
load_log2fc <- function(
    markers_analysis_pd = "microglia_control_vs_pd_nogenetic",
    markers_analysis_gpd = "microglia_control_vs_geneticpd", 
    cluster = FALSE, ...) {

    library("openxlsx")
    
    if (isFALSE(cluster)) {
        message(paste0("loading... ", output_folder, "markers_", markers_analysis_pd,
            "/expressed_markers_all_", markers_analysis_pd, ".xlxs"))
        markers_table_pd <- read.xlsx(paste0(output_folder, "markers_", markers_analysis_pd,
            "/expressed_markers_all_", markers_analysis_pd, ".xlsx"))

        message(paste0("loading... ", output_folder, "markers_", markers_analysis_gpd,
            "/expressed_markers_all_", markers_analysis_gpd, ".xlxs"))
        markers_table_gpd <- read.xlsx(paste0(output_folder, "markers_", markers_analysis_gpd,
            "/expressed_markers_all_", markers_analysis_gpd, ".xlsx"))
    } else {
        message(paste0("loading... ", output_folder, "markers_", markers_analysis_pd,
            "/expressed_markers_", cluster, "_", markers_analysis_pd, ".xlxs"))
        markers_table_pd <- read.xlsx(paste0(output_folder, "markers_", markers_analysis_pd,
            "/expressed_markers_", cluster, "_", markers_analysis_pd, ".xlsx"))

        message(paste0("loading... ", output_folder, "markers_", markers_analysis_gpd,
            "/expressed_markers_", cluster, "_", markers_analysis_gpd, ".xlxs"))
        markers_table_gpd <- read.xlsx(paste0(output_folder, "markers_", markers_analysis_gpd,
            "/expressed_markers_", cluster, "_", markers_analysis_gpd, ".xlsx"))
    }
    
    return(list(
        GPD = tibble::column_to_rownames(as.data.frame(markers_table_gpd), var = "gene"),
        PD = tibble::column_to_rownames(as.data.frame(markers_table_pd), var = "gene")))
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
    warning("The wgcna_subjects list is empty.")
    }

    if (length(regions) == 0) {
    warning("The regions vector is empty.")
    }
    
    # Compute membershio measure
    source("scripts/mri_wgcna.r", local = TRUE)

    if (type == "zscore") {
        traits <- zscore_for_wgcna(data_source_mri, regions, wgcna_subjects)
    } else if (type == "mri") {
        traits <- mri_for_wgcna(data_source_mri, regions, wgcna_subjects)
    } else stop("in wgcna plot load significance wrong load type")

    return(traits)
}

load_significance <- function(norm_counts, type = "zscore", ...) {

    traits <- load_mri_traits(type = type, ...)

    gene_significance <- cor(norm_counts, traits, use = "p")
}

load_log2fc_df <- function(cluster, deg_to_use = "PD", ...) {

    plot_df <- load_log2fc(cluster, ...)[[deg_to_use]]

    plot_df$abs_avg_log2FC <- abs(plot_df$avg_log2FC)

    return(plot_df)
}

load_module_membership <- function(bwnet, norm_counts) {

    module_membership_measure <- cor(norm_counts, bwnet$MEs, use = "p")  # ahh ho invertito e quindo ho una direzione divers
}

make_heatmap <- function(bwnet, excel_filename, traits) {
    
    library(tibble)
    library(openxlsx)

    heatmap_data <- bwnet$MEs %>%
        merge(traits, by = "row.names") %>%
        tibble::column_to_rownames(var = "Row.names")

    # save correlation results
    res <- WGCNA::corAndPvalue(
        y = heatmap_data[, (length(heatmap_data) - ncol(traits) + 1):length(heatmap_data)],
        x = heatmap_data[, 1:(length(heatmap_data) - ncol(traits))], method = "kendall")

    wb <- write_on_excel("correlation", as.data.frame(res$cor), mode = "colorscale")
    wb <- write_on_excel("p.value", as.data.frame(res$p), wb = wb)
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
