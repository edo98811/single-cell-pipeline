
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
            cols = 1:ncol(data)+1, 
            rows = 1:nrow(data)+1, 
            type = "colorScale", 
            rule = c(1, 0, -1),
            style = c("#78b7ff", "#ffffff", "#F8696B"),
            gradient = TRUE
        )
    else if (mode == "pvalue") {
        openxlsx::conditionalFormatting(
            wb, 
            sheet = sheet_name, 
            cols = 1:ncol(data)+1, 
            rows = 1:nrow(data)+1, 
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
    markers_analysisPD = "microglia_control_vs_pd_nogenetic",
    markers_analysisGPD = "microglia_control_vs_geneticpd", 
    cluster = FALSE, ...) {

    library("openxlsx")
    
    if (isFALSE(cluster)) {
        message(paste0("loading... ", output_folder, "markers_", markers_analysisPD,
            "/expressed_markers_all_", markers_analysisPD, ".xlxs"))
        markers_tablePD <- read.xlsx(paste0(output_folder, "markers_", markers_analysisPD,
            "/expressed_markers_all_", markers_analysisPD, ".xlsx"))

        message(paste0("loading... ", output_folder, "markers_", markers_analysisGPD,
            "/expressed_markers_all_", markers_analysisGPD, ".xlxs"))
        markers_tableGPD <- read.xlsx(paste0(output_folder, "markers_", markers_analysisGPD,
            "/expressed_markers_all_", markers_analysisGPD, ".xlsx"))
    } else {
        message(paste0("loading... ", output_folder, "markers_", markers_analysisPD,
            "/expressed_markers_", cluster, "_", markers_analysisPD, ".xlxs"))
        markers_tablePD <- read.xlsx(paste0(output_folder, "markers_", markers_analysisPD,
            "/expressed_markers_", cluster, "_", markers_analysisPD, ".xlsx"))

        message(paste0("loading... ", output_folder, "markers_", markers_analysisGPD,
            "/expressed_markers_", cluster, "_", markers_analysisGPD, ".xlxs"))
        markers_tableGPD <- read.xlsx(paste0(output_folder, "markers_", markers_analysisGPD,
            "/expressed_markers_", cluster, "_", markers_analysisGPD, ".xlsx"))
    }
    
    return(list(GPD = tibble::column_to_rownames(as.data.frame(markers_tableGPD), var = "gene"),
        PD = tibble::column_to_rownames(as.data.frame(markers_tablePD), var = "gene")))
}

merge_dataframes <- function(data1, data2, data3){

    plot_df <- tibble::column_to_rownames(
        merge(
            tibble::column_to_rownames(merge(data1, data2, by = "row.names"), var = "Row.names"), 
            data3, 
            by = "row.names"
        ), 
        var = "Row.names")

}

load_mri_traits <- function(type = "zscore") {

    # Compute membershio measure
    source("scripts/new/mri_analysis.r", local = TRUE)

    if (type == "zscore") {
        traits <- zscore_for_wgcna()
    } else if (type == "mri") {
        traits <- traits <- mri_for_wgcna()
    } else stop("in wgcna plot load significance wrong load type")

    return(traits)
}

load_significance <- function (norm_counts, type= "zscore") {

    traits <- load_mri_traits(type = type)

    gene_significance <- cor(norm_counts, traits, use = "p")
}

load_log2fc_df <- function(cluster, deg_to_use = "PD", ...) {

    plot_df <- load_log2fc(cluster, ...)[[deg_to_use]]

    plot_df$abs_avg_log2FC <- abs(plot_df$avg_log2FC)

    return (plot_df)
}

load_module_membership <- function (bwnet, norm_counts) {

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

    wb <- write_on_excel("correllation", as.data.frame(res$cor), mode = "colorscale")
    wb <- write_on_excel("p.value", as.data.frame(res$p), wb = wb)
    wb <- write_on_excel("t.statistic", as.data.frame(res$t), wb = wb)
    openxlsx::saveWorkbook(wb, excel_filename, overwrite = TRUE)

    message("heatmap saved in: ", output_dir, excel_filename)
}

# In case it is needed
# all'interno c'e una variabile dipendente dall'environment di questa funzione (filename) (non posso usare questa funzione altrove)
volcano_plot <- function(source, filename_in = force(filename), x_name = "NES", y_name="qvalue", 
                        labels = "Description", FCthreshold=0.5, qthreshold=0.05) {
    
    # Check packages
    check_packages(c("EnhancedVolcano"))
    
    # Create dataframe to plot
    if(class(source) == "gseaResult") {x_name = "qvalue"; y_name="NES"}
    
    source_df <- source[,c(x_name, y_name)]
    rownames(source_df) <- source[[labels]]
    if (!nrow(source[source[[y_name]] < qthreshold,])) qthreshold <- source[[y_name]][[ridge_n]]

    # Volcano plot
    EnhancedVolcano(source_df,
                    lab = rownames(source_df),
                    x = x_name,
                    y = y_name,     
                    pCutoff = qthreshold,
                    FCcutoff = FCthreshold) + 
    labs(subtitle = filename_in) + xlab(x_name) + ylab(paste0("log10 ",y_name)) + theme(legend.position="none")
}

# to include again in case
plot_hub_genes <- function(df) {

FCcutoff <- hub_genes_threshold[[2]]
MMcutoff <- hub_genes_threshold[[1]]

save_plot(ggplot(df, aes(x = .data[[module_membership_column]],
y = .data[["avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
    x) + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5,
color = "black", label.x = "right", label.y = "top") + theme_minimal() +
geom_hline(yintercept = c(-FCcutoff, FCcutoff), colour = "black") +
geom_vline(xintercept = c(-MMcutoff, MMcutoff), colour = "black") +
geom_text_repel(data = subset(df, abs(df[["avg_log2FC"]]) > FCcutoff &
    abs(df[[module_membership_column]]) > MMcutoff), aes(label = row.names(subset(df,
    abs(df[["avg_log2FC"]]) > FCcutoff & abs(df[[module_membership_column]]) >
    MMcutoff))), xlim = c(NA, NA), ylim = c(NA, NA)), paste0(output_dir,
deg_to_use, "_scatter_hub/", module, "_scatter_plot", extension_plot))

if (module %in% module_of_interest) {
hub_genes_list <- row.names(subset(submod, abs(df[["avg_log2FC"]]) >
    FCcutoff & abs(df[[module_membership_column]]) > MMcutoff))

# feature_plots(seurat_object_microglia, hub_genes_list,
# name = paste0('WGCNA_', name, '/feature_plot_hubgene_',
# deg_to_use, '_', module, '/'), reduction_name =
# reduction_name)

write.xlsx(data.frame(hub_genes_list), paste0(output_dir, deg_to_use,
    module, "_hubgene.xlsx"))

many_plots(seurat_object_microglia, which = c("dotplot", "heatmap"),
    cluster_column = "microglia_clusters", name = paste0("WGCNA_modules/",
    module, deg_to_use), markers = hub_genes_list)
}
}
