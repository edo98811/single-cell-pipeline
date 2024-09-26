# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


source("R/wgcna_plots_helpers.r", local = TRUE)

plot_functions <- function(bwnet, norm_counts, column_data, output_dir, extension_plot = ".png",
reduction_name = "umap_microglia_harmony_reduction", module_of_interest = c(""),
hub_genes_threshold = c(0.3, 1), which = c(""), cluster = FALSE, ...) {
    
    if (!is.character(which) || !is.vector(which)) stop("which argument must be a character vector")

    
    for (function_name in which) {
        switch(function_name,
            "heatmap_pathology" = try(heatmap_group(bwnet, norm_counts, column_data)),
            "TOM" = try(tom(bwnet, norm_counts, column_data)),
            "dendro" = try(dendro(bwnet, norm_counts, column_data, extension_plot)),
            "heatmap_mri" = try(heatmap_mri(bwnet, norm_counts, column_data)),
            "heatmap_zscore" = try(heatmap_zscore(bwnet, norm_counts, column_data)),
            "violin_plots" = try(violin_plots(bwnet, norm_counts, column_data)),
            "histogram_plot" = try(histogram_plot(bwnet, norm_counts, column_data, cluster, extension_plot, ...)),
            "histogram_plot_significance" = try(histogram_plot_significance(bwnet, norm_counts, column_data, extension_plot, ...)),
            "significance_membership_scatter" = try(significance_membership_scatter(bwnet, norm_counts, column_data, hub_genes_threshold, cluster, extension_plot, ...)),
            "significance_log2fc_scatter" = try(significance_log2fc_scatter(bwnet, norm_counts, column_data, cluster, extension_plot, ...)),
            "correlation_avglog2fc_scatter" = try(correlation_avglog2fc_scatter(bwnet, norm_counts, column_data, cluster, extension_plot, ...)),
            "corr_matrix" = try(corr_matrix(bwnet, norm_counts, column_data)),
            "significance_membership_model" = try(significance_membership_model(bwnet, norm_counts, column_data, ...)),
            warning(function_name, ": Invalid function name"))
    }
}

# Define function for heatmap
heatmap_group <- function(bwnet, norm_counts, column_data) {
    message("heatmap")
    traits_binary <- data.frame(lapply(unique(column_data$pathology), function(x) {
        as.numeric(column_data$pathology == x)
    }))

    colnames(traits_binary) <- unique(column_data$pathology)
    row.names(traits_binary) <- row.names(column_data)

    make_heatmap(bwnet, paste0(output_dir, "trait_correlation_heatmap.xlsx"), traits_binary)

}

# Define function for heatmapMRI
heatmap_mri <- function(bwnet, norm_counts, column_data) {
    message("heatmap mri")
    traits <- load_mri_traits(type = "mri", ...)
    make_heatmap(bwnet, paste0(output_dir, "trait_correlation_heatmap_mri.xlsx"), traits)
}

# Define function for heatmapMRI zscore
heatmap_zscore <- function(bwnet, norm_counts, column_data) {
    message("heatmap zscore")
    traits <- load_mri_traits(type = "zscore", ...)
    make_heatmap(bwnet, paste0(output_dir, "trait_correlation_heatmap_zscore.xlsx"), traits)

}

# Define function for TOM
tom <- function(bwnet, norm_counts, column_data) {
    message("TOM")

    adjacency_matrix <- WGCNA::adjacency(norm_counts, power = soft_power, type = "signed")
    tom <- 1 - WGCNA::TOMsimilarity(adjacency_matrix)
    saveRDS(TOM, paste0(output_dir, "TOM.rds"))
    save_plot(TOMplot(TOM, bwnet$dendrograms[[1]]))
    # 6B. Intramodular analysis: Identifying driver genes ---------------

    # Calculating the module dissimilarity eigengenes
    mediss <- 1 - cor(MEs)

    # Clustering the eigengenes modules
    metree <- hclust(as.dist(mediss), method = "average")
    medissthres <- 0.25
    # Plotting the result sizeGrWindow(7, 6)
    # save_plot(plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "") + 
    #     abline(h = MEDissThres, col = "red"))


}

# Define function for dendro
dendro <- function(bwnet, norm_counts, column_data, extension_plot) {
    message("dendro")

    png(file = paste0(output_dir, "modules_dendrogram", extension_plot))
    plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
        c("unmerged", "merged"), dendroLabels = FALSE, addGuide = TRUE, hang = 0.03,
        guideHang = 0.05)
    dev.off()
    message("plot saved in: ", output_dir, "modules_dendrogram", extension_plot)
}

# Define function for violin_plots
violin_plots <- function(bwnet, norm_counts, column_data, extension_plot) {
    message("violin_plots")

    n_samples <- nrow(norm_counts)

    traits_binary <- data.frame(lapply(unique(column_data$pathology), function(x) {
        as.numeric(column_data$pathology == x)
    })) 
    colnames(traits_binary) <- unique(column_data$pathology)
    row.names(traits_binary) <- row.names(column_data)

    module_gene_mapping <- as.data.frame(bwnet$colors)
    for (color in unique(module_gene_mapping[[1]])) {

        module_genes <- module_gene_mapping %>%
            filter(bwnet$colors == color) %>%
            rownames()

        norm_counts_subset <- norm_counts[, colnames(norm_counts) %in% module_genes]

        # Calculate the gene significance and associated p-values
        # controllare qui come e stato fatto il testo eh (significance)
        gene_signf_corr <- cor(norm_counts_subset, traits_binary, use = "p")
        gene_signf_corr_pvals <- corPvalueStudent(gene_signf_corr, n_samples)

        # For each module creates te violin plots for the 10 genes with the highest correlation to the traits?
        # TODO: check + source
        for (column in colnames(gene_signf_corr_pvals)) {

            genes <- gene_signf_corr_pvals %>%
                as.data.frame() %>%
                arrange(column) %>%
                head(10) %>%
                row.names()
            
            # Calls the violin plot function from seurat_utiles
            violin_plot(seurat_object, genes, name = paste0("WGCNA/test_2", color,
                column), markers_analysis_pd = "microglia_control_vs_pd_nogenetic_all",
                markers_analysis_gpd = "microglia_control_vs_geneticpd_all", extension_plot = extension_plot)

        }

        # Using the gene significance you can identify genes that have a
        # high significance for trait of interest Using the module
        # membership measures you can identify genes with high module
        # membership in interesting modules.
    }
}

# Define function for histogram_plot
histogram_plot <- function(bwnet, norm_counts, column_data, cluster, extension_plot, ...) {
    message("histogram_plot")

    deg_results <- load_log2fc_df(cluster, ...)

    plot_df <- merge(data.frame(bwnet$colors), deg_results,
        by = "row.names", all = FALSE)

    plot_df$abs_avg_log2FC <- abs(plot_df$avg_log2FC)
    plot_df$bwnet.colors <- factor(plot_df$bwnet.colors)

    # Plot (with correct colors)
    color_mapping <- setNames(c(levels(plot_df$bwnet.colors)), levels(plot_df$bwnet.colors))
    save_plot(
        ggplot(
            plot_df, aes(x = bwnet.colors, y = abs_avg_log2FC, color = bwnet.colors))
                + geom_boxplot()
                + scale_color_manual(values = color_mapping) 
                + theme(
                legend.position = "none",      
                axis.text.x = element_text(angle = 45, hjust = 1)  
                )
        , paste0(output_dir, "_boxplot", extension_plot))

    save_plot(
        ggplot(
            plot_df, aes(x = bwnet.colors, y = avg_log2FC, color = bwnet.colors))
                + geom_violin(alpha = 0.2) 
                + geom_jitter(alpha = 0.8, size = 0.2) 
                + scale_color_manual(values = color_mapping) 
                + theme(
                legend.position = "none",           
                axis.text.x = element_text(angle = 45, hjust = 1) 
                ),
        paste0(output_dir, "_violin", extension_plot))

}


histogram_plot_significance <- function(bwnet, norm_counts, column_data, extension_plot, ...) {
    message("histogram_plot significance")

    gene_significance <- load_significance(norm_counts, ...)

    plot_df <- merge(data.frame(bwnet$colors), gene_significance,
        by = "row.names", all = FALSE)
        
    plot_df$bwnet.colors <- factor(plot_df$bwnet.colors)
    for (column in c("superiorfrontal_thickness", "frontal_cortex_thickness")) {
        # Plot (with correct colors)
        color_mapping <- setNames(c(levels(plot_df$bwnet.colors)), levels(plot_df$bwnet.colors))
        save_plot(
            ggplot(
                plot_df, aes(x = .data[["bwnet.colors"]], y = abs(.data[[column]]), color = .data[["bwnet.colors"]]))
                    + geom_boxplot()
                    + scale_color_manual(values = color_mapping) 
                    + theme(
                    legend.position = "none",            # Remove the legend
                    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
                    )
            , paste0(output_dir, column, "_boxplot", extension_plot))

        save_plot(
            ggplot(
                plot_df, aes(x = .data[["bwnet.colors"]], y = .data[[column]], color = .data[["bwnet.colors"]]))
                    + geom_violin(alpha = 0.2) 
                    + geom_jitter(alpha = 0.8, size = 0.2) 
                    + scale_color_manual(values = color_mapping) 
                    + theme(
                    legend.position = "none",            # Remove the legend
                    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
                    ),
            paste0(output_dir, column, "_violin", extension_plot))
    }

}


membership_log2fc_scatter <- function(bwnet, norm_counts, column_data, cluster, extension_plot, ...) {
    message("membership_log2fc_scatter")

    module_membership_measure <- load_module_membership(bwnet, norm_counts) 
    deg_results <- load_log2fc_df(cluster, ...)
    plot_df <- merge_dataframes(data.frame(bwnet$colors), deg_results, module_membership_measure)

    for (module in unique(bwnet$colors)) {

        # Pull out list of genes in that module
        submod <- plot_df %>%
            subset(bwnet.colors == module)

        # membership score
        module_membership_column_abs <- paste0("ME", module, "_abs")
        module_membership_column <- paste0("ME", module)
        submod[[module_membership_column_abs]] <- abs(submod[[module_membership_column]])
        plot_df[[module_membership_column_abs]] <- abs(plot_df[[module_membership_column]])

        # Skip modules that are too small
        if (nrow(submod) < 9)
            next

        # scatter plot of the module genes
        save_plot(ggplot(submod, aes(x = .data[[module_membership_column]],
            y = .data[["avg_log2FC"]])) + xlim(-1, 1) + geom_point(color = module, alpha = 0.5) +
            geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
            x) + ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
            after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
            label.x = "right", label.y = "top") + theme_grey()
            + labs(caption = paste("number of genes in module:", nrow(submod))), paste0(output_dir,
            "membership_log2fc_scatter/", module, "_scatter_plot", extension_plot))

        # scatter plot of all the genes
        save_plot(ggplot(plot_df, aes(x = .data[[module_membership_column]],
            y = .data[["avg_log2FC"]])) + xlim(-1, 1) + geom_point(color = module, alpha = 0.5) +
            geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
            fill = "white") + ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
            after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
            label.x = "right", label.y = "top") + theme_grey(), paste0(output_dir,
            "membership_log2fc_scatter_all/", module, "_scatter_plot", extension_plot))

        # scatter plot of the module genes
        save_plot(ggplot(submod, aes(x = .data[[module_membership_column_abs]],
            y = .data[["abs_avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
            geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
            x) + ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
            after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
            label.x = "right", label.y = "top") + theme_grey()
            + labs(caption = paste("number of genes in module:", nrow(submod))), paste0(output_dir,
            "membership_log2fc_scatter_abs/", module, "_scatter_plot", extension_plot))
        
        # scatter plot of all the genes
        save_plot(ggplot(plot_df, aes(x = .data[[module_membership_column_abs]],
            y = .data[["abs_avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
            geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
            fill = "white") + ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
            after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
            label.x = "right", label.y = "top") + theme_grey(), paste0(output_dir,
            "membership_log2fc_scatter_all_abs/", module, "_scatter_plot", extension_plot))

    }
}


significance_log2fc_scatter <- function(bwnet, norm_counts, column_data, cluster, extension_plot, ...) {
    message("significance_log2fc_scatter")

    gene_significance <- load_significance(norm_counts, ...)
    deg_results <- load_log2fc_df(cluster, ...)
    plot_df <- merge_dataframes(data.frame(bwnet$colors), deg_results, gene_significance)
    

    for (region_significance_column in  c("superiorfrontal_thickness", "frontal_cortex_thickness")) {

        region_significance_column_abs <- paste0(region_significance_column, "_abs")
        plot_df[[region_significance_column_abs]] <- abs(plot_df[[region_significance_column]])

        for (module in unique(plot_df$bwnet.colors)) {

            submod <- plot_df %>%
                subset(bwnet.colors == module)

            # scatter plot of all the genes
            save_plot(ggplot(submod, aes(x = .data[[region_significance_column]],
                y = .data[["avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
                geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
                fill = "white") + ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                label.x = "right", label.y = "top") + theme_grey() + theme(legend.position = "none") 
                + labs(caption = paste("number of genes in module:", nrow(submod))), paste0(output_dir,
                deg_to_use, "significance_log2fc_scatter_all/", region_significance_column, "_", module, "_scatter_plot", extension_plot))
            
            # scatter plot of all the genes
            save_plot(ggplot(submod, aes(x = .data[[region_significance_column_abs]],
                y = .data[["abs_avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
                geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
                fill = "white") + ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                label.x = "right", label.y = "top") + theme_grey() 
                + labs(caption = paste("number of genes in module:", nrow(submod))), paste0(output_dir,
                deg_to_use, "significance_log2fc_scatter_all_abs/", region_significance_column, "_", module, "_scatter_plot", extension_plot))

        }
    }

    
}


significance_membership_scatter <- function(bwnet, norm_counts, column_data, hub_genes_threshold, cluster, extension_plot, ...) {
   
    gene_significance <- load_significance(norm_counts, ...)
    module_membership_measure <- load_module_membership(bwnet, norm_counts) 
    plot_df <- merge_dataframes(data.frame(bwnet$colors), module_membership_measure, gene_significance)

    for (module in unique(bwnet$colors)) {

        submod <- plot_df %>%
            subset(bwnet.colors == module)

        # membership score
        module_membership_column_abs <- paste0("ME", module, "_abs")
        module_membership_column <- paste0("ME", module)

        submod[[module_membership_column_abs]] <- abs(submod[[module_membership_column]])
        plot_df[[module_membership_column_abs]] <- abs(plot_df[[module_membership_column]])

        for (column in c("superiorfrontal_thickness", "frontal_cortex_thickness")) {

            save_plot(ggplot(submod, aes(x = .data[[module_membership_column]],
            y = .data[[column]])) + xlim(-1, 1) + geom_point(color = module, alpha = 0.5) +
            geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
                x) + ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
            after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
            label.x = "right", label.y = "top") + theme_grey() 
            + labs(caption = paste("number of genes in module:", nrow(submod))), paste0(output_dir,
            "significance_membership_scatter/", module, "_", column, "_scatter_plot", extension_plot))

            # scatter plot of all the genes
            save_plot(ggplot(plot_df, aes(x = .data[[module_membership_column]],
            y = .data[[column]])) + xlim(-1, 1) + geom_point(color = module, alpha = 0.5) +
            geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
                fill = "white") + ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
            after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
            label.x = "right", label.y = "top") + theme_grey(), paste0(output_dir,
            "significance_membership_scatter_all/", module, "_", column, "_scatter_plot", extension_plot))

            # scatter plot of the module genes
            save_plot(ggplot(submod, aes(x = .data[[module_membership_column_abs]],
            y = abs(.data[[column]]))) + geom_point(color = module, alpha = 0.5) +
            geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
                x) + ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
            after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
            label.x = "right", label.y = "top") + theme_grey()
            + labs(caption = paste("number of genes in module:", nrow(submod))), paste0(output_dir,
            "significance_membership_scatter_abs/", module, "_", column, "_scatter_plot", extension_plot))
            
            # scatter plot of all the genes
            save_plot(ggplot(plot_df, aes(x = .data[[module_membership_column_abs]],
            y = abs(.data[[column]]))) + geom_point(color = module, alpha = 0.5) +
            geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
                fill = "white") + ggpmisc::stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
            after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
            label.x = "right", label.y = "top") + theme_grey(), paste0(output_dir,
            "significance_membership_scatter_all_abs/", module, "_", column, "_scatter_plot", extension_plot))
        }
    }
}


correlation_avglog2fc_scatter  <- function(bwnet, norm_counts, column_data, cluster, extension_plot, ...) {
    
    message("correlation_avglog2fc_scatter")

    compute_avg_log2fc <- function(module) {

        submod <- plot_df %>%
            subset(bwnet.colors == module)

        submod <- submod[, c("avg_log2FC", "gene")]

        # if (nrow(submod) < 9) next

        return(mean(abs(submod$avg_log2FC)))
    }
        
    traits <- load_mri_traits(type = "zscore", ...)

    heatmap_data <- bwnet$MEs %>%
        merge(traits, by = "row.names") %>%
        tibble::column_to_rownames(var = "Row.names")

    correlation <- WGCNA::corAndPvalue(
        y = heatmap_data[, (length(heatmap_data) - ncol(traits) + 1):length(heatmap_data)],
        x = heatmap_data[, 1:(length(heatmap_data) - ncol(traits))])$cor

    deg_results <- load_log2fc_df(cluster, ...)

    plot_df <- tibble::column_to_rownames(merge(data.frame(bwnet$colors), deg_results,
        by = "row.names", all = FALSE), var = "Row.names")
    plot_df[["gene"]] <- rownames(plot_df)

    avg_log2fc_scores <- sapply(unique(bwnet$colors), compute_avg_log2fc(module))
    
    for (region in colnames(correlation)) {
        submod <- data.frame(abs(correlation[, region]), avg_log2fc_scores)
        colnames(submod)[1] <- "correlation_abs"
        colnames(submod)[2] <- "avg_log2fc_scores_abs"
        submod$color <- substring(rownames(submod), 3)
        save_plot(
            ggplot(
                submod, aes(x = correlation_abs, y = avg_log2fc_scores_abs)) +
                geom_point(alpha = 0.5) +
                geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
                ggpmisc::stat_poly_eq(
                    formula = y ~ x,
                    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
                    parse = TRUE,
                    size = 5,
                    color = "black",
                    label.x = "right",
                    label.y = "top"
                ) +
                ggrepel::geom_text_repel(aes(label = color), size = 3) +
                theme_grey() + 
                theme(legend.position = "none"),
            paste0(output_dir, "correlation_avg_log2fc_scores_scatter/", region, "_scatter_plot", extension_plot))
    }
}

# Define function for corr_matrix
corr_matrix <- function(bwnet, norm_counts, column_data) {
message("corr_matrix")

    # TODO: check if implemented
    cor_matrix <- cor(bwnet$MEs)


    wb <- write_on_excel("correlation", as.data.frame(cor_matrix))
    openxlsx::saveWorkbook(wb, paste0(output_dir, "module_correlation.xlsx"), overwrite = TRUE)

    # Define color palette
    my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

    # Heatmap
    png(file = paste0(output_dir, "modules_correlation_heatmap", extension_plot))
    heatmap(cor_matrix)
    dev.off()
    message("plot saved in: ", output_dir, "modules_correlation_heatmap", extension_plot)
}

significance_membership_model <- function(bwnet, norm_counts, column_data, ...) {
    message("significance_membership_scatter")

    module_membership_measure <- load_module_membership(bwnet, norm_counts)  # ahh ho invertito e quindo ho una direzione divers
    gene_significance <- load_significance(norm_counts, ...)
    plot_df <- merge_dataframes(data.frame(bwnet$colors), gene_significance, module_membership_measure)

    for (column in c("superiorfrontal_thickness", "frontal_cortex_thickness")) {

        formula <- as.formula(paste(column, "~", paste(colnames(plot_df[, colnames(module_membership_measure)]), collapse = " + ")))
        model <- lm(formula, data = plot_df)
        all_model <- data.frame(summary(model)$coefficients)

        single_model <- list()

        for (module in unique(colnames(module_membership_measure))) {

            # Pull out list of genes in that module
            submod <- plot_df %>%
                subset(bwnet.colors == gsub("ME", "", module))

            model <- lm(submod[[column]] ~ submod[[module]])
            single_model[[gsub("ME", "", module)]] <- c("intercept" = summary(model)$coefficients[1, 1], summary(model)$coefficients[2, ], "r.squared" = summary(model)$r.squared) 

        }

        xlsx::write.xlsx(all_model, file = paste0(output_dir, "linear_model_significance_membership_all_module_", column, ".xlsx"))
        xlsx::write.xlsx(data.frame(do.call("rbind", single_model)), file = paste0(output_dir, "linear_model_significance_membership_single_module_", column, ".xlsx"))   
    }
}