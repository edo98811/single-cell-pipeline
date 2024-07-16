# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de

# Main function for statistical tests on mri data
mri_stat_test <- function(filename = "ttest_results_onlyneeded_samseg", data_source = "aparc", regions = c( "superiorfrontal", "caudalmiddlefrontal", "rostralmiddlefrontal")) {

    columns <- create_data_description(regions)
    # save all the column names in a single list
    
    # regions_all <- colnames(data)
    if (data_source == "samseg") data <- load_data_samseg(columns)
    else if (data_source == "aparc") data <- filter_data(load_data_aparc(source, columns), columns)
    else stop("invalid data type")
    # data <- load_data_samseg(columns)

    # adding region of the sensor
    # Try the rostral middle and caudal middle frontal regions
    # data[["frontal_cortex_volume"]] <- data[["superiorfrontal_volume"]] + data[["caudalmiddlefrontal_volume"]] + data[["rostralmiddlefrontal_volume"]]
    # data[["frontal_cortex_thickness"]] <- rowMeans(data[c(paste0("superiorfrontal_thickness"), paste0("caudalmiddlefrontal_thickness"), paste0("caudalmiddlefrontal_thickness"))])
    # data[["frontal_cortex_area"]] <- rowMeans(data[c(paste0("superiorfrontal_area"), paste0("caudalmiddlefrontal_area"), paste0("caudalmiddlefrontal_area"))])

    # Deleting decided outliers
    data <- data[!rownames(data) %in% c("02_073_Michael_Rolle_20210802_131203", 
                                        "02_092Sandra_Pilz_20220506_114256", 
                                        "02_097Hans_Lahr_20220624_113719"), ]
    # get PD and non PD subjects
    dataset_PD <- data[grep("^02", row.names(data)), ]
    dataset_CG <- data[grep("^03", row.names(data)), ]

    test_results <- run_stat_test(dataset_PD, dataset_CG, t.test)
    save_test_results(test_results, filename)

    test_results <- run_stat_test(dataset_PD, dataset_CG, wilcox.test)
    save_test_results(test_results, paste0(filename, "wilc"))

    warning("raincloud to implement")

}

# Main function for anova tests on mri data (to implement)
mri_anova_test <- function(filename = "anova_results", patient_info_file = "MRIdata/patient_info.xlsx", regions =  c( "superiorfrontal", "caudalmiddlefrontal", "rostralmiddlefrontal")) {

    create_data_description(regions)

    data <- filter_data(load_data_aparc(source, columns), columns)
    patient_info <- load_patient_info(patient_info_file)

    # Deleting decided outliers
    data <- data[!rownames(data) %in% c("02_073_Michael_Rolle_20210802_131203", 
                                        "02_092Sandra_Pilz_20220506_114256", 
                                        "02_097Hans_Lahr_20220624_113719"), ]

    # Make the rownames standard (to allow union)
    for (i in seq_along(nrow(data))) 
        rownames(data)[i] <- gsub("^.{3}0?0?", 
                                "",
                                str_extract(
                                    rownames(data)[i], 
                                    "(^[0-9]{2}_[0-9]*)"
                                    )
                                )

    # create unique dataset
    data_w_info <- cbind(
        data,
        patient_info
    )
    
    lapply()
}   

# Runs the statistical test
run_stat_test <- function(dataset_PD, dataset_CG, test = t.test, n = 0, ...) {
    # I dont use regions_all because some were added in the previous step
    # Iterate through all the regions and compute the ttest
    library("broom")

    ttest_results <- data.frame(t(sapply(intersect(colnames(dataset_PD), colnames(dataset_CG)), function(col) {

        ttest_result <- test(remove_outliers(dataset_PD, col, n), remove_outliers(dataset_CG, col, n), alternative = "l")
        ttest_summary <- tidy(ttest_result)

        return(ttest_summary)
    })))

    return(ttest_results)
}

# Runs the anova test
run_anova_test <- function() {

    library("broom")

    test_result <- data.frame(t(sapply(intersect(colnames(dataset_PD), colnames(dataset_CG)), function(col) {

        ttest_result <- test(remove_outliers(dataset_PD, col, n), remove_outliers(dataset_CG, col, n), alternative = "l")
        ttest_summary <- tidy(ttest_result)

        return(ttest_summary)
    })))

    return(test_result)
}

# Define the info to load
create_data_description <- function(regions = c( "superiorfrontal", "caudalmiddlefrontal", "rostralmiddlefrontal")) {

    return(
        list(
            rh_aparc_area.txt = sapply(regions, function(region) {paste0("rh_", region, "_area")}),
            lh_aparc_area.txt = sapply(regions, function(region) {paste0("lh_", region, "_area")}),
            rh_aparc_volume.txt = sapply(regions, function(region) {paste0("rh_", region, "_volume")}),
            lh_aparc_volume.txt =  sapply(regions, function(region) {paste0("lh_", region, "_volume")}),
            rh_aparc_thickness.txt = sapply(regions, function(region) {paste0("rh_", region, "_thickness")}),
            lh_aparc_thickness.txt =  sapply(regions, function(region) {paste0("lh_", region, "_thickness")})
        )
    )
}

# Loads the data on the patients
load_patient_info <- function(filename) {
    library("xlsx")

    return( 
        do.call("rbind", 
            list(
                read.xlsx(filename, 1, header = TRUE),
                read.xlsx(filename, 1, header = TRUE)
                )
            )
        )
}

# To make a raincloud plot with the desired columns( to implement)
raincloud <- function(
    result_folder = "raincloud", 
    data_source = "aparc", 
    regions = c("superiorfrontal", "caudalmiddlefrontal", "rostralmiddlefrontal"),
    wgcna_subjects = list(PD_001 = "02_082", PD_002 = "02_084", PD_005 = "02_074", PD_007 = "02_096",
            PD_008 = "02_097", PD_009 = "02_105", PD_012 = "02_104", PD_016 = "02_108",
            PD_017 = "02_115"),
    columns_to_plot = c("")) {

    # raincloud plot function definition
    plot_raincloud <- function(df, column, folder = "mri_regions") {

        library(ggplot2)
        library(ggdist)

        save_plot(
            ggplot(
                df, 
                aes(x = pathology, y = value, fill = pathology)) 
                + ggdist::stat_halfeye(adjust = 0.5,
                    justification = -0.2, .width = 0) 
                + geom_boxplot(width = 0.15, alpha = 0.6) 
                + theme(legend.position = "none"),
            paste0("MRIdata/",folder , "/", column, ".png")
        )
    }
    columns <- create_data_description(regions)
    # save all the column names in a single list
    
    # regions_all <- colnames(data)
    if (data_source == "samseg") data <- load_data_samseg(columns)
    else if (data_source == "aparc") data <- filter_data(load_data_aparc(source, columns), columns)
    else stop("invalid data type")

    # Deleting decided outliers
    data <- data[!rownames(data) %in% c("02_073_Michael_Rolle_20210802_131203", 
                                        "02_092Sandra_Pilz_20220506_114256", 
                                        "02_097Hans_Lahr_20220624_113719"), ]

    # Make the rownames standard (to allow union with wgcna subjects) (to fix)
    for (i in seq_len(nrow(data))) {
        # try catch block -> if a line encounters an error it is simply deleted (to check if it is ok)
        tryCatch({
            rownames(data)[i] <- str_extract(
                                rownames(data)[i], 
                                "(^[0-9]{2}_[0-9]*)"
                                )
        }, warning = function(w) {
            data <- data[-c(i), ]
        })
    }

    # for normal
    dataset_PD <- data[grep("^02", row.names(data)), ]
    dataset_CG <- data[grep("^03", row.names(data)), ]
    dataset_omics <- data[rownames(data) %in% wgcna_subjects, ]

    # raincloud plot on selected regions
    for (column in colnames(data)) {

        # Create dataframe with value as pathology and patholoy as
        df <- data.frame(
            value = c(
                unlist(dataset_PD[, column]), 
                unlist(dataset_CG[, column]), 
                unlist(dataset_omics[, column])), 
            pathology = c(
                rep("PD", each = length(dataset_PD[, column])), 
                rep("non_PD", each = length(dataset_CG[, column])),
                rep("omics_PD", each = length(dataset_omics[, column])))
            )
        plot_raincloud(df, column, folder = "mri_regions")
    }
    
    # For zscore
    data_zscore <- compute_zscore(data)
    dataset_PD <- data_zscore[grep("^02", row.names(data_zscore)), ]
    dataset_CG <- data_zscore[grep("^03", row.names(data_zscore)), ]
    dataset_omics <- data_zscore[rownames(data_zscore) %in% wgcna_subjects, ]

    for (column in colnames(data_zscore)) {

        # Create dataframe with value as pathology and patholoy as
        df <- data.frame(
            value = c(
                unlist(dataset_PD[, column]), 
                unlist(dataset_CG[, column]), 
                unlist(dataset_omics[, column])), 
            pathology = c(
                rep("PD", each = length(dataset_PD[, column])), 
                rep("non_PD", each = length(dataset_CG[, column])),
                rep("omics_PD", each = length(dataset_omics[, column])))
            )
        plot_raincloud(df, column, folder = "mri_regions_zscore")
    }
}

compute_zscore <- function(data) {

    zscore <- function(value, avg, std_dev) {
        score <- (value-avg)/std_dev
    }

    corrected_data <- data.frame(apply(data, 2, function(region) {
        avg <- mean(region[grep("^03", names(region))])
        std_dev <- sd(region[grep("^03", names(region))])
        
        z_scores <- sapply(region, function(x) zscore(x, avg = avg, std_dev = std_dev))
    }))

}
# Saves the results of the statistical test in an excel file
save_test_results <- function(results, filename) {
    library("xlsx")
    
    # Save data to excel
    write.xlsx(results[order(as.numeric(results$p.value)), ], paste0("MRIdata/",
        filename, ".xlsx"))

    # write.xlsx(data, "MRIdata/regions.xlsx")
    message("results saved in: ", paste0("MRIdata/", filename, ".xlsx"))

}

# only the columns specified in the column argument 
filter_data <- function(data_all, columns) {

    data <- list()
    # Iterate though all the regions that needed to be checked both in left and
    # right (should be all but...)
    regions <- unlist(columns, use.names = FALSE)

    for (region in intersect(lapply(regions[grep("^lh_", regions)], function(x) {
        substring(x, 4, nchar(x))
    }), lapply(regions[grep("^rh_", regions)], function(x) {
        substring(x, 4, nchar(x))
    }))) {
        # If the region is a volumne
        data[[paste0("lh_", region)]] <- data_all[[paste0("lh_", region)]]
        data[[paste0("rh_", region)]] <- data_all[[paste0("rh_", region)]]
        if (grepl("_volume$", region)) {
            data[[region]] <- data_all[[paste0("lh_", region)]] + data_all[[paste0("rh_",
                region)]]
            # if it is an area or thickness
        } else if (grepl("_area$", region) || grepl("_thickness$", region)) {
            data[[region]] <- rowMeans(data_all[c(paste0("lh_", region), paste0("rh_",
                region))])
        } else warning(paste0("column not correctly constructed: ", region))
    }

    if (all(c("superiorfrontal_volume", "caudalmiddlefrontal_volume", "rostralmiddlefrontal_volume") %in% names(data))) {
        data[["frontal_cortex_volume"]] <- data[["superiorfrontal_volume"]] + data[["caudalmiddlefrontal_volume"]] + data[["rostralmiddlefrontal_volume"]]
        data[["frontal_cortex_thickness"]] <- rowMeans(data.frame(data[c("superiorfrontal_thickness","caudalmiddlefrontal_thickness","caudalmiddlefrontal_thickness")]))
        data[["frontal_cortex_area"]] <- rowMeans(data.frame(data[c("superiorfrontal_area","caudalmiddlefrontal_area","caudalmiddlefrontal_area")]))
    }

    return(data.frame(data))
}

# If n is specified, remove outliers that are over or below n devstd
remove_outliers <- function(df, column, n = 0) {
    
    # if n is set to 0 it does not do anything and returns the desired column
    if (n == 0)
        return(df[[column]])

    # else it selects the column and computes the outlier threshold (function
    # of standard dviation)
    column_to_test <- df[[column]]
    outlier_threshold <- n * sd(column_to_test)

    # Deletes the outliers l <- length(column_to_test[column_to_test <
    # mean(column_to_test) - outlier_threshold]) +
    # length(column_to_test[column_to_test > mean(column_to_test) +
    # outlier_threshold])
    column_to_test <- column_to_test[column_to_test > (mean(column_to_test) - outlier_threshold)]
    column_to_test <- column_to_test[column_to_test < (mean(column_to_test) + outlier_threshold)]

    # Prints the rownames of the deleted outliers
    names1 <- rownames(df[df[, column] < (mean(df[, column]) - outlier_threshold),])
    names2 <- rownames(df[df[, column] > (mean(df[, column]) + outlier_threshold),])

    message("ouliers for column ", column, ": ", paste(names1, names2, collapse = ", "))

    # Returns column
    return(column_to_test)
    # return (l)
}

# Load samseg informations
load_data_samseg <- function(columns, source_files = paste0(getwd(), "/TablesFreesurfer_all/stat_files/stat_files")) {
    
    # Get list of files in the folder
    file_list <- list.files(source_files, full.names = TRUE)

    # Initialize an empty list to store the results
    results <- list()

    # Process each file
    for (file in file_list) {
        # Get the file name without the path
        file_name <- basename(file)

        # Read and process the file
        measures <- read_samseg_file(file)

        # Store in the results list
        results[[file_name]] <- measures
    }
    
    return(data.frame(t(do.call("cbind", results))))
}

# Function to read and process a single file
read_samseg_file <- function(file_path) {
    library(stringr)

    # Read the file
    lines <- readLines(file_path)

    # Initialize an empty named numeric vector
    measures <- numeric(length(lines))
    names(measures) <- character(length(lines))

    # Process each line
    for (i in seq_along(lines)) {
        # Remove the leading '# Measure ' and the trailing ', mm^3'
        line <- str_trim(lines[i])
        line <- str_remove(line, "^# Measure ")
        parts <- str_split(line, ",\\s*", simplify = TRUE)

        # Extract region name and measure
        region_name <- parts[1]
        measure <- as.numeric(parts[2])

        # Store in the vector
        measures[i] <- measure
        names(measures)[i] <- region_name
    }

    return(measures)
}

# Load the mrirestuls informations
load_data_aparc <- function(source, columns) {

    # save all the regions data in a single dataframe
    data_all <- do.call("cbind", lapply(names(columns), function(file) {
        read.table(file = paste0(data_folder, "stats_file_redo/", file), header = TRUE,
            row.names = 1)
    }))

    return(data_all)
}

# Load data for the mri needed for the wgcna function
mri_for_wgcna <- function(data_source = "aparc", regions = c("superiorfrontal", "caudalmiddlefrontal", "rostralmiddlefrontal"),
    wgcna_subjects = list(PD_001 = "02_082", PD_002 = "02_084", PD_005 = "02_074", PD_007 = "02_096",
        PD_008 = "02_097", PD_009 = "02_105", PD_012 = "02_104", PD_016 = "02_108",
        PD_017 = "02_115")) {
    
    library(stringr)

    columns <- create_data_description(regions)
    # save all the column names in a single list
    
    # regions_all <- colnames(data)
    if (data_source == "samseg") data <- load_data_samseg(columns)
    else if (data_source == "aparc") data <- filter_data(load_data_aparc(source, columns), columns)
    else stop("invalid data type")
    # data <- load_data_samseg(columns)

    # Make the rownames standard (to allow union with wgcna subjects) (to fix)
    for (i in seq_len(nrow(data))) {
        # try catch block -> if a line encounters an error it is simply deleted (to check if it is ok)
        tryCatch({
            rownames(data)[i] <- str_extract(
                                rownames(data)[i], 
                                "(^[0-9]{2}_[0-9]*)"
                                )

        }, warning = function(w) {
            data <- data[-c(i), ]
        })
    }

    # Selecting only files present in wgcna analysis
    data <- data[rownames(data) %in% wgcna_subjects,]
    rownames(data) <- sapply(rownames(data), function(name) {return(names(wgcna_subjects[wgcna_subjects == name]))})
    return(data)

}

# Load data for the mri needed for the wgcna function
zscore_for_wgcna <- function(data_source = "aparc", regions = c("superiorfrontal", "caudalmiddlefrontal", "rostralmiddlefrontal"),
    wgcna_subjects = list(PD_001 = "02_082", PD_002 = "02_084", PD_005 = "02_074", PD_007 = "02_096",
        PD_008 = "02_097", PD_009 = "02_105", PD_012 = "02_104", PD_016 = "02_108",
        PD_017 = "02_115")) {
    
    library(stringr)

    columns <- create_data_description(regions)
    # save all the column names in a single list
    
    # regions_all <- colnames(data)
    if (data_source == "samseg") data <- load_data_samseg(columns)
    else if (data_source == "aparc") data <- filter_data(load_data_aparc(source, columns), columns)
    else stop("invalid data type")
    # data <- load_data_samseg(columns)

    # Make the rownames standard (to allow union with wgcna subjects) (to fix)
    for (i in seq_len(nrow(data))) {
        # try catch block -> if a line encounters an error it is simply deleted (to check if it is ok)
        tryCatch({
            rownames(data)[i] <- str_extract(
                                rownames(data)[i], 
                                "(^[0-9]{2}_[0-9]*)"
                                )
        }, warning = function(w) {
            data <- data[-c(i), ]
        })
    }

    zscore <- function(value, avg, std_dev) {
        score <- (value-avg)/std_dev
    }

    corrected_data <- data.frame(apply(data, 2, function(region) {
        avg <- mean(region[grep("^03", names(region))])
        std_dev <- sd(region[grep("^03", names(region))])
        
        z_scores <- sapply(region, function(x) zscore(x, avg = avg, std_dev = std_dev))
    }))

    # Selecting only files present in wgcna analysis
    corrected_data <- corrected_data[rownames(corrected_data) %in% wgcna_subjects,]
    rownames(corrected_data) <- sapply(rownames(corrected_data), function(name) {return(names(wgcna_subjects[wgcna_subjects == name]))})
    return(corrected_data)

}
