# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de

# Define the info to load
create_data_description <- function(regions) {

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

# To compute the zscore based on control subjects
compute_zscore <- function(data) {

    zscore <- function(value, avg, std_dev) {
        score <- (value - avg) / std_dev
    }

    corrected_data <- data.frame(apply(data, 2, function(region) {
        avg <- mean(region[grep("^03", names(region))])
        std_dev <- sd(region[grep("^03", names(region))])
        
        z_scores <- sapply(region, function(x) zscore(x, avg = avg, std_dev = std_dev))
    }))
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
        parts <- str_split(line, ",//s*", simplify = TRUE)

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
load_data_aparc <- function(columns, source_files = "C:/Users/Edoardo/Desktop/Seurat/input/stats_file_redo") {

    # save all the regions data in a single dataframe
    data_all <- do.call("cbind", lapply(names(columns), function(file) {
        read.table(file = paste0(source_files, file), header = TRUE,
            row.names = 1)
    }))

    return(data_all)
}