# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


source("scripts/mri_utils.r", local = TRUE)

# Load data for the mri needed for the wgcna function
mri_for_wgcna <- function(data_source = "aparc", regions = c(), wgcna_subjects = list()) {
    
    library(stringr)

    columns <- create_data_description(regions)
    # save all the column names in a single list
    
    # regions_all <- colnames(data)
    if (data_source == "samseg") data <- load_data_samseg(columns)
    else if (data_source == "aparc") data <- filter_data(load_data_aparc(columns), columns)
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
zscore_for_wgcna <- function(data_source = "aparc", regions = c(), wgcna_subjects = list()) {
    
    library(stringr)

    columns <- create_data_description(regions)
    # save all the column names in a single list
    
    # regions_all <- colnames(data)
    if (data_source == "samseg") data <- load_data_samseg(columns)
    else if (data_source == "aparc") data <- filter_data(load_data_aparc(columns), columns)
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
    corrected_data <- compute_zscore(data)

    # zscore <- function(value, avg, std_dev) {
    #     score <- (value - avg) / std_dev
    # }

    # corrected_data <- data.frame(apply(data, 2, function(region) {
    #     avg <- mean(region[grep("^03", names(region))])
    #     std_dev <- sd(region[grep("^03", names(region))])
        
    #     z_scores <- sapply(region, function(x) zscore(x, avg = avg, std_dev = std_dev))
    # }))

    # Selecting only files present in wgcna analysis
    corrected_data <- corrected_data[rownames(corrected_data) %in% wgcna_subjects,]
    rownames(corrected_data) <- sapply(rownames(corrected_data), function(name) {return(names(wgcna_subjects[wgcna_subjects == name]))})
    return(corrected_data)

}
