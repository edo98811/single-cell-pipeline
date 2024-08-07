
#finds the directories that match the pattern
find_matching_directories <- function(root_directory, pattern, message=TRUE) {
  # Use list.dirs to recursively list all subdirectories
  all_directories <- list.dirs(root_directory, full.names = TRUE, recursive = TRUE) 

  # Create an empty list to store the matching directory paths
  matching_directories <- list()

  # Iterate through the directory paths and filter matching names
  for (dir_path in all_directories) {
    dir_name <- basename(dir_path)
    if (grepl(pattern, dir_name, fixed = TRUE)) {
      dir_path <- gsub("//", "/", dir_path)
      matching_directories <- c(matching_directories, list(dir_path))
      if (message) message(paste0("found ", dir_path))
    }
  }

  return(matching_directories)
}

# retuns the number of principal components that explain more than a given variabce
analyze_explained_variance <- function(seurat_object, desired_variance, reduction_to_inspect = "pca", default=16) {

  stdev <- seurat_object[[reduction_to_inspect]]@stdev
  for (i in seq_along(stdev)) {
    if (stdev[i] > desired_variance) {
      dimensions <- i  # Store the index if the condition is met
    }
  }
  if (!is.numeric(dimensions)) {cat("Failed to identify correct number of dimensions, using default: 16"); 
    dimensions <- default}
  else cat("Dimensions selected: ", dimensions, " from reduction: ", reduction_to_inspect,"\n")
  
  return(dimensions)
}

# remove parts of the string to keep only te subject name
remove_parts <- function(input_string) {
  # Remove "/input/" from the input string
  output_string <- gsub(paste0(project_folder, "input/"), "", input_string)
  
  # Remove "_filtered_feature_bc_matrix" from the updated string
  output_string <- gsub("_filtered_feature_bc_matrix", "", output_string)
  
  return(trimws(output_string))
}

# saves a plot in the desired dimention in the desired path
save_plot <- function(plotting_function, plotname, x=7, y=7, title=NA){
  check_packages(c("ggplot2", "stringr"))
  
  plot <- force(plotting_function)
  
  # Conditions  on specific plots
  plot_func_name <- str_split(as.character(substitute(plotting_function)), "[(]")[[1]]
  if (plot_func_name =="DotPlot") {
    plot <- plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    plot <- plot + theme(axis.text.x = element_text(size = 6)) 
  }
  if (plot_func_name =="DoHeatmap") {
    plot <- plot + theme(axis.text.y = element_text(size = 6)) 
  }
  #if (is.na(title)) plot + ggtitle(CreateTitle(ggtitle(plot)))
  #else plot + ggtitle(CreateTitle(title))

  # Save plot
  ggsave(plotname, plot = plot, width = x, height = y, bg = "white", create.dir=TRUE)
  message(paste0("plot saved in: ", plotname))
}

# reads from a file in a defined format (for more info...) and saves the patient info i a dataframe containing ID - patology
find_conditions <- function(file_path) {

  # Read lines from the text file
  lines <- readLines(file_path)

  # Split lines based on the separator '-'
  split_lines <- strsplit(lines, " - ", fixed = TRUE)

  # Create a dataframe from the split lines
  df <- do.call(rbind, lapply(split_lines, function(x) {
    data.frame(
      subject = x[1],
      condition = trimws(x[2]),
      stringsAsFactors = TRUE
    )
  }))
  return(df)
}

# delete files according to pattern
delete_files <- function(folder_path, word) {
  # List all files in the specified folder
  files <- list.files(folder_path, full.names = TRUE)
  
  # Filter files containing the word "dimensionality"
  files_to_delete <- files[grep(word, files)]
  
  # Delete files containing "dimensionality"
  if (length(files_to_delete) > 0) {
    file.remove(files_to_delete)
    cat("Files deleted.\n")
  } else {
    cat("No files found.\n")
  }
}

# Function to start timing
start_timing <- function() {
  start_time <<- Sys.time()  # Assign current time to 'start_time' using <<- to assign in the global environment
}

# Function to end timing and save the result to a text file
end_timing <- function(filename = paste0(project_folder, "output/time_log.txt")) {
  end_time <- Sys.time()  # Record the end time
  elapsed_time <- end_time - start_time  # Calculate elapsed time
  
  # Write elapsed time to a text file
  cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), " End_time:",
   format(end_time, "%Y-%m-%d %H:%M:%S"), " Elapsed time:", elapsed_time,
    "\n", file = filename, append = TRUE)
}
#format(elapsed_time, "%M:%OS", trim = TRUE)

#to define the folder in which the results are
define_output_folder <- function(runID){
  output_folder <- paste0(project_folder, "output_4/", runID, "/")
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive=TRUE)
  return(output_folder)
}

#old plot function
unified_plots <- function(seurat_objects) {
  plot_list <- lapply(seurat_objects, function(seurat_object) {
    DimPlot(seurat_object, reduction = "umap")  # Modify group.by as needed
  })
  plot_grid <- do.call("grid.arrange", plot_list)
  ggsave("/output/cluster_plot_all.png", plot = plot_grid, width = 40, height = 40, limitsize = FALSE)  
  
  
  
  # varibale features plot
  plot_list <- lapply(seurat_objects, function(seurat_object) {
    
    top10 <- head(VariableFeatures(seurat_object), 10)
    
    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(seurat_object)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    plot1 + plot2
    
  })
  
  plot_grid <- do.call("grid.arrange", plot_list)
  ggsave("/output/variable_features_plot_all.png", plot = plot_grid, width = 40, height = 40, limitsize = FALSE) 
  
  
  # elbow plot
  plot_list <- lapply(seurat_objects, function(seurat_object) {
    pca_results <- seurat_object[["pca"]]
    ElbowPlot(pca_results, ndims = 30)
  })
  
  plot_grid <- do.call("grid.arrange", plot_list)
  
  # Save the plot to an image file
  # Replace "plot_image.png" with the desired file name and format (e.g., "plot_image.png" or "plot_image.pdf")
  ggsave("/output/elbowplot_plot_all.png", plot = plot_grid, width = 40, height = 40, limitsize = FALSE)  
}


DoBarplot <- function() {
  if (length(names(df)) > 9) stop("too many markers provided")
  for (c in names(df)) {
    barplot(t(as.matrix(df[[c]])), beside = TRUE, col = rainbow(length(unique(df$Cluster))),
            names.arg = unique(df$Cluster), xlab = "Cluster", ylab = "Expression Value",
            main = "Grouped Bar Plot by Cluster")
  }
}

CreateTitle <- function(string) {
  
  words <- strsplit(string[[1]], "_")
  capitalized_words <- toupper(substring(words, 1, 1))
  capitalized_words <- paste0(capitalized_words, substring(words, 2))
  title_case_string <- paste(capitalized_words, collapse = " ")
  return(title_case_string)
}

# To update the log od the experiment
update_text_file <- function(string1, string2) {
  if (FALSE) {
  string1 <- unlist(strsplit(string1, "/"))[length(strsplit(string1, "/")[[1]])]
  
  message <- paste0(string1, ' - content: ', string2)
  
  file_path <- paste0(output_folder, "folder_content.txt")
  
  if (!file.exists(file_path)) {
    # Create the file if it doesn't exist and write the desired line of text
    writeLines(message, file_path)
    writeLines("\n", file_path)
    return(invisible())
  }
  # Read the content of the file
  file_content <- readLines(file_path)
  
  # Check if string1 is present in the file
  string1_index <- grep(string1, file_content)
  
  if (length(string1_index) > 0) {
    # Replace the line containing string1 with the updated content
    file_content[string1_index] <- message
  } else {
    # Append a new line with the updated content
    file_content <- c(file_content, message, "\n")
  }
  
  # Write the updated content back to the file
  writeLines(file_content, file_path)
  return(invisible())
  }
}

# To set up the output folder
set_up_output <- function(output_dir, message = "") {
  
  # update_text_file(output_dir, message)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
  return(output_dir)
}

# Helper function for wgcna old
load_mri_data <- function() {
  
  
  subjects <- list(PD_001 = "02_082",
                   PD_002 = "02_084",
                   PD_005 = "02_074",
                   PD_007 = "02_096",
                   PD_008 = "02_097",
                   PD_009 = "02_105",
                   PD_012 = "02_104",
                   PD_016 = "02_108",
                   PD_017 = "02_115")
  
  columns <- list(rh_aparc_area.txt = c("rh_superiorfrontal_area"),
                  lh_aparc_area.txt = c("lh_superiorfrontal_area"),
                  rh_aparc_volume.txt = c("rh_superiorfrontal_volume"),
                  lh_aparc_volume.txt = c("lh_superiorfrontal_volume"),
                  rh_aparc_thickness.txt = c("rh_superiorfrontal_thickness"),
                  lh_aparc_thickness.txt = c("lh_superiorfrontal_thickness"))
  
  data_list <- list()
  for (file in names(columns)) {
    data <- rownames_to_column(read.table(file = paste0(data_folder, "stats_file_redo/", file), 
                                          header = TRUE, row.names=1), 
                               var ="Row.names")
    
    # Find the real name of the subjects (for each file, even though the should be the same)
    for (i in 1:nrow(data)) {
      for (subject in names(subjects)) {
        if(grepl(subjects[[subject]], data[i,"Row.names"])) 
          subjects[[subject]] <- data[i,"Row.names"]
      }
    }
    
    # create the corrected subject id dataframe
    subject_df <- data.frame(t(data.frame(subjects)))
    colnames(subject_df)[1] <- "Row.names"
    subject_df <- rownames_to_column(subject_df, var ="subject_id")
    
    # 2 non ci sono 
    
    ## select only wanted subjects and colnames ----
    # select the subjects to keep
    to_keep <- data[data[, "Row.names"] %in% unlist(subjects, use.names = FALSE), ]
    
    # selects the columns (traits) to keep
    to_keep <- to_keep[,colnames(to_keep) %in% c(columns[[file]], "Row.names")]
    
    # set rownames to null and merged on Row.names column (subject names), 
    # then delete the Row.names column (the subject names are not needed)
    row.names(to_keep) <- NULL
    to_keep <- merge(to_keep, subject_df, by = "Row.names")
    to_keep$Row.names <- NULL
    
    # add to the list, setting subject id as rownames
    data_list[[tools::file_path_sans_ext(file)]] <- column_to_rownames(to_keep, var="subject_id")
    
  }
  
  # create only one dataframe and compute the metrics
  traits <- bind_cols(data_list)
  traits$superiorfrontal_area <- traits$lh_superiorfrontal_area + traits$rh_superiorfrontal_area
  traits$superiorfrontal_thickness <- rowMeans(traits[c("lh_superiorfrontal_thickness", "rh_superiorfrontal_thickness")])
  traits$superiorfrontal_volume <- traits$lh_superiorfrontal_volume + traits$rh_superiorfrontal_volume
  
  return(traits)
}

# old
ttest_mridata_new <- function(filename = "ttest_results_onlyneeded_without_outliers_1") {
  
  check_packages(c("broom", "tidyverse", "openxlsx"))
  
  columns <- list(rh_aparc_area.txt = c(
    "rh_superiorfrontal_area"
  ),
  lh_aparc_area.txt = c(
    "lh_superiorfrontal_area"
  ),
  rh_aparc_volume.txt = c(
    "rh_superiorfrontal_volume"
  ),
  lh_aparc_volume.txt =  c(
    "lh_superiorfrontal_volume"
  ),
  rh_aparc_thickness.txt = c(
    "rh_superiorfrontal_thickness"
  ),
  lh_aparc_thickness.txt =  c(
    "lh_superiorfrontal_thickness"
  ))
  
  # save all the regions data in a single dataframe
  data_all <- do.call("cbind", 
                  lapply(names(columns),function(file) {
                    read.table(file = paste0(data_folder, "stats_file_redo/", file), 
                               header = TRUE, row.names=1)
                  })
  )
  
  # save all the column names in a single list
  regions_all <- unlist(columns, use.names=FALSE)
  # regions_all <- colnames(data)
  
  data <- list()
  # Iterate though all the regions that needed to be checked both in left and right (should be all but...)
  for(region in 
      intersect(
        lapply(regions_all[grep("^lh_", regions_all)], function(x) {substring(x, 4, nchar(x))}) ,
        lapply(regions_all[grep("^rh_", regions_all)], function(x) {substring(x, 4, nchar(x))})
      )
  ) {
    # If the region is a volumne
    data[[paste0("lh_", region)]] <- data_all[[paste0("lh_", region)]]
    data[[paste0("rh_", region)]] <- data_all[[paste0("rh_", region)]]
    if (grepl("_volume$", region)) {
      data[[region]] <- data_all[[paste0("lh_", region)]] + data_all[[paste0("rh_", region)]]
      # if it is an area or thickness
    } else if (grepl("_area$", region) || grepl("_thickness$", region)) {
      data[[region]] <- rowMeans(data_all[c(paste0("lh_", region), paste0("rh_", region))])
    } else warning(paste0("column not correctly constructed: ", region))
  }
  rm(data_all)
  data <- data.frame(data)
  
  data <- data[!rownames(data) %in% c("02_073_Michael_Rolle_20210802_131203", 
                                      "02_092Sandra_Pilz_20220506_114256", 
                                      "02_097Hans_Lahr_20220624_113719"), ]
  
  # get PD and non PD subjects
  dataset_PD <- data[grep("^02", row.names(data)), ] 
  dataset_CG <- data[grep("^03", row.names(data)), ]
  
  
  # I dont use regions_all because some were added in the previous step
  # Iterate through all the regions and compute the ttest
  ttest_results <- data.frame(
    t(
      sapply(colnames(data), function(col) {
        #message("for column: ", col, "n outliers PD ", remove_outliers(dataset_PD[[col]]))
        #message("for column: ", col, "n outliers CG ", remove_outliers(dataset_CG[[col]]))
        ttest_result <- t.test(remove_outliers(dataset_PD, col), remove_outliers(dataset_CG, col), alternative = "l")
        ttest_summary <- tidy(ttest_result)
        #ttest_summary["outlier_PD"] <- remove_outliers(dataset_PD[[col]])
        #ttest_summary["outlier_CG"] <- remove_outliers(dataset_CG[[col]])
        
        return(ttest_summary)
      }))
  )

  # Save data to excel
  library("xlsx")
  write.xlsx(ttest_results[order(as.numeric(ttest_results$p.value)),], paste0("MRIdata/", filename, ".xlsx"))
  write.xlsx(data, "MRIdata/regions.xlsx")
  message("results saved in: ", paste0("MRIdata/", filename, ".xlsx"))
  
  warning("raincloud to implement")
  return()
  
  # raincloud plot function definition
  plot_raincloud <- function(df, column) {
    check_packages(c("ggdist"))
    
    save_plot(
      ggplot(df, aes(x = pathology, y = value, fill=pathology)) +
        
        ggdist::stat_halfeye(
          adjust = 0.5,
          justification = -.2, 
          .width = 0) +
        
        geom_boxplot(
          width= 0.15, 
          alpha = 0.6) +
        
        theme(legend.position="none"), 
      paste0("mri_data/",column,".png"))
  }
  
  # raincloud plot on selected regions
  for (column in c(unlist(columns, use.names = FALSE), "superiorfrontal_area", "superiorfrontal_thickness", "superiorfrontal_volume")) {
    
    df <- data.frame(value = c(unlist(dataset_02_complete[,column]), unlist(dataset_03_complete[,column])),
                     pathology = c(rep("PD", each = length(dataset_02_complete[,column])),
                                   rep("non_PD", each = length(dataset_03_complete[,column]))))
    
    plot_raincloud(df, column)
    
  }
  
}

# old
remove_outliers <- function(df, column, n = 0) {
  
  # if n is set to 0 it does not do anything and returns the desired column
  if (n == 0) return(df[[column]])
  
  # else it selects the column and computes the outlier threshold (function of standard dviation)
  column_to_test <- df[[column]]
  outlier_threshold <- n * sd(column_to_test)
  
  # Deletes the outliers
  # l <- length(column_to_test[column_to_test < mean(column_to_test) - outlier_threshold]) +  length(column_to_test[column_to_test > mean(column_to_test) + outlier_threshold])
  column_to_test <- column_to_test[column_to_test > (mean(column_to_test) - outlier_threshold)]
  column_to_test <- column_to_test[column_to_test < (mean(column_to_test) + outlier_threshold)]
  
  # Prints the rownames of the deleted outliers
  names1 <- rownames(df[df[,column] < (mean(df[,column]) - outlier_threshold), ])
  names2 <- rownames(df[df[,column] > (mean(df[,column]) + outlier_threshold), ])
  message("ouliers for column ", column, ": ", paste(names1, names2, collapse = ", "))
  
  # Returns column
  return(column_to_test)
  # return (l)
}

# Function wrapper with futures and sink old
run_asyncronously <- function(func, ..., output_file = "output.txt") {
  
  # load liraries 
  check_packages("future")
  
  # Open a connection to the output file
  sink(output_file)
  
  # Se max object limit
  options(future.globals.maxSize= 4194304000)
  
  # Plan to use a multicore backend for futures
  plan(multicore)
  
  # Run the function asynchronously
  future(func(...))
  
  # Close the connection to the output file
  sink()
  
  # Give a message to indicate that the function is running in the background
  cat("Function is running in the background. Check ", paste0(getwd(), "output.txt"), " for progress./n")
}


load_rcc_data <- function(file_path) {

  # library(data.table)
  # file_path <- "C:/Users/Edoardo/Desktop/Seurat/input/open  access data/GSE227990_RAW/GSM7111949_20201203_30102529261121-01_1_01_01.RCC/GSM7111949_20201203_30102529261121-01_1_01_01.RCC"
  file_path <- "C:/Users/Edoardo/Desktop/Seurat/input/open  access data/GSE216281_counts.txt/GSE216281_counts.txt"
  # Read the .rcc file
  # rcc_data <- fread(file_path, skip = 11, fill = TRUE)  # Skip the header part if necessary

  # Read in entire file
  # tirosh <- read.delim(file_path, header = T, stringsAsFactors = F)
  data <- read.table(file_path, header = TRUE, sep = " ")

  # Extract the count data
  # gene_names <- rcc_data$ProbeName  # Adjust according to actual column name
  # counts <- rcc_data$Count          # Adjust acording to actual column name

  # Create a matrix
  count_matrix <- data.frame(Gene = gene_names, Count = counts)
}

load_txt_data <- function() {
    data <- read.table(file_path, header = TRUE, sep = " ")

    # Load the packages
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    
    # Convert Ensembl IDs to gene names
    rownames(data) <- mapIds(org.Hs.eg.db, keys = rownames(data), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
    seurat_obj <- CreateSeuratObject(counts = data)
}