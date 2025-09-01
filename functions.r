# %% [markdown]
# # Introduction
# 
# This notebooks holds all the functions used during the pipeline. <br>
# During the pipeline, at the start of each notebook, we source function.r to load all the functions at once. <br>
# This notebook provides a more readable file, corresponding to the function.r file.

# %% [markdown]
# # Known bugs :
# 
# - Metaclustering with louvain or other clustering that have a limited number of cluster can fail if the max_k is too high compared to the number of cluster. The exact number where it will fail is hard to predict as it seem to depend of random subsetting inside 

# %% [markdown]
# ## rename_files
# 
# Simple function to rename files from a folder. <br>
# Uses str_split_i to identify the code name based on a separator and a position
# 
# **ARGUMENTS**
# - `folder` : character string, path to the folder which contains the files to rename. All files matching the file_type will be renamed with the same method.
# -  `new_name` : character string, new name to be implemented at the start of all the files in the folder. `""` will not add anything before the retrieve code
# -  `file_type` : character string, pattern to identify the file to rename. Expected to contain the file extension, as it will be used at the end of the new files.
# -  `separator` : character string, pattern to identify the separator for str_split_i (usually a special character)
# -  `code_position` : integer value, used to identify the code position compared to `separator`. 1 means before the first separator, 2 after the first and so on. Negative integer can be used to start from the right hand side.
# 
# **VALUE**<br>
# Used to rename files in folder. <br>
# Returns folder invisibly.

# %%
rename_files <- function(folder,
                         new_name = "",
                         file_type = ".tsv",
                         separator = "_",
                         code_position = 2) {
  #Preparing the file list
  files <- list.files(path = folder,
                      pattern = file_type)
  #For each file
  for (fil in files) {
    #Extract the sample code
    code_name <- stringr::str_split_i(fil,
                                      pattern = separator,
                                      i = name_position)
    #Prepare a new file name
    new_file <- paste0(folder, new_name, "_", code_name, file_type)
    #Rename the file
    file.rename(paste0(folder, fil), new_file)
  }
  
  invisible(folder)
}

# %% [markdown]
# ## Prepare_region_prop
# 
# This function takes a folder as an input, and for each file in that folder creates regionprops files in the appropriate new directory.
# 
# **ARGUMENTS**
# - `in_dir` : character string indicating the path to the input folder
# - `out_dir` : character string indicating the path to the regionprops output folder
# - `pattern` : character string indicating the file extension. Only ".tsv" and ".csv" are supported
# - `keep` : Vector containing character strings used to identify, via regular expression the column to keep in the regionprop files
#   
# **VALUE** <br>
# Creates the region prop csv files in prop_dir. <br>
# Returns in_dir invisibly.

# %%
prepare_region_prop <- function(in_dir = in_dir,
                                out_dir = prop_dir,
                                pattern = c(".tsv", ".csv"),
                                keep = keep_prop,
                                use_roi = FALSE) {

  pattern <- match.arg(pattern)
  ###I Reading files:
  #Listing files in folder
  files <- list.files(
    path = in_dir,
    pattern = pattern
  )
  purrr::walk(
    files,
    record_region_prop,
    in_dir = in_dir,
    out_dir = out_dir,
    pattern = pattern,
    use_roi = use_roi
  )

  invisible(in_dir)
}

# %% [markdown]
# ### record_region_prop
# 
# Helper function that reads one file and creates a corresponding regionprops file in the regionprop folder

# %%
record_region_prop <- function(file, in_dir, out_dir, pattern, use_roi) {
  if (pattern == ".csv") {
    df_i <- readr::read_csv(
      paste0(in_dir, file),
      show_col_types = FALSE
    )
  } else if (pattern == ".tsv") {
    df_i <- readr::read_tsv(
      paste0(in_dir, file),
      show_col_types = FALSE
    )
  } else { #returns a error if unsupported extension
    stop(
      "File extension not supported : Supported pattern '.csv' ou '.tsv'"
    )
  }

  if (use_roi) {
    df_i <- extract_roi(df_i, file)
  }

  #We prepare the output full name,
  ## Including the prefix (path to folder),
  ## The name of the original file without the extension (pattern)
  ## and the the new output extension (.csv)
  output_name <- paste0(
    out_dir,
    stringr::str_split_i(
      file,
      pattern = pattern,
      i = 1
    ),
    ".csv"
  )

  ### II Fusion and export
  #Keep only the relevant columns

  #Finding the column
  col_props <- unlist(lapply(keep_prop, function(x) {
    grep(pattern = x, colnames(df_i), value = TRUE)
  }))
  #Subset the dataframe based on thos columns
  df_i <- df_i[, col_props]

  # We create a new dataframe, with names matching the expected output
  df_prop_i <- data.frame(
    "Object" = rownames(df_i),
    "area" = df_i[, 3],
    "centroid.1" = df_i[, 1],
    "centroid.2" = df_i[, 2],
    "axis_major_length" = df_i[, 4],
    "axis_minor_length" = df_i[, 5],
    "eccentricity" = df_i[, 6]
  )
  #We write the new file in the regionprops folder :
  readr::write_csv(
    df_prop_i,
    file = output_name,
    quote = "none"
  )
}

# %% [markdown]
# ## Prepare_intensity
# 
# This function takes a folder as an input, and for each file in that folder creates a intensity file in the intenisties folder.
# 
# **ARGUMENTS**
# - `in_dir` : character string indicating the path to the input folder
# - `out_dir`: character string indicating the path to the regionprops output folder
# - `pattern` : character string indicating the file extension. Only ".tsv" and ".csv" are supported
# - `subset`: Character strings used to identify, via regular expression, the columns contain the cell part to keep for the intensity file
# - `tag_used` : Vector containing character strings used to identify, via regular expression, the columns containing the intensities
# - `panel` : Vector containing character strings used to rename the columns
# 
# **VALUE**<br>
# Creates the corresponding intensity files in the intens_dir. <br>
# Returns in_dir invisibly.

# %%
prepare_intensity <- function(in_dir,
                              out_dir = intens_dir,
                              pattern = c(".csv", ".tsv"),
                              subset = "Entire.cell",
                              tag_used,
                              exclude = "",
                              panel_used,
                              use_roi = FALSE) {
                  
  pattern <- match.arg(pattern)
  ###I Reading files
  #Listing files:
  files <- list.files(path = in_dir,
                      pattern = pattern)

  purrr::walk(
    files,
    record_intensity,
    in_dir = in_dir,
    out_dir = out_dir,
    pattern = pattern,
    subset = subset,
    tag_used = tag_used,
    exclude = exclude,
    panel_used = panel_used,
    use_roi = use_roi
  )

  invisible(in_dir)
}

# %% [markdown]
# ### record_intensity
# 
# Helper function that reads one file and creates a corresponding intensity file in the intensity folder

# %%
record_intensity <- function(
  file,
  pattern,
  in_dir,
  out_dir,
  use_roi,
  tag_used,
  panel_used,
  subset,
  exclude
) {
  if (pattern == ".csv") {
    df_i <- readr::read_csv(
      paste0(in_dir, file),
      show_col_types = FALSE
    )
  } else if (pattern == ".tsv"){
    df_i <- readr::read_tsv(
      paste0(in_dir, file),
      show_col_types = FALSE
    )
  }else {#Returns an error if unsupported extension
    stop("Error = File extension not supported : Supported pattern '.csv' ou '.tsv'") # nolint: line_length_linter.
  }

  if (use_roi) {
    df_i <- extract_roi(df_i, file)
  }

  #We prepare the name of the output file
  output_name <- paste0(
    out_dir,
    stringr::str_split_i(
      file,
      pattern = pattern,
      i = 1
    ),
    ".csv"
  )
  #II Subsetting columns :
  col_tag <- lapply(
    tag_used,
    function(x) {
    #We keep the columns that :
    #1 correspond to the element x of the pannel 
    #2 correspond to the selected cellpart subset
    correct <- grepl(pattern = x, colnames(df_i)) &
      grepl(pattern = subset, colnames(df_i))
    if (!is.null(exclude)) {
      correct <- correct & !grepl(pattern = exclude, colnames(df_i))
    }
    colnames(df_i)[correct]
    }
  ) |>
    unlist()

  #III : PREPARING AND EXPORTING DATAFRAME
  #Subsetting dataframe
  df_intens_i <- df_i[, col_tag]
  #Rename column to match the tag
  colnames(df_intens_i) <- panel_used
  #We add a column Object that contain the "name" of the cell
  df_intens_i$Object <- rownames(df_intens_i)
  #And reorder the dataframe for clarity
  df_ordered <- df_intens_i[, c("Object", panel_used)]

  #And we store the new file as csv file
  readr::write_csv(
    df_ordered,
    file = output_name,
    quote = "none"
  )
}

# %% [markdown]
# ## prepare_panel
# 
# This function is used to prepare a panel.csv file. <br>
# It is also used for its side effect of verifying correct match between panel and tag columns.
# 
# **Arguments**
# - `in_dir` : character string indicating the path to the input folder
# - `out_dir`: character string indicating the path to the regionprops output folder
# - `pattern` : character string indicating the file extension. Only ".tsv" and ".csv" are supported
# - `subset`: Character strings used to identify, via regular expression, the columns contain the cell part to keep for the intensity file
# - `tag` : Vector containing character strings used to identify, via regular expression, the columns containing the intensities
# - `panel` : Vector containing character strings used to rename the columns
# - `accept_incorrect_match` : Logical, set to TRUE if you want to bypass the column check
#   
# **VALUE**<br>
# Creates a panel.csv file in out_dir <br>
# Returns in_dir invisibly.
# 
# **DETAILS**
# Since we define both the tags and the panel, this function does not do a lot. <br>
# We use it to try to do an ultimate verification that the name of the panel is found in the corresponding column before simply writing a csv. <br>
# This functionnality can be bypassed by the argument `accept_incorrect_match`, although it is not advisable unless you are absolutely sure you have identified why the match is incorrect and you are absolutely sure that no column is identified twice by one tag. <br>
# **Skipping this verification can lead to serious misattribution of the data that will not be easy to spot**

# %%
prepare_panel <- function(in_dir = in_dir,
                          out_dir = "./",
                          pattern = c(".tsv", ".csv"),
                          subset = "Entire.cell",
                          tag = tag,
                          exclude = "",
                          panel = panel,
                          accept_incorrect_match = FALSE) {
  pattern <- match.arg(pattern)
  ###I Reading files:
  #Listing files :
  files <- list.files(path = in_dir,
                      pattern = pattern)
  j <- 1 #We only work on the first file
  if (pattern == ".csv") {
    df_i <- read.delim(paste0(in_dir, files[j]),
                       sep = ",")
  } else if (pattern == ".tsv") {
    df_i <- read.delim(paste0(in_dir, files[j]),
                       sep = "\t")
  }
  ### II Select used colums
  col_channel <- unlist(lapply(tag, function(x) {
    #Matching the element of the panel and the cell_part
    correct <- grepl(pattern = x, colnames(df_i)) &
      grepl(pattern = subset, colnames(df_i))
    if (!is.null(exclude)) {
      correct <- correct & !grepl(pattern = exclude, colnames(df_i))
    }

    colnames(df_i)[correct]
  }))
  # Making sure the colnames contains both the panel name and the tag name
  # Bypassable exception
  if (!accept_incorrect_match) {
    for (i in seq_along(tag)) {
      if (!grepl(panel[[i]], col_channel[[i]])) {
        stop("Error, could not find ", panel[[i]], " in matching column :",
             col_channel[[i]], ".\n 
             Please make sure the panel and tag variables are correct")
      }
    }
  }
  #Create the dataframe
  panel_df <- data.frame("channel" = tag,
                         "name" = panel,
                         "keep" = 1) #To be loaded into R
  #And save the dataframe
  write.csv(panel_df,
            file = paste0(out_dir, "panel.csv"),
            quote = FALSE,
            row.names = FALSE)

  invisible(in_dir)
}

# %% [markdown]
# ## prepare_metadata 
# 
# This function takes metadata specified either manually or automatically, and creates a table from it which can be either stored or returned.
# 
# **ARGUMENTS**
# - `output_dir` : character string indicating the path to the folder where the metadata file will be stored.
# - `use_threshold` : logical value, indicating whether to add thresholds in the metadata.
# - `threshold_dir` : character string indicating the path to the folder where the treshold files are located.
# - `keep_lower` and `keep_upper` : logical values, indicating wheter to keep the lower and upper threshold in the metadata.
# - `intens_dir` : path to the intesities folder. Required if `use_threshold = FALSE`
# - `file_extension` : character string indicating the extension of the threshold files. Only ".tsv" or ".csv" are currently supported.
# - `subset` : character string used to recognize, with regular expression, the subset you want to select. "" will include all the files.
# - `tag ` and `panel`: character vectors used to identify and rename the column respectively.
# - `indication_in_name` : logical value, indicating whether to search for the indication in the name.
# - `indication_position` : integer value, indicating the position of the indication within the sample_id (relative to the separator).
# - `separator` : character string, used to identify the separator for the indication using regular expression.
# - `indication` : character vector, used to add manually the indications if "indication_in_name" is TRUE. Ignored if "indication_in_name" is FALSE.
# - `additional_metadata` : list containing vectors with any additional metadata you want to provide. The vector must contain a value for each sample.
# - `store_csv` : logical value, indicating whether to store the output as a ".csv" file. If FALSE, the output will only be returned as a dataframe.
# 
# **VALUE** <br>
# If store_csv is TRUE : write a csv file and returns output_dir invisibly. <br>
# If store_csv is FALSE : returns the data.frame containing the metadata.
# 
# **DETAILS**<br>
# This function can prepare metadata in a number of ways :
# - It can infer an indication from the file names if indication_in_name is TRUE, using str_split_i
# - It can use a threshold file to add the thresholds to the metadata
# - It can also take any additionnal manual metadata from the additional_metadata list

# %%
prepare_metadata <- function(out_dir,
                             use_threshold,
                             threshold_dir = NULL,
                             keep_lower = TRUE,
                             keep_upper = TRUE,
                             intens_dir,
                             file_extension = ".tsv",
                             subset = "",
                             panel = "",
                             no_threshold = "",
                             tag = "",
                             exclude = NULL,
                             indication_in_name = TRUE,
                             indication_position = 2,
                             separator = "_",
                             indication = NULL,
                             additional_metadata = NULL,
                             store_csv = TRUE,
                             use_roi = FALSE,
                             roi_in_name = TRUE,
                             roi_position = 1,
                             sample_id_position = 3,
                             num_threads = 1) {
  #Listing files
  if (use_threshold) {
    files <- list.files(path = threshold_dir,
                        pattern = subset,
                        full.names = TRUE)

    files_name <- list.files(path = threshold_dir, pattern = subset) |>
      stringr::str_split_i(file_extension, 1)

  } else {
    files <- list.files(path = intens_dir,
                        pattern = subset,
                        full.names = TRUE)
    # explicitly sort the files alphabetically
    files_name <- list.files(path = intens_dir, pattern = subset) |>
      stringr::str_split_i(file_extension, 1)
  }
  #Handling incorrect length of indication
  if (!indication_in_name &&
        length(indication) != length(files)) {
    cli::cli_abort(c(
      "Incorrect number of value for manual indication.",
      "x" = "Either set indication_in_name to TRUE or provide an indication for each sample.", # nolint: line_length_linter.
      "i" = "Number of provided indications : {length(indication)}",
      "i" = "Number of provided files : {length(files)}"
    ))
  }
  #Preparing sample id

  sample_id <- purrr::map_chr(
    files_name,
    function(x) stringr::str_split_i(x, separator, sample_id_position)
  )

  #preparing automatic indication :
  if (indication_in_name) {
    indication <- purrr::map_chr(
      files_name,
      function(x) stringr::str_split_i(x, separator, indication_position)
    )
  }

  #prepare roi name :
  if (use_roi && roi_in_name) {
    roi <- purrr::map_chr(
      files_name,
      function(x) stringr::str_split_i(x, separator, roi_position)
    )
  }

  #preparing addition metadata dataframe :
  if (length(additional_metadata > 0)) {
    add_df <- data.frame("sample_id" = sample_id)
    for (i in seq_along(additional_metadata)) {
      add_df <- cbind(add_df, additional_metadata[[i]])
    }
  }

  #Using prepare_threshold to get the thr_dataframe
  if (use_threshold && !is.null(threshold_dir)) {
    #initialize thr_df
    thr_df <- NULL

    thr_df <- purrr::reduce(
      threshold_dir,
      function(thr_df, x) {
        new_df <- prepare_threshold(
          threshold_dir = x,
          subset = subset,
          tag = tag,
          exclude = exclude,
          panel = panel,
          no_threshold = no_threshold,
          threshold_extension = file_extension,
          keep_lower = keep_lower,
          keep_upper = keep_upper,
          num_threads = num_threads
        )
        thr_df <- rbind(thr_df, new_df)

        return(thr_df)
      },
      .init = thr_df
    )
  }

  #Merge dataframes :
  final_df <- cbind(sample_id, indication)
  if (use_roi) final_df <- cbind(final_df, roi)
  if (exists("add_df")) final_df <- cbind(final_df, add_df)
  if (exists("thr_df")) final_df <- cbind(final_df, thr_df)
  #Write it on the disk :
  if (store_csv) {
    readr::write_csv(
      final_df,
      file = paste0(out_dir, "metadata.csv"),
      quote = "none",
      num_threads = num_threads
    )

    invisible(final_df)
  } else {
    return(final_df)
  }
}

# %% [markdown]
# ## Prepare_threshold
# 
# This functions takes the files within a threshold folder as an input, and outputs a table with a row for each sample and a column for  the lower threshold of each marker. <br>
# This function is called within `prepare_metadata`
# 
# **ARGUMENTS**
# - `threshold_dir ` : path to folder where the threshold files are located
# - `subset` : character string used to identify, with regular expression, the subset of the files to consider. Setting `subset = ""` will load all files.
# - `panel` : character vector used to identify, with regular expression, the markers from the column names.
# - `threshold_extension` : charcater string, indicating the extension the threshold files. Only ".csv" and ".tsv" are supported
# - `keep_lower` and `keep_upper` : Logical values, wheter to keep the lower and upper threshold in the metadata respectively
# 
# **VALUE** <br>
# Returns a dataframe that contains the positivity threshold for all markers (column) and individual (rows).

# %%
prepare_threshold <- function(
  threshold_dir,
  subset,
  tag,
  exclude = NULL,
  panel,
  no_threshold,
  threshold_extension,
  keep_lower = TRUE,
  keep_upper = TRUE,
  num_threads = 1
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))

  thr_files <- list.files(path = threshold_dir,
                          pattern = subset)

  thr_df <- NULL
  #will select column
  used <- tag[!tag %in% no_threshold] %>%
    .[!grepl(exclude, .)]
  #For colnames att the end
  used_panel <- panel[!tag %in% no_threshold] %>%
    .[!grepl(exclude, .)]
  #will select rows
  row_to_keep <- rlang::chr()
  if (keep_lower) row_to_keep <- c(row_to_keep, "Lower Threshold")
  if (keep_upper) row_to_keep <- c(row_to_keep, "Upper Threshold")

  thr_df <- purrr::reduce(
    thr_files,
    function(thr_df, x) {
      df <- load_file(
        file = paste0(threshold_dir, x),
        extension = threshold_extension,
        num_threads = num_threads
      ) %>%
        dplyr::filter(rowName %in% row_to_keep) %>%
        select(contains(used)) %>%
        select(!contains(exclude))

      #transform to wide if 2 rows are kept
      if (nrow(df) == 2) df <- cbind(df[1, ], df[2, ])

      thr_df <- rbind(thr_df, df)
    },
    .init = thr_df
  )

  #fix names :
  rownames(thr_df) <- thr_files %>%
    str_split_i(threshold_extension, 1)

  lower_colnames <- paste0("lower_threshold_", used_panel)
  upper_colnames <- paste0("upper_threshold_", used_panel)
  #Add the the colnames that correspond to the one used
  total_colnames <- NULL
  if (keep_lower) total_colnames <- c(total_colnames, lower_colnames)
  if (keep_upper) total_colnames <- c(total_colnames, upper_colnames)
  #Set colnames and return thr_df
  colnames(thr_df) <- total_colnames
  return(thr_df)
}

# %% [markdown]
# ### load_file
# 
# Simple helper used to call read_csv or read_tsv accordingly

# %%
load_file <- function(file, extension = c(".tsv", ".csv"), num_threads = 1) {
  extension <- match.arg(extension)

  if (extension == ".tsv") {
    readr:::read_tsv(file, show_col_types = FALSE, num_threads = num_threads)
  } else if (extension == ".csv") {
    readr::read_csv(file, show_col_types = FALSE, num_threads = num_threads)
  }
}

# %% [markdown]
# # Data preprocessing

# %% [markdown]
# ## create_threshold_assay
# 
# This function is used to create new assay transformed by a threshold. 
# Thresholds are expected to be recorded within meta, and set individually for each marker and each sample_id
# 
# **ARGUMENTS :**
# - `spe` : SpatialExperiment object or SingleCellExperiment Object
# - `meta` : dataframe containing the metadata including the thresholds
# - `panel` : character vector, names of the markers, should match the rownames of spe. The names of the collumn of meta are expected to contain the corresponding string
# - `name` and `logical_name` : character values, names given to the new thresholded and logical assay respectively
# - `lower_threshold` and `upper_threshold` : logicles, whether to use the upper or lower threhsold
# - `no_threshold` : character vector, names of the markers for which no threshold was set
# - `logical` : logical (sic) value, whether to also record an assay transformed as 0 (outside threshold range) or 1 (inside)
# - `discard_no_threshold` : logical, whether to set every unthresholded markers to 0 ONLY FOR THE LOGICAL ASSAY
# 
# **VALUE** <br>
# Return the spe object with the new assay
# 
# **DETAILS** <br>
# Assay will be available at`assay(spe, name)` and `assay(spe, logical_name)`

# %%
create_threshold_assay <- function(
  spe,
  param_name = "preprocessing",
  meta = default_parameter(spe, param_name, "meta"),
  panel = default_parameter(spe, param_name, "panel"),
  base_assay = "counts",
  threshold_name = "threshold_counts",
  lower_threshold = FALSE,
  upper_threshold = FALSE,
  no_threshold = NULL,
  logical = FALSE,
  logical_name = "threshold_logical",
  discard_no_threshold = TRUE
) {
  env <- environment()
  param_list <- get_called_param(env)
  spe <- spe |>
    address_parameters(
      param_name = param_name,
      parameters = param_list,
      envir = env
    )

  #Initialise assay
  used_threshold <- panel[!panel %in% no_threshold]
  samples <- unique(spe$sample_id)

  new_assay <- NULL
  new_assay <- purrr::reduce(
    samples,
    function(new_assay, x) {
      subset_assay <- apply_threshold(
        spe,
        sample = x,
        panel = used_threshold,
        meta = meta,
        assay_name = base_assay,
        lower_threshold = lower_threshold,
        upper_threshold = upper_threshold
      )
      #Still need to deal with the rows that do not have a threshold : 
      kept_assay <- SummarizedExperiment::assay(spe, base_assay) %>%
        .[rownames(.) %in% no_threshold, spe$sample_id == x]
      #combines then reorders the assay
      subset_assay <- rbind(subset_assay, kept_assay) %>%
        .[rownames(spe), ]
      new_assay <- cbind(new_assay, subset_assay)
    },
    .init = new_assay
  )
  
  #Integrate the completed new_assay within the SPE object.
  SummarizedExperiment::assay(spe, threshold_name) <- new_assay
  if (logical) {
    new_assay[new_assay != 0] <- 1
    if (discard_no_threshold) {
      new_assay[rownames(new_assay) %in% no_threshold, ] <- 0
    }
    SummarizedExperiment::assay(spe, logical_name) <- new_assay
  }

  return(spe)
}

# %% [markdown]
# ### apply_threshold
# 
# Tocomment
# returns an assay modified by the threshold

# %%
apply_threshold <- function( #returns an assay modified by the threshold for that sample
  spe,
  sample,
  panel, #is expected to be the panel with threshold
  meta, #tibble containing the metadata, including the thresholds
  assay_name = "counts",
  lower_threshold = FALSE,
  upper_threshold = FALSE
) {
  subset_assay <- spe[rownames(spe) %in% panel, spe$sample_id == sample] %>%
    SummarizedExperiment::assay(assay_name)

  if (lower_threshold) {
    #Extract the threshold
    lower_meta <- meta %>%
      dplyr::filter(sample_id == sample) %>%
      dplyr::select(
        dplyr::contains("lower_threshold")
      )

    #matching name format with the rownames of the assay
    colnames(lower_meta) <- str_split_i(
      colnames(lower_meta),
      pattern = "_",
      i = 3
    )
    #ensures same order and create a vector for map2
    lower_meta <- lower_meta[, rownames(subset_assay)] %>%
      as.integer()

    old_names <- rownames(subset_assay)

    subset_assay <- purrr::map2(
      split(subset_assay, row(subset_assay)),
      lower_meta,
      function(x, y) {
        x[x < y] <- 0
        x
      }#returns list of rows
    ) %>%
      do.call(rbind, .) #recreates the matrix
    #and fix the names
    rownames(subset_assay) <- old_names
  }
  if (upper_threshold) {
    upper_meta <- meta %>%
      filter(sample_id == sample) %>%
      select(contains("upper_threshold"))

    #matching name format with the rownames of the assay
    colnames(upper_meta) <- str_split_i(
      colnames(upper_meta),
      pattern = "_",
      i = 3
    )
    #ensures same order and create a vector for map2
    upper_meta <- upper_meta[, rownames(subset_assay)] %>%
      as.integer()

    old_names <- rownames(subset_assay)

    subset_assay <- purrr::map2(
      split(subset_assay, row(subset_assay)),
      upper_meta,
      function(x, y) {
        x[x > y] <- 0
        x
      }#returns list of rows
    ) %>%
      do.call(rbind, .) #recreates the matrix
    #and fix the names
    rownames(subset_assay) <- old_names
  }

  return(subset_assay)
}

# %% [markdown]
# ## transform_assay
# 
# This function takes a spatialExperiment object, modify an assay with a transformation and can generate quality control plot at the same time <br>
# 
# **ARGUMENTS**
# - `spe` : An object inheriting SummarizedExperiment
# - `param_name` : string indicating where parameters are stored. See parameters for details.
# - `assay_name` : a matrix like object, like the assay slot of SingleCellExperiment objects
# - `transform_method` : character string, the method used for transformation.
# - `transform_parameter` : numeric value, the parameter used for transformation. The use depends on the method.
# - `plot_transformation` : logical, whether to also produce associated plots (ridge plots and signal noise ratio plots)
# - `markers_to_plot` : character vector, subset of rownames of spe to plot
# - `downsample_plot` : numerical value between 0 and 1, the proportion of cells sampled for the plots
# - `seed` : : numerical value, for RNG control
# - `patient_id` : character string, name of the colData that contains the patient or samples identification
# - `store_plot` : logical, wether to store plot as a file, or in the metadata
# - `separate` : logical, wether to store the ridgeplots as a single file or separate plots
# - `device_type` : string, which graphical device to use,
# - `base_height`, `base_width` : numerical values in inches, dimension of the plots
# - `plot_dir` : path to directory where plots are stored
# - `quiet` : logical, wether to hide additional messages
# - `...` : passed to the graphical device, for additional control over plots
# 
# **VALUE**<br>
# Return the spe with the modified assay in assay(spe, transformed_name). <br>
# If plot_transformation is TRUE, will also generate ridgeplots and signal to noise ratio plots. <br>
# If store_plot is true, plots will be located in plot_dir/transformed_plot/{transformed_method}_{transformed_parameters} <br>
# If store_plot is FALSE, plots will be located the metadata of spe in the "transformation" category
# 
# **DETAILS**<br>
# *Formulas:*
# - asinh : `y = asinh(x/parameter)` with asinh(x) = log(x + sqrt(x^2 + 1))
# - asinh_minus : `y = asinh((x - parameter)/parameter)`
# - log : `y = log(x + parameter)`
# - roots : `y = x^(1/parameter)`
# - standard : `y = (x - mean(x))/sd(x)` 
# 
# Plots within the metadata can be accessed with :
# get_metadata(spe, category = "transformation", name = paste0(transform_method, "_", transform_parameter))

# %%
transform_assay <- function(
  spe,
  param_name = "preprocessing",
  assay_name = default_parameter(spe, param_name, "assay_name"),
  transform_method = c("asinh", "log", "standard", "roots", "asinh_minus"),
  transform_parameter = 1,
  transformed_name = "transformed",
  plot_transformation = TRUE,
  markers_to_plot = rownames(spe),
  downsample_plot = 0.1,
  seed = 42,
  patient_id = "patient_id",
  store_plot = FALSE,
  separate = FALSE,
  device_type = c("pdf", "postscript", "svg",
                  "bitmap", "png", "jpeg", "bmp", "tiff"),
  base_height = 5,
  base_width = 5,
  plot_dir = "./",
  quiet = FALSE,
  called = FALSE,
  ...
) {
  if(!called) {
    env <- environment()
    param_list <- get_called_param(env)
    spe <- spe |>
      address_parameters(
        param_name = param_name,
        parameters = param_list,
        envir = env
      )
  }

  #Checking arguments :
  transform_method <- match.arg(transform_method)
  assay <- SummarizedExperiment::assay(spe, assay_name)

  SummarizedExperiment::assay(spe, transformed_name) <- switch(
    transform_method,
    "asinh" = asinh(assay / transform_parameter),
    "asinh_minus" = asinh((assay - transform_parameter) / transform_parameter),
    "log" = log(assay + transform_parameter),
    "standard" = (assay - mean(assay)) / sd(assay),
    "roots" = assay^(1 / transform_parameter)
  )

  if (plot_transformation) {
    if (store_plot) {
      transform_dir <- paste0(plot_dir, "transform_plot/")
      called_dir <- paste0(
        transformed_dir, transform_method, "_", transform_parameter, "/"
      )
      if (!dir.exists(plot_dir)) dir.create(plot_dir)
      if (!dir.exists(transform_dir)) dir.create(transform_dir)
      if (!dir.exists(called_dir)) dir.create(called_dir)
      metadata_name <- NULL
    } else {
      called_dir <- NULL
      metadata_name <- paste0(
        transform_method, "_", transform_parameter
      )
    }

    spe <- plot_transformation(
      spe,
      called = TRUE,
      assay_name = assay_name,
      transformed_name = transformed_name,
      markers_to_plot = markers_to_plot,
      downsample_plot = downsample_plot,
      seed = seed,
      patient_id = patient_id,
      metadata_name = metadata_name,
      store_plot = store_plot,
      separate = separate,
      device_type = device_type,
      plot_dir = called_dir,
      quiet = quiet,
      base_height = base_height,
      base_width = base_width,
      ...
    )
  }
  return(spe)
}

# %% [markdown]
# ### plot_transformation
# 
# Helper function to create plots related to the transformation (ridgeplots and snr plots), create on a subset of spe. <br>
# It also address them at the same time based on store_plot :
# - if TRUE, will store the plots in the plot_dir, in a format decided by device_type, using the organize_and_store function
# - if FALSE, it will store the plots in the metadata, under the "transformation" category

# %%
plot_transformation <- function(
  spe,
  param_name = "preprocessing",
  assay_name = "counts",
  transformed_name = "transformed",
  markers_to_plot = rownames(spe),
  downsample_plot = 0.1,
  seed = 42,
  patient_id = patient_id,
  metadata_name = transformed_name,
  store_plot = FALSE,
  separate = FALSE,
  device_type = c("pdf", "postscript", "svg",
                  "bitmap", "png", "jpeg", "bmp", "tiff"),
  base_height = 5,
  base_width = 5,
  plot_dir = "./",
  quiet = FALSE,
  called = FALSE,
  ...
) {
  if (!called) {
    env <- environment()
    param_list <- get_called_param(
      env,
      excluded = c("spe", "param_name", "register_param", "register_env")
    )
    spe <- spe |>
      address_parameters(
        param_name = param_name,
        parameters = param_list,
        envir = env
      )
  }
  # Prepare subset
  set.seed(seed)
  n_cells <- ncol(spe)
  # Prepare the indices that will serve for downsampling
  cell_indices <- n_cells |>
    seq_len() |>
    sample(size = downsample_plot * n_cells)

  subset <- spe[markers_to_plot, cell_indices]

  # Create plots
  ridge_plots <- plot_expression(
    subset,
    called = TRUE,
    assay_name = transformed_name,
    markers_to_plot = markers_to_plot,
    plot_type = "ridgeplot",
    group.by = patient_id,
    quiet = quiet,
    store_plot = FALSE, # addressed separately
    register_param = FALSE
  ) |>
    setNames(markers_to_plot)

  snr_plot <- plot_snr_cells(
    subset,
    trans_assay = transformed_name,
    count_assay = assay_name,
    seed = seed
  )

  # Adress plots depending on store_plot :
  if (store_plot) {
    if (!quiet) message("Storing plots")
    # Plotting ridges
    organize_and_store(
      ridge_plots,
      plot_name = names(ridge_plots),
      plot_dir = plot_dir,
      separate = separate,
      device_type = device_type,
      base_height = base_height,
      base_width = base_width,
      ...
    )

    organize_and_store(
      list(snr_plot),
      plot_name = list("snr_plot"),
      plot_dir = plot_dir,
      separate = separate,
      device_type = device_type,
      base_height = base_height,
      base_width = base_width,
      ...
    )
  } else {
    spe <- spe |>
      store_multiple_metadata(
        ridge_plots,
        category = "transformation",
        name = metadata_name,
        type = names(ridge_plots)
      ) |>
      store_metadata(
        snr_plot,
        category = "transformation",
        name = metadata_name,
        type = "snr_plot"
      )
  }
  
  invisible(spe)
}

# %% [markdown]
# ## preprocess_data

# %%
preprocess_data <- function(
  spe,
  param_name = "preprocessing",
  assay_name = "counts",
  transformed_name = "transformed",
  integration = "none",
  group_by = "sample_id",
  downsampling = 1,
  plot_transformation = FALSE,
  use_channel = rownames(spe),
  seed = 42,
  quiet = FALSE
) {
  env <- environment()
  param_list <- get_called_param(env)
  spe <- spe |>
    address_parameters(
      param_name = param_name,
      parameters = param_list,
      envir = env
    )

  if (!quiet) message("Subsampling")
  set.seed(seed)
  n_cells <- ncol(spe)
  # Prepare the indices that will serve for downsampling
  cell_indices <- n_cells %>%
    seq_len() %>%
    sample(size = downsampling * n_cells)

  spe <- spe[, cell_indices]

  if (!quiet) message("Transforming assay")

  spe <- transform_assay(spe, param_name)

  if (!quiet) message("Building PCA")
  spe <- scater::runPCA(
    spe,
    name = "PCA",
    subset_row = use_channel,
    exprs_values = transformed_name,
    BSPARAM = BiocSingular:::ExactParam()
  )

  if (!quiet) message("Building UMAP")
  spe <- scater::runUMAP(
    spe,
    name = "UMAP",
    ncomponents = 2,
    subset_row = use_channel,
    exprs_values = transformed_name
  )

  integration <- tolower(integration)

  if ("harmony" %in% integration) {

    if (!quiet) message("Integrating with Harmony")
    spe <- harmony::RunHarmony(
      spe,
      group.by.vars = group_by,
      assay.type = transformed_name
    )

    if (!quiet) message("Building UMAP from Harmony")
    spe <- scater::runUMAP(
      spe,
      dimred = "HARMONY",
      name = "UMAP_harmony_corrected"
    )
  }

  if ("fastmnn" %in% integration) {
    if (!quiet) message("Integrating with fastMNN")
    out <- batchelor::fastMNN(
      spe,
      batch = spe[[group_by]],
      auto.merge = TRUE,
      subset.row = use_channel,
      assay.type = transformed_name,
      BSPARAM = BiocSingular::ExactParam()
    )

    SingleCellExperiment::reducedDim(spe, "fastMNN") <-
      SingleCellExperiment::reducedDim(out, "corrected")

    if(!quiet) message("Building UMAP from fastMNN")
    spe <- scater::runUMAP(
      spe,
      dimred = "fastMNN",
      name = "UMAP_fastmnn_corrected"
    )
  }

  return(spe)
}

# %% [markdown]
# # Clustering

# %% [markdown]
# ## perform_clustering
# 
# This function takes a spatialExperiment object (or a SingleCellExperiment object) and returns it with one clustering done.
# 
# **ARGUMENTS**
# 
# - `spe`: a spatialExperiment object (or singleCellExperiment)
# - `param_name` : name of parameters to use, see parameters for details.
# - `assay_name` : character string indicating the assay containing the transformed expression values
# - `dim_red_for_clusters` : character value indicating the name of the reducedDim of spe to use for clustering
# - `cluster_method` : character value indicating which clustering function to use.<br>
# Currently supported "FlowSOM", "louvain" and "BLUSPARAM" <br>
# Using "BLUSPARAM" also requires specification of the list_BLUSPARAM argument
# - `flowsom_parameter` : integer value, corresponding to max_K the maximum number of metacluster to consider
# - `xdim` and `ydim` : integer values, corresponding to the dimension of the initial grid for the Self Organizing Map
# - `louvain_parameter` : integer vector, corresponding to the k number of nearest neighbor to use for intiating the graph
# - `weight_type_parameter` : character vector, indicating the function used for weight to give to edges during graph initiation. <br>
# See makeSNNGraph from bluster package for more information.
# - `metacluster` : logical, wether to perform metaclustering using consensusClusterPlus. Always happen for FlowSOM (as part of the standard workflow), but can be toggled for the other clustering_methods
# - `seed` : integer value, the seed used for random number generation
# - `BLUSPARAM` : one BlusterParam objects, passed to clusterRows if "BLUSPARAM" is in `cluster_method`
# - `BPPARAM` : BiocParallelParam object, to control parallel processing
# - `numthreads` : integer, number of threads used for parallel processing. Is overtakeng by BPPARAM if nworker is specified.
# - `quiet` : logical value, wheter to silence the verbose.
# - `called` : logical, for developper, see parameters.
# 
# **VALUE**<br>
# Return the same spatialExperiment object with all clustering done.<br>
# 
# **DETAILS**<br>
# The clustering will be stored in the colData slot, which can be accessed via spe$name_your_cluster. <br>
# For flowSOM, only the clusters are stored in the colData slot, the rest of the informations is store in the metadata slot. <br>
# 
# The information about every clustering performed is stored in the metadata slot, accessible via metadata(spe)$clustering_methods <br>
# It is recorded in a tibble with 5 columns :
# - "name" which match the name of the colData
# - "fun" containing the function used
# - "dim_red" containing the dimension reduction used
# - "metaclsutering" : logical, wether metaclustering has been done
# - "parameters" containing a list of all the parameters used. <br>

# %%
perform_clustering <- function(
  spe,
  param_name = "cluster",
  assay_name = "transformed",
  dim_red_for_clusters = "none",
  cluster_method = c("flowsom", "louvain", "blusparam"),
  xdim = 10,
  ydim = 10,
  return_fsom = TRUE,
  louvain_parameter = 10,
  weight_type_parameter = c("jaccard", "rank"),
  metacluster = FALSE,
  max_k = 10,
  seed = 42,
  BLUSPARAM = list(),
  quiet = FALSE,
  num_thread = 1,
  BPPARAM = BiocParallel::SerialParam(),
  metadata_method = "clustering_methods",
  called = FALSE,
  ...
) {
  # Adress parameters
  if (!called) {
    env <- environment()
    param_list <- get_called_param(env)
    spe <- spe |>
      address_parameters(
        param_name = param_name,
        parameters = param_list,
        envir = env
      )
  }


  validate_perform_cluster(
    spe,
    assay_name = assay_name,
    dim_red_for_clusters = dim_red_for_clusters,
    cluster_method = cluster_method,
    weight_type_parameter = weight_type_parameter,
    BLUSPARAM = BLUSPARAM,
    BPPARAM = BPPARAM
  )

  if (cluster_method == "flowsom") {
    spe <- perform_flowsom(
      spe,
      param_name = param_name,
      assay_name = assay_name,
      max_k = max_k,
      flowsom_red_dim = dim_red_for_clusters,
      metadata_method = metadata_method,
      quiet = quiet,
      return_fsom = TRUE,
      xdim = xdim,
      ydim = ydim,
      called = TRUE,
      ...
    )
  } else if (cluster_method == "louvain") {
    spe <- perform_louvain(
      spe,
      param_name = param_name,
      dim_red_for_clusters = dim_red_for_clusters,
      weight_type_parameter = weight_type_parameter,
      louvain_parameter = louvain_parameter,
      metadata_method = metadata_method,
      metacluster_louvain = metacluster,
      num_thread = num_thread,
      assay_name = assay_name,
      max_k = max_k,
      seed = seed,
      quiet = quiet,
      called = TRUE,
      ...
    )
  } else if (cluster_method == "blusparam") {
    spe <- perform_blusparam(
      spe,
      param_name = param_name,
      dim_red_for_clusters = dim_red_for_clusters,
      list_blusparam = BLUSPARAM,
      name_blusparam = name,
      metadata_method = "clustering_methods",
      metacluster_blusparam = metacluster,
      num_thread = num_thread,
      assay_name = assay_name,
      max_k = max_k,
      seed = seed,
      quiet = quiet,
      called = TRUE,
      ...
    )
  }

  return(spe)
}

# %% [markdown]
# ## perform_flowsom 
# 
# This function takes a spatialExperiment (or SingleCellExperiment) object and returns it after performing flowSOM and ConsensusClusterPlus
# 
# **ARGUMENTS** <br>
# - `spe`: a spatialExperiment object (or SingleCellExperiment)
# - `assay_name` : character string, indicating which assay contains the transformed expression values
# - `max_k` : character value, indicating the maximum number of metaclusters to consider <br>
# MetaclusterPlus will perform all metaclustering from 2 metaclusters to max_k. <br>
# - `flowsom_red_dim` : character string (non vectorized), indicating which dimension reduction to use for flowsom. <br>
# Currently supported : "none", "harmony", "fastmnn"
# - `metadata_method` character string indicating the name under which the method will be recorded inside the metadata of the spe. <br>
# Setting it to NULL will not record the clustering in the metadata.
# -`quiet` logical value, wheter to silence the verbose.
# - `return_fsom` : logical value, wheter to return the FlowSOM object inside the metadata of the spe
# - `xdim` and `ydim` : integer values, corresponding to the dimension of the initial grid for the Self Organizing Map
# - `...` : any other argument to be passed to BuildMST, which itselfs passes it to flowSOM::SOM <br>
# Some of the relevant arguments are :
# - `rlen` : integer value, number of times to loop over training data
# - `mst` : integer value, number of times to build a Minimal Spanning Tree
# - `alpha` : numeric value, start and end learning rate
# - `radius` : numeric value, start and end radius,
# - `importance` : array with numeric values. Each parameters will be scaled according to importance
# 
# **VALUE**<br>
# The spatialExperiment (or SingleCellExperiment) with cluster based on flowSOM and metaclusters based on ConsensusClusterPlus
# 
# **DETAILS** <br>
# This function is heavily inspired by and reuse some parts of the CATALYST::cluster function. <br>
# Compared to CATALYST::cluster, this function proposes two main features :
# - It can use dimension reduction, such as those created by batch effect correction algorithm,
# - It can return the flowSOM object, inside the metadata slot of the spe, to be used with the plotting function from the flowSOM package<br>
# 
# The flowSOM results are returned thusly :
# - the clusters are stored in the colData slot, named `flowsom_{flowsom_red_dim}_{xdim}x{ydim}`
# - the metacluster codes are stored in `metadata(spe)$cluster_codes`
# - the flowsom codes (position on the grid) are stored in `metadata(spe)$SOM_codes`
# - the delta area plot is stored in `metadata(spe)$delta_area`
# - the flowSOM object is stored in `metadata(spe)$fsom`
# 
# > All metadata are stored with the name followed by _{flowsom_red_dim}_{xdim}x{ydim}, which we omit here for clarity. <br>

# %%
perform_flowsom <- function(spe,
                            param_name = "cluster",
                            assay_name = "counts",
                            max_k = 20,
                            flowsom_red_dim = c("none", "harmony", "fastmnn"),
                            metadata_method = "clustering_methods",
                            quiet = FALSE,
                            return_fsom = TRUE,
                            xdim = 10,
                            ydim = 10,
                            BPPARAM = BiocParallel::SerialParam(),
                            seed = 42,
                            called = FALSE,
                            ...) {
  ### CHECKING PARAMETER ###############################################
  if (!called) {
    env <- environment()
    param_list <- get_called_param(
      env,
      excluded = c("spe", "param_name", "called")
    )
    spe <- spe |>
      address_parameters(
        param_name = param_name,
        parameters = param_list,
        envir = env
      )
  }

  if (is.null(flowsom_red_dim)) flowsom_red_dim <- "none"

  flowsom_red_dim <- tolower(flowsom_red_dim)
  flowsom_red_dim <- match.arg(flowsom_red_dim)

  if (flowsom_red_dim == "harmony" &&
        !"HARMONY" %in% SingleCellExperiment::reducedDimNames(spe)) {
    abort_not_found(
      "flowsom_red_dim",
      valid = "Dimension reduction for flowSOM",
      choices = SingleCellExperiment::reducedDimNames(spe),
      not = flowsom_red_dim
    )
  } else if (flowsom_red_dim == "fastmnn" &&
               !"fastMNN" %in% SingleCellExperiment::reducedDimNames(spe)) {
    abort_not_found(
      "flowsom_red_dim",
      valid = "Dimension reduction for flowSOM",
      choices = SingleCellExperiment::reducedDimNames(spe),
      not = flowsom_red_dim
    )
  }
  ### BUILDING FLOWSOM TREE ########################################
  if (is.null(flowsom_red_dim)) flowsom_red_dim <- "none"

  transposed_assay <- switch(
    flowsom_red_dim,
    "none" =  t(assay(spe, assay_name)[rowData(spe)$use_channel, ]), #nolint
    "harmony" = reducedDim(spe, "HARMONY"),
    "fastmnn" = reducedDim(spe, "fastMNN"),
    default = t(assay(spe, assay_name)[rowData(spe)$use_channel, ]) #nolint
  )
  #Need to provide colname to fastMNN for ReadInput
  if (flowsom_red_dim == "fastmnn"){
    n_col <- ncol(transposed_assay)
    col_names <- paste0("MNN_", seq_len(n_col))
    colnames(transposed_assay) <- col_names
  }
  
  #Prepare cluster names
  cluster_name <- paste0("flowsom_", flowsom_red_dim, "_", xdim, "x", ydim) %>%
    tolower()

  if (!quiet) message("Building FlowSOM tree : ", cluster_name, ".")

  #Create the flowSOM object step by step :
  fsom <- FlowSOM::ReadInput(transposed_assay)
  fsom <- FlowSOM::BuildSOM(fsom,
                   silent = TRUE,
                   xdim = xdim,
                   ydim = ydim,
                   ...)
  fsom <- FlowSOM::BuildMST(fsom,  silent = TRUE, tSNE = FALSE)
  #Integrate the clusters into colData(spe)
  spe[[cluster_name]] <- fsom$map$mapping[, 1] %>%
    as.factor()

  param <- list(
    name = cluster_name,
    type = "flowsom",
    reducedDim = flowsom_red_dim,
    parameters = list(
      "max_k" = max_k,
      "xdim" = xdim,
      "ydim" = ydim
    ),
    metaclustering = FALSE
  )
  spe <- spe %>%
    store_metadata(
      param,
      category = metadata_method,
      name = cluster_name
    )

  #incorporate the fsom object into the metadata
  if (return_fsom) {
    spe <- spe %>%
      store_metadata(
        fsom,
        category = "clustering_data",
        name = cluster_name,
        type = "fsom"
      )
  }

  ### METACLUSTERING #################################################################
  dim_red <- switch(
    flowsom_red_dim,
    "harmony" = "HARMONY",
    "fastmnn" = "fastMNN",
    "none" = ,
    .default = NULL
  )
  spe <- build_metacluster(
    spe,
    param_name = paste0(param_name, "_metacluster"),
    cluster_name = cluster_name,
    cluster_type = "flowsom",
    dim_red = dim_red,
    assay_name = assay_name,
    clustering_metadata = metadata_method,
    max_k = max_k,
    distance = "euclidean",
    BPPARAM = BPPARAM,
    quiet = quiet,
    seed = seed,
    called = TRUE,
    ...
  )

  return(spe)
}

# %% [markdown]
# ## perform_louvain
# 
# Designed to be called inside a multicluster, this functi

# %%
perform_louvain <- function(
  spe,
  param_name = "cluster",
  dim_red_for_clusters = "none",
  weight_type_parameter = "jaccard",
  louvain_parameter = 10,
  metadata_method = "clustering_methods",
  metacluster_louvain = FALSE,
  assay_name = "counts",
  max_k = 3,
  BPPARAM = BiocParallel::SerialParam(),
  num_thread = BiocParallel::bpnworkers(BPPARAM),
  quiet = TRUE,
  called = FALSE,
  ...
) {
  if (!called) {
    env <- environment()
    param_list <- get_called_param(
      env,
      excluded = c("spe", "param_name", "called")
    )
    spe <- spe |>
      address_parameters(
        param_name = param_name,
        parameters = param_list,
        envir = env
      )
  }

  cluster_name <- paste0(
    "louvain_",
    dim_red_for_clusters, "_",
    weight_type_parameter, "_",
    louvain_parameter
  ) %>%
    tolower()

  if (dim_red_for_clusters == "none") dim_red_for_clusters <- NULL

  if (!quiet) message("Running ", cluster_name, "\n")
  #Performs louvain clustering
  spe[[cluster_name]] <- scran::clusterCells(
    spe,
    assay.type = assay_name,
    use.dimred = dim_red_for_clusters,
    BLUSPARAM = bluster::SNNGraphParam(
      k = louvain_parameter,
      cluster.fun = "louvain",
      type = weight_type_parameter,
      num.threads = num_thread
    )
  ) %>%
    as.factor()
  #Add the name to the clustering methods
  param <- list(
    name = cluster_name,
    type = "louvain",
    reducedDim = dim_red_for_clusters,
    parameters = list(
      "weight_type_parameter" = weight_type_parameter,
      "k_neighbors" = louvain_parameter
    ),
    metaclustering = FALSE
  )

  spe <- spe %>%
    store_metadata(
      param,
      category = metadata_method,
      name = cluster_name
    )

  if (metacluster_louvain) {
    if (dim_red_for_clusters == "none") {
      dim_red <-  NULL
    } else {
      dim_red <- dim_red_for_clusters
    }

    spe <- build_metacluster(
      spe,
      param_name = paste0(param_name, "_metacluster"),
      cluster_name = cluster_name,
      cluster_type = "louvain",
      dim_red = dim_red,
      assay_name = assay_name,
      clustering_metadata = metadata_method,
      max_k = max_k,
      distance = "euclidean",
      BPPARAM = BPPARAM,
      quiet = quiet,
      called = TRUE,
      ...
    )

  }
  return(spe)
}

# %% [markdown]
# ## perform_blusparam
# 
# NB : num_thread will affect metaclustering but not clustering, as the thread usage is defined inside each BLUSPARAM object. If you want to make use of multi-threading, make sure to modify the BLUSPARAM accordingly.

# %%
perform_blusparam <- function(
  spe,
  param_name = "cluster",
  dim_red_for_clusters = "none",
  BLUSPARAM = NULL,
  name_blusparam = NULL,
  metadata_method = "clustering_methods",
  metacluster_blusparam = FALSE,
  BPPARAM = BiocParallel::SerialParam(),
  num_thread = BiocParallel::bpnworkers(BPPARAM),
  assay_name = "counts",
  max_k = 3,
  quiet = TRUE,
  called = FALSE,
  ...
) {
  if (!called) {
    env <- environment()
    param_list <- get_called_param(
      env,
      excluded = c("spe", "param_name", "called")
    )
    spe <- spe |>
      address_parameters(
        param_name = param_name,
        parameters = param_list,
        envir = env
      )
  }

  if (is.null(list_blusparam)){ #Short circuit if called without blusparam
    return(spe)
  }
  
  # prepare arguments
  cluster_name <- paste0(
    "custom_",
    dim_red_for_clusters, "_",
    name_blusparam
  ) %>%
    tolower()
  
  if (dim_red_for_clusters == "none") dim_red_for_clusters <- NULL

  if (!quiet) message("Running", cluster_name, "\n")
  spe[[cluster_name]] <- scran::clusterCells(
    spe,
    assay.type = assay_name,
    use.dimred = dim_red_for_clusters,
    BLUSPARAM = BLUSPARAM
  ) %>%
    as.factor()

  #Record metadata
  param <- list(
    name = cluster_name,
    type = "custom",
    reducedDim = dim_red_for_clusters,
    parameters = list(
      "BLUSPARAM" = BLUSPARAM
    ),
    metaclustering = FALSE
  )
  spe <- spe %>%
    store_metadata(
      param,
      category = metadata_method,
      name = cluster_name
    )

  if (metacluster_blusparam) {
    if (dim_red_for_clusters == "none") {
      dim_red <-  NULL
    } else {
      dim_red <- dim_red_for_clusters
    }

    spe <- build_metacluster(
      spe,
      param_name = paste0(param_name, "_metacluster"),
      cluster_name = cluster_name,
      cluster_type = "custom",
      dim_red = dim_red,
      assay_name = assay_name,
      clustering_metadata = metadata_method,
      max_k = max_k,
      distance = "euclidean",
      BPPARAM = BPPARAM,
      quiet = quiet,
      called = TRUE,
      ...
    )
  }
  return(spe)
}

# %% [markdown]
# ## build_metacluster
# 
# This function builds metacluster using ConsensusClusterPlus. <br>
# 
# **ARGUMENTS**
# - `spe` : A signleCellExperiment or SpatialExperiment object. Clustering must be recorded in the metadata (see multi_cluster)
# - `cluster_name` : character string indicating the name of the clustering to consider. Must be present in the metadata and in the colData
# - `dim_red` : character string, indicating the reducedDim to use (PCA/UMAP/HARMONY etc.) if NULL, will use an assay instead.
# - `assay_name` : character string, indicating which assay to use (if dim_red = NULL)
# - `max_k` : integer value, indicating the maximum number of metaclusters. Will test all k value from 2 to max_k
# - `clustering_metadata` : character string, indicating the name of metadata(spe) which contains the information about the clustering performed. Usual value within multi_cluster is "clustering_method"
# - `clusterAlg` : Character string, algorithm to use for clustering. <br>Possible values are "hc" for hierarchical clustering, "pam" for partitioning around medioids, "km" for k-means
# - `distance` : character string, how to compute distance. <br>Possible values are "pearson" (1- pearson correlation), "spearman" (1 - spearman correlation), "euclidean", "binary", "maximum", "minkowski", "camberra".
# - `seed` : positive integer value, seed for random generation
# - `quiet` : logical value : whether to hide or display messages during metaclustering
# - `...` : additional arguments passed to ConsensusClusterPlus
# 
# **VALUE**<br>
# Return the spe with the following :
# - **metaclustering set to TRUE** for `cluster_name` column in `metadata(spe)[[clustering_metadata]]`
# - **A new or updated "meta_code" element in `metadata(spe)`** name with `cluster_name`, that contains a dataframe of factors translating the clusters into metaclusters. For example the column meta20 of this dataframe will contain a factor translating the clusters (as rows) into metaclusters (as value of the factor) for k = 20. If needed, the cluster names are stored in the cluster column
# - **A new or updated "delta_area element in `metadata(spe)`**, containing the delta_area plot (created thanks to CATALYST) under the name "cluster_name"
# 
# **DETAILS**<br>
# This function has been highly motivated by the implementation of ConsensusClusterPlus by CATALYST::cluster. Compared to CATALYST's function, it is expanded to be compatible with any clustering made and recorded into a singleCellExperiment. It also stores nicely the metadata into lists that are accessible with the cluster_name. <br>
# For more details about ConsensusClusterPlus, please refer to their documentation. <br>
# For more details about the Consensus Clustering method, please see : <br>Monti, S., Tamayo, P., Mesirov, J. *et al.* Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data. *Machine Learning* 52, 91–118 (2003). https://doi.org/10.1023/A:1023949509487

# %%
build_metacluster <- function(
  spe,
  param_name = "metacluster",
  cluster_name = default_parameter(spe, param_name, "cluster_name"),
  cluster_type= default_parameter(spe, param_name, "cluster_type"),
  dim_red = NULL,
  assay_name = "counts",
  max_k = 3,
  clustering_metadata = "clustering_methods",
  clusterAlg = c("hc", "pam", "km"),
  distance = c("pearson", "spearman", "euclidean",
               "binary", "maximum", "minkowski", "camberra"),
  quiet = FALSE,
  BPPARAM = BiocParallel::SerialParam(),
  seed = 42,
  called = FALSE,
  ...
) {
  if (!called) {
    env <- environment()
    param_list <- get_called_param(
      env,
      excluded = c("spe", "param_name", "called")
    )
    spe <- spe |>
      address_parameters(
        param_name = param_name,
        parameters = param_list,
        envir = env
      )
  }
  if(!quiet) message("Metaclustering :", cluster_name, ".")
  #Attach packages
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  
  validate_build_metacluster(
    spe = spe,
    cluster_name = cluster_name,
    assay_name = assay_name,
    dim_red = dim_red,
    max_k = max_k,
    clustering_metadata = clustering_metadata,
    quiet = quiet,
    BPPARAM = BPPARAM
  )
  
  #clusterAlg and distance
  clusterAlg <- match.arg(clusterAlg)
  distance <- match.arg(distance)

  if (!quiet) message("Building matrix for :", cluster_name, ".")
  #NB : reducedDim natively in col = feature and row = observation
  #So need to transpose it

  if (is.null(dim_red)) {
    assay <- SummarizedExperiment::assay(spe, assay_name)
  } else {
    assay <- t(
      SingleCellExperiment::reducedDim(
        spe,
        dim_red
      )
    )
  }
  factor <- as.factor(
    SummarizedExperiment::colData(spe)[[cluster_name]]
  )
  

  if (!quiet) message("Aggregating cluster : ", cluster_name, ".")
  #summarizeAssayByGroup returns a summarizedExperiment object
  #With only one assay (if given a matrix), containing the mean expression matrix
  mat <- scuttle::summarizeAssayByGroup(
    assay,
    ids = factor,
    statistics = "mean",
    store.number = NULL,
    BPPARAM = BPPARAM
  ) %>%
    SummarizedExperiment::assay() #Extracts the only assay
  if (!quiet) message("Metaclustering : ", cluster_name, ".")
  # Warn when max_k is dangerous
  # The exact number that fail is hard to predict as it depends on random subsetting
  # inside consensusClusterPlus
  # Could maybe be technically handled via condition handlers
  # But the biological relevance of metaclustering a very small number of cluster is dubious
  # So warning the user may be better
  safe_k <- length(levels(factor)) - 10
  if(max_k > safe_k) {
    cli::cli_alert_warning(c(
      "Warning : Max_k close to the number of total clusters \n",
      "Metaclustering may fail and might not be relevant"
    ))
  }
  #Running ConsensusClusterPlus
  #The NULL pdf captures the unwanted plots (plots will be stored later)
  pdf(NULL)
  metaclusters <- suppressMessages(
    ConsensusClusterPlus::ConsensusClusterPlus(
      mat,
      maxK = max_k,
      reps = 100,
      distance = "euclidean",
      seed = seed,
      plot = NULL,
      ...
    )
  )
  graphics.off()


  if (!quiet) message("Storing metaclusters for :", cluster_name, ".")
  #prepare storage :
  #n_cluster tracks the number of initial cluster
  n_cluster <- spe %>%
    SummarizedExperiment::colData() %>%
    .[[cluster_name]] %>%
    unique() %>%
    length()
  #label correspond to all the k done (from 2 to max_k)
  label <- seq(from = 2, to = max_k)
  #Meta_code will contain the factor for translating a cluster to the corresponding metacluster
  meta_code <- data.frame(
    seq_len(n_cluster),
    purrr::map(metaclusters[-1], "consensusClass")
  ) %>%
    dplyr::mutate_all(\(x) factor(x, levels = sort(unique(x))))
  colnames(meta_code) <- c("cluster", sprintf("meta%s", label))

  delta_area <- metaclusters %>%
    CATALYST:::.plot_delta_area()

  #Store results in spe
  #Prepare meta_code metadata :
  spe <- spe %>%
    store_metadata(
      meta_code,
      category = "metaclustering",
      name = cluster_name,
      type = "meta_code"
    ) %>%
    store_metadata(
      delta_area,
      category = "metaclustering",
      name = cluster_name,
      type = "delta_area"
    ) %>%
    store_metadata( #Will set the metaclustering value to TRUE
      data = TRUE,
      category = clustering_metadata,
      name = cluster_name,
      type = "metaclustering"
    )

  if (!quiet) message("Metaclustering done :", cluster_name, ".")
  return(spe)
}

# %% [markdown]
# # Plotting
# 
# General consideration :
# - These functions are for the most part inspired by and/or wrapper around functions of the {dittoSeq} packages (dittoPlot, dittoDimPlot, dittoHeatmap)
# - With some rare exceptions these functions take a store_plot arguments can either write the plot(s) on a file by calling organize_and_store (if store_plot = TRUE) or return a plot/plot_list (if store_plot = FALSE)
# - Functions that do not take a store_plot arguments always returns the plot or plot_list

# %% [markdown]
# ## plot_expression
# 
# Wrapper around dittoPlot for the dittoSeq package, to perform violin, ridge, box or jitter plot. <br>
# 
# **ARGUMENTS**
# - `spe` : a SingleCellExperiment or SpatialExperiment
# - `assay_name` : character string, name of the assay to use for intensities/expression
# - `variable` : character string, name of the variable to plot. Currently only marker is supported, which will plot for each marker
# - `markers_to_plot` : character vector, indicating which markers to plot. If left to the default NULL, will plot all the markers in rownames(spe)
# - `plot_type` : character vector, indicating the type of plot to use, possible values are : ("ridgeplot", "vlnplot", "boxplot", "jitter"). <br>
# Order is important, dictates layer order. The first type will be in the back and the last in the front of the plot.
# - `group.by` : character string, corresponding to the colData to group the plot by.
# - `legend.show` : logical, whether to add the legend to the plot
# - `quiet` : logical, if FALSE will display progression messages
# - `store_plot` : logical, whether to store the plots. If FALSE, will return the list of plots
# - `plot_dir` : character string, path to the folder where the plots will be stored
# - `separate` : logical, when storing plot, whether to store each in a individual file or to store one file organized with grid.arrange
# - `device_type` : character string corresponding to the type of graphical device to use for storing. <br> Possible values are ("pdf", "postscript", "svg", "bitmap", "png", "jpeg", "bmp", "tiff")
# - `file_extension` : character string, extension of the file to use (like ".pdf" for pdf files)
# - `...` : additional arguments passed to dittoPlot and to device function
#   
# **VALUE**<br>
# If store_plot is TRUE, create new files and returns spe invisibly. <br>
# If store_plot is FALSE, returns a list of ggplot objects

# %%
plot_expression <- function(
  spe,
  param_name = "plots",
  assay_name = "counts",
  variable = c("marker"),
  markers_to_plot = rownames(spe),
  plot_type = c("vlnplot", "jitter"),
  group.by = "sample_id",
  legend.show = FALSE,
  quiet = FALSE,
  store_plot = FALSE,
  plot_dir = "./",
  separate = FALSE,
  device_type = "pdf",
  base_height = 5,
  base_width = 5,
  register_param = TRUE,
  register_env = parent.frame(),
  called = FALSE,
  ...
) {
  if (dev.cur() > 0) {
    on.exit(graphics.off(), add = TRUE)
  }
  if (!called) {
    env <- environment()
    param_list <- get_called_param(
      env,
      excluded = c("spe", "param_name", "register_param", "register_env")
    )
    spe <- spe |>
      address_parameters(
        param_name = param_name,
        parameters = param_list,
        envir = env
      )
  }
  

  if (!quiet) message("Preparing variables")
  #argument check :
  variable <- match.arg(variable)
  plot_type <- match_arg_multiple(
    plot_type,
    valid = "Plot types",
    choices = c("ridgeplot", "vlnplot", "boxplot", "jitter"),
    arg_name = "plot_type"
  )

  spe |>
    verify_assay(assay_name) |>
    verify_markers(markers_to_plot)

  if (!quiet) message("Preparing mapper")
  mapper <- prepare_mapper(
    spe,
    variable = variable,
    subset = markers_to_plot
  )


  #Prepare plot_lists
  if (!quiet) message("Preparing list of plots")
  plot_list <- purrr::map(
    mapper,
    function(x) {
      dittoSeq::dittoPlot(
        spe,
        var = x,
        assay = assay_name,
        group.by = group.by,
        plots = plot_type,
        legend.show = legend.show,
        ...
      )
    }
  )
  #Store plots

  if (store_plot) {
    #Store plots
    if (!quiet) message("Storing plots")

    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    expression_dir <- paste0(plot_dir, "expression/")
    if (!dir.exists(expression_dir)) dir.create(expression_dir)

    organize_and_store(
      plot_list,
      plot_name = mapper,
      plot_dir = expression_dir,
      separate = separate,
      device_type = device_type,
      base_height = base_height, #Always in inches
      base_width = base_width,
      ...
    )
    invisible(spe)
  } else {
    # If we return the plot list, we offer to reassign spe
    # To take advantage of the recorderded parameter
    if (register_param) assign("spe", spe, envir = register_env)
    return(plot_list)
  }
}

# %% [markdown]
# ### plot_multiple_expression
# 
# COnvinient wrapper around plot_expression to plot expression for multiple clusters, such as when benchmarking <br>
# Uniquely, does not use parameter itself, rather passes it to plot_expression, so the call will not be recorded
# 
# store_plot, here plural, will override the store_plot parameter.
# 
# Multiple subfolder will be created inside the plot_dir.

# %%
plot_multiple_expression <- function(
  spe,
  param_name = "plots",
  by = "colData",
  select = "patient_id",
  store_plot = FALSE,
  plot_dir = "./",
  quiet = FALSE,
  ...
) {
  mapper <- prepare_mapper(spe, variable = by, subset = select)

  if (store_plot) {
    # For store plot, we will create subfolder based on the mapper
    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    subfolders <- paste0(plot_dir, mapper, "/")

    spe <- purrr::reduce2(
      mapper,
      subfolders,
      function(spe, x, folder) {
        # We call repeatedly plot_expression, grouping by the mapper
        # And ploting in the subfolder
        if (!quiet) message(paste0("Plotting ", x))
        if (!dir.exists(folder)) dir.create(folder)
        spe <- plot_expression(
          spe,
          param_name,
          group.by = x,
          plot_dir = folder,
          store_plot = TRUE,
          ...
        )
        return(spe)
      },
      .init = spe
    )

    return(spe)
  } else {
    plot_list <- purrr::map(
      mapper,
      function(x) {
        plot_list <- plot_expression(
          spe,
          param_name,
          group.by = x,
          store_plot = FALSE,
          ...
        )
      }
    ) %>%
      setNames(mapper)
    return(plot_list)
  }
}

# %% [markdown]
# ## Plot_marker_flowsom

# %%
plot_marker_flowsom <- function(
  spe,
  param_name = "plots",
  cluster_name = default_parameter(spe, param_name, "cluster_name"),
  markers_to_plot = NULL,
  store_plot = FALSE,
  plot_dir = "./",
  separate = TRUE,
  device_type = "pdf",
  base_height = 5,
  base_width = 5,
  quiet = TRUE,
  register_param = TRUE,
  register_env = parent.frame(),
  ...
) {
  #Wrapper aroung FlowSOM::PlotMarker to fetch the correct fsom object
  #Returns a list of ggplot, one for each marker in markers_to_plot
  # ... is passed to PlotMarker for additional customization.
  #Checks fsom_name argument
  env <- environment()
  param_list <- get_called_param(
    env,
    excluded = c("spe", "param_name", "register_param", "register_env")
  )
  spe <- spe |>
    address_parameters(
      param_name = param_name,
      parameters = param_list,
      envir = env
    )

  fsom <- spe |>
    verify_metadata(
      cluster_name,
      category = "clustering_data",
      type = "fsom"
    ) |>
    get_metadata(
      category = "clustering_data",
      name = cluster_name,
      type = "fsom"
    )

  #Checks if fsom is a FlowSOM object
  if (!inherits(fsom, "FlowSOM")) {
    msg <- c(
      "Incorrect fsom object",
      "x" = "fsom class is : {class(fsom)}",
      "i" = "fsom must be a fsom object"
    )
    cli::cli_abort("error_fsom", message = msg)
  }
  if (!quiet) message("Preparing mapper")
  mapper <- prepare_mapper(spe, variable = "marker", subset = markers_to_plot)
  #Plot each marker in markers_to_plot

  if(!quiet) message("Building plots")
  plot_list <- purrr::map(
    mapper,
    function(x) {
      FlowSOM::PlotMarker(fsom, x, ...)
    }
  )

  if (store_plot) {
    #Store plots
    if (!quiet) message("Storing plots")

    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    flowsom_dir <- paste0(plot_dir, "flowsom/")
    if (!dir.exists(flowsom_dir)) dir.create(flowsom_dir)
    
    organize_and_store(
      plot_list,
      plot_name = mapper,
      plot_dir = flowsom_dir,
      separate = separate,
      device_type = device_type,
      base_height = base_height, #Always in inches
      base_width = base_width,
      ...
    )
    invisible(spe)
  } else {
    # If we return the plot list, we offer to reassign spe
    # To take advantage of the recorderded parameter
    if (register_param) assign("spe", spe, envir = register_env)
    names(plot_list) <- mapper
    return(plot_list)
  }
}

# %% [markdown]
# ## plot_mst

# %%
plot_mst <- function(
  spe,
  param_name = "plots",
  cluster_name,
  assay_name = "counts",
  color_by = c("none", "marker", "metacluster"),
  markers = NULL,
  col_pal = viridis(100),
  background_by = c("none", "metacluster"),
  metacluster_k = 20,
  background_pal = rainbow(metacluster_k),
  store_plot = FALSE,
  plot_dir = "./",
  separate = TRUE,
  device_type = "pdf",
  base_height = 5,
  base_width = 5,
  quiet = TRUE,
  register_param = TRUE,
  register_env = parent.frame()
) {
  ### PREPARATION ####################################
  #Check arguments :
  env <- environment()
  param_list <- get_called_param(
    env,
    excluded = c("spe", "param_name", "register_param", "register_env")
  )
  spe <- spe |>
    address_parameters(
      param_name = param_name,
      parameters = param_list,
      envir = env
    )
  fsom <- spe |>
    verify_metadata(
      cluster_name,
      category = "clustering_data",
      type = "fsom"
    ) |>
    get_metadata(
      category = "clustering_data",
      name = cluster_name,
      type = "fsom"
    )
  #Checks if fsom is a FlowSOM object
  if (!inherits(fsom, "FlowSOM")) {
    msg <- c(
      "Incorrect fsom object",
      "x" = "fsom class is : {class(fsom)}",
      "i" = "fsom must be a fsom object"
    )
    cli::cli_abort("error_fsom", message = msg)
  }

  #Color_by and background_by
  color_by <- match.arg(color_by)
  background_by <- match.arg(background_by)

  if (color_by == "marker") {
    spe |>
      verify_assay(assay_name)

    #Prepare the mean_assay :
    mean_assay <- scuttle::summarizeAssayByGroup(
      spe,
      ids = SummarizedExperiment::colData(spe)[[cluster_name]],
      statistics = "mean",
      assay.type = assay_name,
      subset.row = markers
    ) %>%
      SummarizedExperiment::assay("mean")
  }

  #prepare the metacluster factor :
  if (background_by == "metacluster" || color_by == "metacluster") {
    metacluster <- spe |>
      verify_metadata(
        cluster_name,
        category = "clustering_data",
        type = "meta_code"
      ) |>
      get_metadata(
        category = "clustering_data",
        name = cluster_name,
        type = "meta_code"
      ) |>
      dplyr::pull(
        dplyr::all_of(paste0("metacluster_", metacluster_k))
      )
  } else {
    metacluster <- NULL
  }

  if (color_by == "marker") {
    plot_list <- purrr::map(
      markers,
      function(x) {
        plot_mst_layers(
          fsom = fsom,
          color_by = color_by,
          mean_assay = mean_assay,
          marker = x,
          col_pal = col_pal,
          background_by = background_by,
          metacluster = metacluster,
          background_pal = background_pal
        ) |>
          setNames(x)
      }
    )
  } else {
    plot_list <- plot_mst_layers(
      fsom = fsom,
      color_by = color_by,
      mean_assay = mean_assay,
      marker = markers,
      col_pal = col_pal,
      background_by = background_by,
      metacluster = metacluster,
      background_pal = background_pal
    ) %>%
      list() %>%
      setNames(paste0(cluster_name, "_", color_by))
    #Returning as a list (of 1) to be compatible with gridExtra
  }

  if (store_plot) {
    #Store plots
    if (!quiet) message("Storing plots")

    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    mst_dir <- paste0(plot_dir, "mst/")
    if (!dir.exists(mst_dir)) dir.create(mst_dir)

    organize_and_store(
      plot_list,
      plot_name = names(plot_list),
      plot_dir = mst_dir,
      separate = separate,
      device_type = device_type,
      base_height = base_height, #Always in inches
      base_width = base_width,
      ...
    )
    invisible(spe)
  } else {
    # If we return the plot list, we offer to reassign spe
    # To take advantage of the recorderded parameter
    if (register_param) assign("spe", spe, envir = register_env)
    return(plot_list)
  }
}

# %% [markdown]
# ### plot_mst_layers

# %%
plot_mst_layers <- function(
  fsom,
  color_by = c("none", "marker", "metacluster"),
  mean_assay = NULL,
  marker = NULL,
  col_pal = viridis::viridis(100),
  background_by = c("none", "metacluster"),
  metacluster = NULL,
  background_pal = if(!is.null(metacluster)) rainbow(nlevels(metacluster)) else NULL #nolint
){ 
  #Check arguments :
  #Check color arguments :
  color_by <- match.arg(color_by)

  if (color_by == "marker") {
    mean_assay_must <- list(
      "be_found" = \(x) !is.null(x),
      "be_numerical" = \(x) is.numeric(x),
      "have_one_value_per_cluster" = \(x) ncol(x) == max(fsom$map$mapping[, 1])
    )
    mean_assay <- must_multiple("mean_assay", mean_assay_must, mean_assay)
    if (!marker %in% rownames(mean_assay)) {
      abort_not_found(arg_name = "marker",
                      valid = "marker",
                      choices = rownames(mean_assay),
                      not = marker)
    }
  }

  #Check background arguments :
  background_by <- match.arg(background_by)
  #Check metacluster
  if (background_by == "metacluster" ||
        color_by == "metacluster") {
    metacluster_must = list(
      "be_found" = \(x) !is.null(x),
      "be_factor" = \(x) is.factor(x),
      "have_one_value_per_cluster" = \(x) length(x) == max(fsom$map$mapping[, 1]) #nolint
    )
    metacluster <- must_multiple("metacluster", metacluster_must, metacluster)
  }
  
  #Prepare col_value and back_value
  col_values <- switch(
    color_by,
    "marker" = mean_assay[marker, ],
    "metacluster" = metacluster,
    "none" = ,
    .default = NULL
  )
  back_values <- switch(
    background_by,
    "metacluster" = metacluster,
    "none" = ,
    .default = NULL
  )
  #Create the plot layer by layer
  label <- switch(
    color_by,
    "marker" = marker,
    "metacluster" = "metaclusters",
    .default = NULL
  )

  FlowSOM::PlotFlowSOM(fsom, view = "MST") |>
    FlowSOM::AddBackground(
      backgroundValues = back_values,
      backgroundColors = background_pal
    ) |>
    FlowSOM::AddNodes(
      values = col_values,
      colorPalette = col_pal,
      label = label
    )
  #Will return the ggplot object.
}

# %% [markdown]
# ## plot_umap

# %%
plot_umap <- function(
  spe,
  param_name = "plots",
  by = c("colData", "reducedDim", "marker"),
  assay_name = "counts",
  reducedDim = "UMAP",
  reducedDim_col = "patient_id",
  address_umap = FALSE,
  select = NULL,
  store_plot = FALSE,
  plot_dir = "./",
  separate = TRUE,
  device_type = "pdf",
  base_height = 5,
  base_width = 5,
  quiet = TRUE,
  register_param = TRUE,
  register_env = parent.frame(),
  ...
) {
  env <- environment()
  param_list <- get_called_param(
    env,
    excluded = c("spe", "param_name", "register_param", "register_env")
  )
  spe <- spe |>
    address_parameters(
      param_name = param_name,
      parameters = param_list,
      envir = env
    )

  by <- match.arg(by)
  if (by == "marker") verify_assay(spe, assay_name)
  if (!is.null(select)) {
    switch(
      by,
      "colData" = spe |> verify_multiple_colData(select),
      "reducedDim" = purrr::map(select, \(x) verify_reducedDim(spe, x)),
      "marker" = spe |> verify_markers(select)
    )
  }
  if (by != "reducedDim") verify_reducedDim(spe, reducedDim)

  # prepare mapper
  mapper <- prepare_mapper(spe, variable = by, subset = select)
  # For colData also prepare which UMAP to plot depending on address_to_umap
  if(by == "colData") {
    if(address_umap) {
      umaps <- address_to_umap(mapper, spe)
    } else {
      umaps <- rep_len(reducedDim, length(mapper))
    }
  }

  if (by == "reducedDim") {
    plot_list <- purrr::map(
      mapper,
      function(x) {
        dittoSeq::dittoDimPlot(
          spe,
          var = reducedDim_col,
          reduction.use = x,
          ...
        ) +
          ggtitle(x)
      }
    ) |> setNames(mapper)
  } else if (by == "marker") {
    plot_list <- purrr::map(
      mapper,
      function(x) {
        dittoSeq::dittoDimPlot(
          spe,
          var = x,
          assay = assay_name,
          reduction.use = reducedDim,
          ...
        ) +
          ggtitle(x)
      }
    ) |> setNames(mapper)
  } else if (by == "colData") {
    plot_list <- purrr::map2(
      mapper,
      umaps,
      function(x, y) {
        dittoSeq::dittoDimPlot(
          spe,
          var = x,
          assay = assay_name,
          reduction.use = y,
          ...
        ) +
          ggtitle(x)
      }
    ) |> setNames(mapper)
  }

  if (store_plot) {
    #Store plots
    if (!quiet) message("Storing plots")
    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    umap_dir <- paste0(plot_dir, "umap/")
    if (!dir.exists(umap_dir)) dir.create(umap_dir)

    organize_and_store(
      plot_list,
      plot_name = names(plot_list),
      plot_dir = umap_dir,
      separate = separate,
      device_type = device_type,
      base_height = base_height, #Always in inches
      base_width = base_width,
      ...
    )
    invisible(spe)
  } else {
    # If we return the plot list, we offer to reassign spe
    # To take advantage of the recorderded parameter
    if (register_param) assign("spe", spe, envir = register_env)
    return(plot_list)
  }
}

# %% [markdown]
# ## plot_expr_heatmap

# %%
plot_expr_heatmap <- function(
  spe,
  param_name = "plots",
  select = default_parameter(spe, param_name, "select"),
  assay_name = "counts",
  exclude_row = NULL,
  color_pal = viridis::viridis(100),
  store_plot = FALSE,
  plot_dir = "./",
  separate = TRUE,
  device_type = "pdf",
  base_height = 5,
  base_width = 5,
  quiet = TRUE,
  register_param = TRUE,
  register_env = parent.frame(),
  ...
) {
  #stop device on exit
  if (dev.cur() > 0) {
    on.exit(graphics.off(), add = TRUE)
  }
  env <- environment()
  param_list <- get_called_param(
    env,
    excluded = c("spe", "param_name", "register_param", "register_env")
  )
  spe <- spe |>
    address_parameters(
      param_name = param_name,
      parameters = param_list,
      envir = env
    )

  spe |>
    verify_assay(assay_name) |>
    verify_multiple_colData(select)
  #Check arguments :
  if (!is.null(exclude_row)) verify_markers(spe, exclude_row)

  mapper <- prepare_mapper(spe, variable = "colData", subset = select)
  row_subset <- rownames(spe)[!rownames(spe) %in% exclude_row]

  #Aggreagates the mean expression for each cluster
  if (!quiet) message("Computing expression values")
  mean_list <- purrr::map(
    mapper,
    function(x) {
      scater::aggregateAcrossCells(
        as(spe, "SingleCellExperiment"),
        ids = SummarizedExperiment::colData(spe)[[x]],
        use.assay.type = assay_name,
        subset.row = row_subset,
        statistics = "mean"
      )
    }
  )
  #Prepare plot list
  if (!quiet) message("Preparing plots")
  plot_list <- purrr::map2(
    mean_list,
    mapper,
    function(x, y) {
      ggplotify::as.ggplot(
        dittoSeq::dittoHeatmap(
          x,
          assay = assay_name,
          cluster_cols = TRUE,
          scale = "row",
          heatmap.colors.max.scaled = color_pal,
          annot.by = c(y, "ncells"),
          ...
        )
      ) + ggplot2::ggtitle(y)
    }
  )
  names(plot_list) <- mapper

  if (store_plot) {
    #Store plots
    if (!quiet) message("Storing plots")
    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    ht_dir <- paste0(plot_dir, "expression_heatmap/")
    if (!dir.exists(ht_dir)) dir.create(ht_dir)

    organize_and_store(
      plot_list,
      plot_name = names(plot_list),
      plot_dir = ht_dir,
      separate = separate,
      device_type = device_type,
      base_height = base_height, #Always in inches
      base_width = base_width,
      ...
    )
    invisible(spe)
  } else {
    # If we return the plot list, we offer to reassign spe
    # To take advantage of the recorderded parameter
    if (register_param) assign("spe", spe, envir = register_env)
    return(plot_list)
  }
}

# %% [markdown]
# ## Plot_cluster_tracker
# 
# Function that plots the cluster_tracker plot, equivalent to what ConsensusClusterPlus does, but using pheatmap instead.
# 
# **ARGUMENTS** <br>
# - `spe` : a SingleCellExperiment or SpatialExperiment object
# - `cluster_name` : Character string, indicating which cluster to plot. The cluster must have metaclustering done, with metadata matching the output of `build_metacluster`. <br> See that function for more information.
# - `color` : Vector of colors to use.
# - `cluster_cols` and  `cluster_rows` : Logicles, whether to cluster columns or rows of the heatmap. Column will represent the clusters values in cluster_name, while rows will be the metaclusters. <br> Leaving those parameters to defaults allows to reproduce a similar plot as ConsensusClusterPlus output.
# - `...` : additional arguments passed to pheatmap.
#   
# **VALUE** <br>
# A pheatmap object. By default pheatmap also prints plot to screen on top of returning the object.
# 
# **DETAILS**<br>
# The plot is made to visualize the stability of the metaclustering. For each cluster in `cluster_name` (columns) you will be able to visualize the attributed metacluster for each k values (rows), from  2 to maxK. The value of maxK (max_k in build_metacluster) is retrieved directly from the metadata. <br>
# This function is motivated to recreate a plot done by ConsensusClusterPlus. However, ConsensusClusterPlus creates it's plot directly while running the clustering algorithm, which might not be the desired behavior, if you plan to cluster first and visualize later down the pipeline. <br>
# ConsensusClusterPlus also takes a convoluted way around making this plot, which makes it difficult to reproduce, customize or adapt it. Here we make use a heatmap (powered by pheatmap) to achieve a similar results with more flexibility, albeit with a bit less visual consistency. <br>
# For more details about ConsensusClusterPlus, please refer to their documentation.

# %%
plot_cluster_tracker <- function(spe,
                                 cluster_name,
                                 color = rainbow(max_k),
                                 cluster_cols = TRUE,
                                 cluster_rows = FALSE,
                                 metadata_name = "clustering_methods",
                                 ...) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))

  stopifnot(inherits(spe, "SingleCellExperiment"))
  metacluster_done <- get_metadata(
    spe,
    category = metadata_name,
    name = cluster_name,
    type = "metaclustering"
  )

  if (!metacluster_done) {
    msg <- c(
      "Error : no metaclustering recorded for {cluster_name}",
      "x" = "Please make sure metaclustering has been performed",
      "i" = "Appropriate metadata also need to be present",
      "i" = "See build_metacluster for more information"
    )
    cli::cli_abort(msg)
  }
  #Extract codes from metadata
  code <- get_metadata(
    spe,
    category = "metaclustering",
    name = cluster_name,
    type = "meta_code"
  )
  #prepare matrix for plot
  mat <- code %>%
    purrr::map(S4Vectors::unfactor) %>%
    purrr::map(as.integer) %>%
    as.data.frame() %>%
    t()
  #prepare colors
  #First row being 2 metaclusters
  max_k <- nrow(mat) + 1
  color <- color[seq_len(max_k)]

  pheatmap::pheatmap(
    mat[-1, ],
    color = color,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    ...
  )
}

# %% [markdown]
# ## Plot_delta_area

# %%
plot_delta_area <- function(
  spe,
  cluster_names,
  store_plot = FALSE,
  plot_dir = "./",
  separate = TRUE,
  device_type = "pdf",
  base_height = 5,
  base_width = 5,
  quiet = TRUE
) {
  #Prepare the names of delta area plots
  #Either select all delta area plots
  plot_list <- purrr::map(
    cluster_names,
    function(x) {
      spe |>
        verify_metadata(
          x,
          category = "clustering_data",
          type = "delta_area"
        ) |>
        get_metadata(
          category = "clustering_data",
          name = x,
          type = "delta_area"
        )
    }
  ) |>
    setNames(cluster_names)
  
  if (store_plot) {
    #Store plots
    if (!quiet) message("Storing plots")

    if (!dir.exists(plot_dir)) dir.create(plot_dir)
    delta_dir <- paste0(plot_dir, "delta_area/")
    if (!dir.exists(delta_dir)) dir.create(delta_dir)

    organize_and_store(
      plot_list,
      plot_name = names(plot_list),
      plot_dir = delta_dir,
      separate = separate,
      device_type = device_type,
      base_height = base_height, #Always in inches
      base_width = base_width,
      ...
    )
    invisible(spe)
  } else {
    return(plot_list)
  }
}

# %% [markdown]
# ## plot_snr_cells
# 
# 

# %%
plot_snr_cells <- function(
  spe,
  trans_assay = "transformed",
  count_assay = "counts",
  seed = 42
) {
  if (dev.cur() > 0) {
    on.exit(graphics.off(), add = TRUE)
  }
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  
  stopifnot(inherits(spe, "SingleCellExperiment"))
  # Create a tibble to use in the ggplot
  cur_snr <- compute_snr(spe, trans_assay, count_assay, seed)
  
  p <- cur_snr %>% ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(log2(signal), log2(snr))) +
    ggrepel::geom_label_repel(
      ggplot2::aes(
        log2(signal),
        log2(snr),
        label = marker
      )
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::ylab("Signal-to-noise ration [log2]") +
    ggplot2::xlab("Signal intensity [log2]")

  return(p)
}

# %% [markdown]
# ### compute_snr
# 
# Function used to compute a signal/noise ratio based on a predictive model powered by Mclust. Returns a tibble with the results for each marker.
# 
# **ARGUMENTS**
# - `spe` : a SpatialExperiment or SingleCellExperiment object
# - `transfromed_assay` and `count_assay` : Character strings, respectively the name of the transformed assay (log, asihn, etc.) and the name of the counts assay (raw counts)
# - `seed` : numerical value used to seed random number generation
# 
# **VALUE** <br>
# Returns a tibble with 3 columns :
# - "snr" : the signal/noise ratio
# - "signal" : the mean intensity of the signal
# - "marker" : the name of the marker

# %%
compute_snr <- function(
  spe,
  trans_assay = "transformed",
  count_assay = "counts",
  seed = 42
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  suppressPackageStartupMessages(stopifnot(require("mclust")))

  spe %>%
    verify_assay(trans_assay) %>%
    verify_assay(count_assay)

  set.seed(seed)

  signals <- purrr::map(#for each marker
    rownames(spe),
    function(x) {
      #extract the transformed and non transformed values
      cur_exprs <- SummarizedExperiment::assay(spe, trans_assay)[x, ]
      cur_counts <- SummarizedExperiment::assay(spe, count_assay)[x, ]
      #Train a model on transformed:
      cur_model <- mclust::Mclust(cur_exprs, G = 2)
      #calculate the mean of counts for each class :
      mean1 <- mean(cur_counts[cur_model$classification == 1])
      mean2 <- mean(cur_counts[cur_model$classification == 2])
      #Attributes the highest mean to signal, the lowest to noise
      signal <- ifelse(mean1 > mean2, mean1, mean2)
      noise <- ifelse(mean1 > mean2, mean2, mean1)
      #return the ratio and the signal
      return(
        list(
          snr = signal / noise,
          signal = signal,
          x
        )
      )
    }
  )
  # Returns a tibble to use in the ggplot
  tibble::tibble(
    "snr" = purrr::map_dbl(signals, \(x) x[[1]]),
    "signal" = purrr::map_dbl(signals, \(x) x[[2]]),
    "marker" = purrr::map_chr(signals, \(x) x[[3]])
  )
}

# %% [markdown]
# # Benchmarking

# %% [markdown]
# ## Benchmark_clustering
# 
# Wrapper function aroud preprocessing, clustering and the primary plots. 
# 
# Supported plots : heatmap, violin, and umap
# 
# This function utilises param_name (default "benchmarking") for its function call, and provided preprocessing_param to preprocessing. <br>
# It will create a set of cluster parameters based on benched_param. Common_param can be used to apply some parameters across all clustering (on top of global_parameters). <br>
# Finally, it will call plot function using plot_param. Plot_param can be a single value, which will be used for all the plots in get_plots, or one name for each plot within get_plots
# 
# If store_plot is true, will plot in subfolders of plot_dir <br>
# Else the plots will be store in the metadata of spe, category = benchmark_plots

# %%
benchmark_clustering <- function(
  spe,
  param_name = "benchmarking",
  preprocessing_param = "preprocessing",
  skip_preprocess = FALSE,
  plot_param = "plots",
  benched_names = default_parameter(spe, param_name, "benched_names"),
  benched_param = default_parameter(spe, param_name, "benched_param"),
  common_param = default_parameter(spe, param_name, "common_param"),
  get_plots = "heatmap", # possible values : heatmap, violin, umap, mst (only for flowsom)
  quiet = FALSE,
  use_channel = rownames(spe),
  store_plot = FALSE,
  plot_dir = "./",
  separate = FALSE,
  device_type = "pdf",
  ... #additional arguments passed to plot and device functions
) {
  # adress parameters
  env <- environment()
  param_list <- get_called_param(
    env
  )
  spe <- spe |>
    address_parameters(
      param_name = param_name,
      parameters = param_list,
      envir = env
    )

  # check arguments :
  get_plots <- tolower(get_plots)
  get_plots <- match_arg_multiple(
    get_plots,
    valid = "Plot types",
    choices = c("heatmap", "umap", "violin"),
    arg_name = "get_plots"
  )

  if (length(get_plots) > 1 &&
        length(plot_param) > 1 &&
        length(get_plots) != length(plot_param)) {
    cli::cli_abort(c(
      "Incorrect plot_param argument",
      "x" = "plot_param must be length 1 or the same length as get_plots",
      "i" = "Length of get_plots = {length(get_plots)}",
      "i" = "Length of plot_param = {length(plot_param)}"
    ))
  }
  benched_names <- tolower(benched_names)
  # Prepare a vector of plot name to iterate along get_plots, if necessary
  if (length(plot_param) == 1) {
    plot_param <- rep_len(plot_param, length(get_plots))
  }
  # Transforms the logical into character if required
  if (is.logical(use_channel)) use_channel <- rownames(spe)[use_channel]

  # Prepare parameters:
  if (!quiet) message("Preparing benchmark parameters")

  spe <- set_benchmark_parameters(
    spe,
    param_names = benched_names,
    benched_param = benched_param,
    common_param = common_param
  )

  # preprocessing
  if (!skip_preprocess) {
    if(!quiet) message("Preprocessing data")
    spe <- preprocess_data(spe, preprocessing_param)
  }


  #Clustering
  if(!quiet) message("Clustering")
  spe <- purrr::reduce(
    benched_names,
    function(spe, x) perform_clustering(spe, x),
    .init = spe
  )

  # Plotting
  if (!quiet) message("Plotting")

  # If the plot are stored, call the functions directly
  if (store_plot) {
    spe <- purrr::reduce2(
      get_plots,
      plot_param,
      function(spe, plot, param) {
        # Calls corresponding function depending on plot_type
        if(plot == "heatmap") {
          spe <- plot_expr_heatmap(
            spe,
            param,
            select = benched_names,
            ...
          )
        } else if(plot == "umap") {
          spe <- plot_umap(
            spe,
            param,
            by = "colData",
            select = benched_names,
            ...
          )
        } else if (plot == "violin") {
          spe <- plot_multiple_expression(
            spe,
            param,
            markers_to_plot = use_channel,
            by = "colData",
            select = benched_names,
            store_plot = TRUE,
            plot_dir = plot_dir,
            quiet = quiet,
            ...
          )
        } 
      },
      .init = spe
    )
  } else {
    # if we don't store the plot, we will prepare a list of plot
    # And then address it to the metadata of spe


    plot_list <- purrr::map2(
      get_plots,
      plot_param,
      function(plot, param) {
        # Calls corresponding function depending on plot_type
        if(plot == "heatmap") {
          plot_expr_heatmap(
            spe,
            param,
            select = benched_names,
            ...
          )
        } else if(plot == "umap") {
          plot_umap(
            spe,
            param,
            by = "colData",
            select = benched_names,
            ...
          )
        } else if (plot == "violin") {
          plot_multiple_expression(
            spe,
            param,
            markers_to_plot = use_channel,
            by = "colData",
            select = benched_names,
            store_plot = FALSE,
            quiet = quiet,
            ...
          )
        } 
      }
    ) %>% #Fix names of the list
      setNames(get_plots)
    
    # Transpose the list to be organized by clustering rather than by plot
    plot_list <- purrr::transpose(plot_list)
    # We have a list named after each benched_names
    # Each element contain every plot required in get_plot
    
    spe <- purrr::reduce2(
      plot_list,
      names(plot_list),
      function(spe, plot, cluster) {
        spe <- spe %>%
          store_metadata(
            plot,
            category = "benchmark_plots",
            name = cluster
          )
      },
      .init = spe
    )
  }
  
  return(spe)
}

# %% [markdown]
# ### set_benchmark_parameters
# 
# Helper to create several parameters sets to benchmark. <br>
# You need to provide a list of names and a list of parameter sets. <br>
# benched_param is a list organised as arg1 = c(A, B, C), arg2 = c(1, 2, 3) <br>
# will create sets of parameters where (arg1 = A, arg2 = 1) etc. <br>
# A summary of the benched parameters will be available in parameters(spe, "benchmarking")
# 
# Common_param are parameters shared by all the benchmarked.

# %%
set_benchmark_parameters <- function(
  spe,
  param_names,
  benched_param,
  common_param
) {
  param_list <- transpose(benched_param) |>
    setNames(param_names)
  stopifnot(length(param_names) == length(param_list))
  
  spe <- purrr::reduce2(
    param_names,
    param_list,
    function(spe, name, param) {
      spe <- spe |>
        set_parameters(name, param) |>
        set_parameters(name, common_param)
      
      return(spe)
    },
    .init = spe
  )
  
  parameters(spe, "benchmarking") <- benched_param
  
  return(spe)
}

# %% [markdown]
# # Stats

# %% [markdown]
# ## compare_indication
# 
# Wrapper to perform a statistical comparison of abundance (count, proportion or density) of categories stored in colData between two indications.
# The results are printed on the console and stored in the metadata of the object
# 
# **Arguments**
# - `spe` : a SpatialExperiment object (compatible with SingleCellExperiment and SummarizedExperiment)
# - `sol_to_compare` : character string, name of the colData for which the comparison is to be made
# - `indication_col` : character string, name of the colData which contains the information about the indication (patient groups)
# - `patient_id_col` : character string, name of the colData which contains the identifier of the patients
# - `paired` : logical, whether to perform a paired or unpaired test. 
# - `sample_id_col` : character string, name of the colData which contains the identifier of the samples for each patient. Is ignored for unpaired test, as only patient_id is used for unpaired tests
# - `clever_test` : logical, whether to use an adaptative test based on the results of the normality assessment
# - `force_test` : character string, test to perform if clever_test is FALSE. Supports wilcoxon, t test and Welch (t test with unequal variance)
# - `count_type` : character string, what metric to consider : counts, proportion or density. See below
# - `tissue_area_name` : character string, name of metadata(spe) that contains information on the tissue area. See below for format. 
# - `quiet` : logical, whether to hide additional information
# - `shapiro_alpha` : numerical, value to use as an alpha risk for shapiro test. See normality assessment.
# - `skew_threshold` : numerical, threshold of absolute skew to reject normal distribution. See normality assessment.
# - `kurtosis_range` : numerical vector of length 2, boundaries of absolute kurtosis beyond which to reject normal distribution. See normality assessment.
# - `adjust` : character string, which function to use or not ("none") for p_value adjustement for multiple testing. See p.adjust() documentation for more detail
# - `signif_threshold` : numerical, what threshold to consider the adjusted pvalue significant 
# 
# **Value** <br>
# Return the spe object with the following :
# - in metadata(spe)$normality_check[[col_to_compare]] a list of 5 metrics for normality assessment (see below)
# - in metadata(spe)$statistical_test[[col_to_compare]] a tibble containing the results of the test with for each feature, the mean value for either group, the p value, the adjusted p value, whether it is significant and the statistical test used.
# 
# This function also prints a data.frame with the test results unless quiet = TRUE
# 
# **Details** <br>
# 
# *Normality assessment*
# Parametric test such as t test can only be used for normally distributed data, otherwise non parametric test such as wilcoxon test should be used. <br>
# This function always attempts to estimate whether the data follows a normal distribution, even if clever_test is FALSE, so that you can verify your assumption after the test is performed.
# The information is recorded in metadata(spe)$normality_check[[col_to_compare]] with the following :
# - $normality_data : a dataframe containing the values of shapiro test pvalue, skewedness and kurtosis of the data
# - $normality_check : a dataframe containing the assessment of each metric. FALSE denotes that normal distribution has been excluded
# - $both_normal : a dataframe containing for each feature whether both groups have been estimated as a normal distribution (TRUE) or not (FALSE)
# - $qq_plot : a ggplot object containing a qqplot, which is a representation of data distribution compared to a normal distribution
# - $equal_var : a dataframe containing for each feature, whether equal variance has been rejected (FALSE) or not (TRUE) using bartlett test
# 
# Of note, there is no absolute metric that can individually garantee to recognize a normal distribution, and the evaluation should be made with a combination of :
# - The shapiro test checks if the distribution is different from the normal distribution (null hypothesis : it is not different). That test can reject the null hypothesis for very small and non relevant deviation to the normal distribution for very large sample size, and at the opposite can fail to reject the null hypothesis for very small sample size
# - Skewedness is the evaluation of the assymetry of the data. Normal distribution has a skewedness = 0. If your data is skewed, it may not follow a normal distribution.
# - Kurtosis evaluate the curve of the distribution, and the amount of tailing values. Normal distribution has a kurtosis of 3. If your data has a very different kurtosis, it might not follow a normal distribution
# - Visual inspection of the distribution compared to a normal distribution with the qqplots are usefull to verify the normality assumption, in complement to previous metrics, none of which are perfect
# 
# Student t test also assumes equal variance between the two groups. If equal variance is rejected, a Welch test (generalized t test for unequal variance) should be performed
# This function assess equal variance using bartlett test, which similarly to shapiro's test can reject equal variance, but the sensibility and specificity will depend on the sample size.
# 
# *count_type* <br>
# This function can consider 3 metrics for comparison of the categories :
# - The absolute counts of cell of one category
# - The proportion of cells of one category, relative to the total cells
# - The density of cells of one category, relative to the tissue area (/mm², /µm², etc.)
# The results are not stored inside spe (only the mean value for each indication), but can be computed using the helper function prepare_counts()
# 
# *tissue_area_name* <br>
# This information must be a named vector with a value for each patient_id (unpaired) or sample_id(paired)

# %%
compare_indication <- function(
  spe,
  col_to_compare,
  indication_col = "indication",
  patient_id_col = "patient_id",
  paired = FALSE,
  sample_id_col = if (!paired) patient_id_col else "sample_id",
  skip_normality = FALSE,
  clever_test = TRUE,
  force_test = c("wilcox", "t_test", "welch"),
  count_type = c("count", "proportion", "density"),
  tissue_area_name = "tissue_area",
  quiet = FALSE,
  shapiro_alpha = 0.05,
  skew_threshold = 1,
  kurtosis_range = c(2, 4), #Absolute not excess
  adjust = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), #nolint
  signif_threshold = 0.05
) {
  #Could generalize later but for now we assume : 
  ##If paired, you have mutliple sample_id (images) by patient_id
  ## If unpaired, you have one sample by patient
  ##In reality you can have mutliple ROIs for each patient, 
  ##but if unpaired analysis, it is the same as considering them different patient
  ## Hence, you can bypass this difficulty by setting something like "patient_id_col = 'sample_id'"


  validate_compare_indication(
    spe, col_to_compare, indication_col, patient_id_col,
    paired, sample_id_col, clever_test, force_test,
    count_type, tissue_area_name, quiet, adjust,
    signif_threshold, skip_normality
  )

  if (paired) {
    #if paired, we need to have both patient_id and sample_id
    cluster_df <- SummarizedExperiment::colData(spe) %>%
      tibble::as_tibble() %>%
      dplyr::select(
        {{patient_id_col}},
        {{sample_id_col}},
        {{indication_col}},
        {{col_to_compare}},
      )
  } else {
    #If unpaired, patient_id is supposed to be the study level
    #We create a fictive sample_id col that will be ignored for the most part
    #Because group_by will not create new group with it
    #This makes the rest of the code consistent for both option
    if (!quiet) {
      message("Unpaired test : study assumed to be at the patient_id level. \n
              'sample_id_col' argument will be ignored.")
    }
    cluster_df <- colData(spe) %>%
      tibble::as_tibble() %>%
      dplyr::select(
        {{patient_id_col}},
        {{indication_col}},
        {{col_to_compare}},
      ) %>% #add a fictive column
      dplyr::mutate(sample_id = !!rlang::sym(patient_id_col))
    #Reallocate sample_id
    #Not done in the default value to avoid
    # throwing an error if sample_id col does not exists
    sample_id_col <- "sample_id"
  }
  
  cluster_df <- cluster_df %>%
    dplyr::group_by(
      !!rlang::sym(patient_id_col), #For unpaired sample_id = patient_id
      !!rlang::sym(sample_id_col),  #So only grouped once by patient_id
      !!rlang::sym(col_to_compare),
      !!rlang::sym(indication_col)
    ) %>%
    dplyr::count() %>%
    tidyr::pivot_wider(
      id_cols = c({{patient_id_col}}, {{indication_col}}, {{sample_id_col}}),
      names_from = {{col_to_compare}},
      values_from = n
    ) %>%
    dplyr::ungroup() %>%
    prepare_counts(#Adapts depending on count_type
      patient_id_col = patient_id_col,
      sample_id_col = sample_id_col,
      indication_col = indication_col,
      count_type = count_type,
      tissue_area_name = tissue_area_name,
      spe = spe
    ) %>%
    tidyr::pivot_longer(
      cols = -c({{patient_id_col}}, {{indication_col}}, {{sample_id_col}}),
      names_to = "feature"
    )

  #Performs tests to estimate whether normal distribution could be assumed
  #Results store in metadata(spe)$normality_check
  if (!skip_normality) {
    spe <- spe %>%
      assert_normal(
        data = cluster_df,
        name = col_to_compare,
        feature = "feature",
        value = "value",
        indication_col = indication_col,
        alpha = shapiro_alpha,
        skew_threshold = skew_threshold,
        kurtosis_range = kurtosis_range
      ) %>%
      assert_equal_var(
        data = cluster_df,
        name = col_to_compare,
        feature = "feature",
        value = "value",
        indication_col = indication_col
      )
  }


  if (clever_test) {
    #Propose a statistical test depending on previous assertion
    spe <- spe %>%
      pick_clever_test(name = col_to_compare) %>%
      perform_clever_test(
        data = cluster_df,
        paired = paired,
        col_to_compare = col_to_compare,
        indication_col = indication_col,
        patient_id_col = patient_id_col,
        adjust = adjust,
        signif_threshold = signif_threshold,
        print_results = !quiet
      )
  } else {
    spe <- spe %>%
      perform_force_test(
        data = cluster_df,
        force_test = force_test,
        paired = paired,
        col_to_compare = col_to_compare,
        indication_col = indication_col,
        patient_id_col = patient_id_col,
        adjust = adjust,
        signif_threshold = signif_threshold,
        print_results = !quiet
      )
  }

  return(spe)
}

# %% [markdown]
# ### prepare_count

# %%
prepare_counts <- function(
  data,
  patient_id_col,
  sample_id_col,
  indication_col,
  count_type,
  tissue_area_name,
  spe
) {
  #designed to be called for a dataframe in the wide format
  #returns the dataframe with the counts adapted to count_type
  suppressPackageStartupMessages(stopifnot(require("dplyr")))

  if (count_type == "count") {
    #Short circuit for raw counts
    return(data)
  } else if (count_type == "density") {
    tissue_area <- S4Vectors::metadata(spe)[[tissue_area_name]]

    data <- data %>%
      tibble::column_to_rownames(sample_id_col) %>%
      .[names(tissue_area), ] %>% #Sorting in the order of tissue_area
      dplyr::mutate(tissue_area = tissue_area) %>% #Add the tissue area data
      dplyr::mutate(
        dplyr::across(#divide all counts by tissue area
          .cols = -dplyr::all_of(c("tissue_area", indication_col, patient_id_col)),
          .fns = ~ .x / tissue_area
        )
      ) %>%
      dplyr::select(-"tissue_area") %>%
      tibble::rownames_to_column(var = sample_id_col) #And clean up the data

    return(data)
  } else if (count_type == "proportion") {
    colnames(data)
    data <- data %>%
      dplyr::mutate(total = rowSums(#add a total column
        dplyr::across(
          -dplyr::all_of(
            c(sample_id_col, patient_id_col, indication_col)
          )
        ),
        na.rm = TRUE
      )) %>%
      dplyr::mutate(
        dplyr::across(#divide the count by that total
          -dplyr::all_of(
            c(patient_id_col, sample_id_col, indication_col)
          ),
          .fns = ~ .x / total
        )
      ) %>%
      dplyr::select(-total) %>%
      tidyr::replace_na(replace = 0)

    return(data)
  } else {
    cli::cli_abort("Error, unrecognized count_type")
  }
}

# %% [markdown]
# ### assert_normal

# %%
assert_normal <- function(
  spe,
  data,
  name = "cluster",
  feature = "feature",
  value = "value",
  indication_col = "indication",
  alpha = 0.05,
  skew_threshold = 1,
  kurtosis_range = c(2, 4)
){
  #Expects a long df, with groups in indication_col,
  #features (cluster) in column feature, value in value
  #WIll perform shapiro test, and compute skewdness and kurtosis and a qq_plot
  #Will also check for equal variance with bartlett test
  # will store the results inside metadata(spe)$normality_check[[name]]
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  suppressPackageStartupMessages(stopifnot(require("rlang")))
  p <- plot_qq(data, feature, indication_col)

  df <- data %>%
    dplyr::group_by(!!sym(feature), !!sym(indication_col)) %>%
    dplyr::summarize(
      shapiro_p = shapiro.test(!!sym(value))$p.value,
      skewness = moments::skewness(!!sym(value)),
      kurtosis = moments::kurtosis(!!sym(value)),
      .groups = "drop"
    )

  spe <- spe %>%
    store_metadata(
      df,
      category = "normality_check",
      name = name,
      type = "normality_data"
    )

  normal_df <- df %>%
    dplyr::group_by(!!sym(feature)) %>%
    dplyr::summarize(
      shapiro_null = all(shapiro_p > alpha, na.rm = TRUE),
      no_skew = all(abs(skewness) < skew_threshold, na.rm = TRUE),
      no_kurt = all(
        dplyr::between(
          kurtosis,
          kurtosis_range[[1]],
          kurtosis_range[[2]]
        ),
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  spe <- spe %>%
    store_metadata(
      normal_df,
      category = "normality_check",
      name = name,
      type = "normality_test"
    )


  normal_df <- normal_df %>%
    dplyr::group_by(!!sym(feature)) %>%
    dplyr::summarize(
      normal_distribution = all(shapiro_null, no_skew, no_kurt)
    )
  
  spe <- spe %>%
    store_metadata(
      normal_df,
      category = "normality_check",
      name = name,
      type = "both_normal"
    ) %>%
    store_metadata(
      p,
      category = "normality_check",
      name = name,
      type = "qq_plot"
    )
  
  return(spe)
}

# %% [markdown]
# #### plot_qq

# %%
plot_qq <- function(data, feature, indication) {
  wrap_formula <- as.formula(paste("~", feature))
  grid_formula <- as.formula(paste(indication, "~", feature))

  data |>
    dplyr::group_by(!!sym(feature)) |>
    dplyr::mutate(sample = row_number()) |>
    ggplot2::ggplot(ggplot2::aes(sample = value)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(color = "red") +
    ggplot2::facet_wrap(wrap_formula, scales = "free") +
    ggplot2::facet_grid(grid_formula) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Q-Q Plots by Feature",
      y = "Sample Quantiles",
      x = "Theoretical Quantiles"
    )
}

# %% [markdown]
# ### assert_equal_var

# %%
assert_equal_var <- function(
  spe,
  data,
  name = "cluster",
  feature = "feature",
  value = "value",
  indication_col = "indication"
) {
  suppressPackageStartupMessages(stopifnot(require("rlang")))
  df <- data |>
    dplyr::group_by(!!sym(feature)) |>
    dplyr::summarize(
      bartlett_p = bartlett.test(!!sym(value) ~ !!sym(indication_col))$p.value,
      equal_var = bartlett_p > 0.05
    )

  spe <- spe |>
    store_metadata(
      df,
      category = "normality_check",
      name = name,
      type = "equal_var"
    )

  return(spe)
}

# %% [markdown]
# ### pick_clever_test

# %%
pick_clever_test <- function(spe, name = "cluster", feature = "feature"){
  both_normal <- spe |>
    get_metadata(
      category = "normality_check",
      name = name,
      type = "both_normal"
    )
  equal_var <- spe |>
    get_metadata(
      category = "normality_check",
      name = name,
      type = "equal_var"
    )

  df <- merge.data.frame(
    both_normal,
    equal_var,
    by.y = feature
  ) %>%
    dplyr::select(-bartlett_p) |>
    dplyr::mutate(test = dplyr::case_when(
      !normal_distribution ~ "wilcox",
      normal_distribution & equal_var ~ "t_test",
      normal_distribution & !equal_var ~ "welch"
    ))

  spe <- spe |>
    store_metadata(
      df,
      category = "normality_check",
      name = name,
      type = "clever_test"
    )

  return(spe)
}

# %% [markdown]
# ### perform_clever_test

# %%
perform_clever_test <- function(
  spe,
  data,
  paired,
  col_to_compare,
  patient_id_col = "patient_id",
  indication_col = "indication",
  feature = "feature",
  value = "value",
  print_results = TRUE,
  adjust = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
  signif_threshold = 0.05
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  suppressPackageStartupMessages(stopifnot(require("rlang")))
  groups <- unique(data[[indication_col]])
  mean_name_1 <- paste0("mean_", groups[[1]])
  mean_name_2 <- paste0("mean_", groups[[2]])

  clever_df <- spe |>
    get_metadata(
      category = "normality_check",
      name = col_to_compare,
      type = "clever_test"
    )

  if (paired) {
    df <- data %>%
      tidyr::pivot_wider(
        id_cols = c(!!sym(patient_id_col), !!sym(feature)),
        names_from = !!sym(indication_col),
        values_from = !!sym(value)
      ) %>%
      dplyr::left_join(clever_df, by = {{feature}}) %>%
      dplyr::group_by(!!sym(feature)) %>%
      dplyr::summarize(
        !!mean_name_1 := mean(!!sym(groups[[1]]), na.rm = TRUE),
        !!mean_name_2 := mean(!!sym(groups[[2]]), na.rm = TRUE),
        p_value = eval(
          call_test(
            x = !!sym(groups[[1]]),
            y = !!sym(groups[[2]]),
            test = test,
            paired = TRUE
          )
        )$p.value,
        performed_test = dplyr::first(test),
        .groups = "drop"
      )
  } else {
    df <- data %>%
      dplyr::left_join(clever_df, by = {{feature}}) %>%
      dplyr::group_by(!!sym(feature)) %>%
      dplyr::summarize(
        !!mean_name_1 := mean(value[indication == groups[[1]]], na.rm = TRUE),
        !!mean_name_2 := mean(value[indication == groups[[2]]], na.rm = TRUE),
        p_value = eval(
          call_test(
            formula = !!sym(value) ~ !!sym(indication_col),
            test = test,
            paired = FALSE
          )
        )$p.value,
        performed_test = dplyr::first(test)
      )
  }

  df <- df %>%
    dplyr::mutate(
      adjusted_pv = p.adjust(p_value, method = adjust),
      signif = adjusted_pv < signif_threshold
    )
  if (print_results) print(df)
  spe %>%
    store_metadata(
      df,
      category = "statistical_test",
      name = col_to_compare,
      type = "clever_test"
    )
}

# %% [markdown]
# ### perform_forced_test

# %%
perform_force_test <- function(
  spe,
  data,
  force_test,
  paired,
  col_to_compare,
  patient_id_col = "patient_id",
  indication_col = "indication",
  feature = "feature",
  value = "value",
  print_results = TRUE,
  adjust = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
  signif_threshold = 0.05
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  suppressPackageStartupMessages(stopifnot(require("rlang")))
  groups <- unique(data[[indication_col]])
  mean_name_1 <- paste0("mean_", groups[[1]])
  mean_name_2 <- paste0("mean_", groups[[2]])

  if (paired) {
    df <- data %>%
      tidyr::pivot_wider(
        id_cols = c(!!sym(patient_id_col), !!sym(feature)),
        names_from = !!sym(indication_col),
        values_from = !!sym(value)
      ) %>%
      dplyr::group_by(!!sym(feature)) %>%
      dplyr::summarize(
        !!mean_name_1 := mean(!!sym(groups[[1]]), na.rm = TRUE),
        !!mean_name_2 := mean(!!sym(groups[[2]]), na.rm = TRUE),
        pvalue = eval(
          call_test(
            x = !!sym(groups[[1]]),
            y = !!sym(groups[[2]]),
            test = force_test,
            paired = TRUE
          )
        )$p.value,
        performed_test = !!force_test,
        .groups = "drop"
      )
  } else {
    df <- data %>%
      dplyr::group_by(!!sym(feature)) %>%
      dplyr::summarize(
        !!mean_name_1 := mean(value[indication == groups[[1]]]),
        !!mean_name_2 := mean(value[indication == groups[[2]]]),
        p_value = eval(
          call_test(
            formula = !!sym(value) ~ !!sym(indication_col),
            test = force_test,
            paired = FALSE
          )
        )$p.value,
        performed_test = !!force_test
      )
  }
  df <- df %>%
    dplyr::mutate(
      adjusted_pv = p.adjust(p_value, method = !!adjust),
      signif = adjusted_pv < !!signif_threshold
    )
  if (print_results) print(df)
  spe %>%
    store_metadata(
      df,
      category = "statistical_test",
      name = col_to_compare,
      type = force_test
    )
}

# %% [markdown]
# #### call_test
# 
# This somewhat convoluted function is used within perform_clever_test to build the expression of the function call for the test.
# It's lack of simplicity is motivated to accomodate for both formula or x and y arguments, which is needed to account for both paired and unpaired test

# %%
call_test <- function(
  formula = NULL,
  x = NULL,
  y = NULL,
  test,
  paired
) {
  #For unpaired, supports by formula and by x
  #For paired, only by x/y is supported
  test <- validate_call_test(formula, x, y, paired, test)

  if (!is.null(formula) && !paired) {
    if(test == "t_test"){
      call <-  rlang::expr(
        t.test(formula = !!formula, var.equal = TRUE)
      )
    } else if (test == "welch") {
      call <- rlang::expr(
        t.test(formula = !!formula, var.equal = FALSE)
      )
    } else if (test == "wilcox") {
      call <- rlang::expr(
        wilcox.test(formula = !!formula)
      )
    }
  } else if (!is.null(x)) {

    if (test == "t_test"){
      call <-  rlang::expr(
        t.test(x = !!x, y = !!y, var.equal = TRUE, paired = !!paired)
      )
    } else if (test == "welch") {
      call <- rlang::expr(
        t.test(x = !!x, y = !!y, var.equal = FALSE, paired = !!paired)
      )
    } else if (test == "wilcox"){
      call <- rlang::expr(
        wilcox.test(x = !!x, y = !!y, paired = !!paired)
      )
    }
  }
  return(call)
}

# %% [markdown]
# # Integrate supervised

# %% [markdown]
# ## integrate_supervised

# %%
integrate_supervised <- function(
  spe,
  folder,
  file_extension = c(".csv", ".tsv"),
  colData_name = "supervised_celltype",
  col_to_select = dplyr::everything(),
  patient_id_col = "patient_id",
  filter_name = "OmiqFilter",
  quiet = FALSE
) {
  #col_to_select : an expression compatible with tidyselect
  #See dplyr_tidy_select for more information
  #Of note, select happen after removing filter_name from the name of the column
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  file_extension <- match.arg(file_extension)

  files <- list.files(folder, full.names = TRUE)
  file_names <- list.files(folder) %>% stringr::str_split_i(file_extension, 1)

  if (!all(
    file_names %in% unique(SummarizedExperiment::colData(spe)[[patient_id_col]])
  )) {
    missing <- file_names %>%
      .[!file_names %in% unique(SummarizedExperiment::colData(spe)[[patient_id_col]])] %>% #nolint
      paste0(file_extension)

    cli::cli_abort(c(
      "Error, some files do not match a patient_id",
      "x" = "All files in folder must be named after one patient_id",
      "i" = "The following files could not be matched : {missing}"
    ))
  }

  if (file_extension == ".csv") {
    read <- function(file) readr::read_csv(file, show_col_types = FALSE)
  } else if (file_extension == ".tsv") {
    read <- function(file) readr::read_tsv(file, show_col_types = FALSE)
  }

  #Initialize a new column if it did not exist before
  if (!colData_name %in% names(SummarizedExperiment::colData(spe))) {
    SummarizedExperiment::colData(spe)[[colData_name]] <- "unknown"
  }

  spe <- purrr::reduce2(
    files,
    file_names,
    function(spe, file, file_name) {
      if (!quiet) message("Reading ", file_name)
      data <- read(file)

      spe <- assign_supervised(
        spe,
        data,
        patient_id = file_name,
        patient_id_col = patient_id_col,
        colData_name = colData_name,
        col_to_select = col_to_select,
        quiet = quiet,
        filter_name = filter_name
      )
    },
    .init = spe
  )

  return(spe)
}

# %% [markdown]
# ### assign_supervised
# 
# Function that takes a SpatialExperiment (or SingleCellExperiment) object and a dataframe corresponding to a supervised gating for one patient, and return the SpatialExperiment object with that supervised gating integrated in the colData_name colData.
# 
# **ARGUMENTS**
# - `spe` : SpatialExperiment (or SingleCellExperiment) object
# - `data` : a dataframe corresponding to the gating for one patient. It is expexted that each cell is a row, each column a filter (for example B cells, CD4+ T cells, etc.), and the value recorded as 1 or 0, 1 indicating that the cell matches that filter
# - `keep_gate` : character vector indicating the colnames of the column (gating filters) to keep, after removing filter_name. See details
# - `col_to_select` : expression to select which column to keep using dplyr::select. See details.
# - `patient_id` : character string, the name of the considered patient, as recorded in the patient_id_col
# - `patient_id_col` : character string, name of the colData(spe) that records the names of the patient/sample
# - `colData_name` : character string, name of the colData where the results will be stored
# - `exclusive` : logical, whether to consider the gates exclusive. See details.
# - `filter_name` : character string indicating which part of the column names to remove, see details.
# - `quiet` : logical, whether to hide or display messages
# 
# **DETAILS**
# 
# *Filter_name* : Some software, like OMIQ, add a constant part to the name of the column of the table it creates (for example, OMIQ might create a table where the column are "OmiqFilter CD4+ T cells"), which just creates clutter while selecting column. <br>
# To help with this, this function proposes to remove that part using string_split_i(), keeping only what is after the first filter_name string. <br>
# If nothing is to be removed, set `filter_name = ""`
# 
# *Selecting features* : For selecting the column of data, col_to_select has priority over keep_gate.  <br>
# Technically, only col_to_select if used by calling dplyr::select(!!col_to_select), but if col_to_select is null (default), it will default to expr(all_of(keep_gate)), providing a more user friendly way to select column for users unfamiliars with expression, while retaining the flexibility of unquoting expressions of advanced users.
# 
# *Exclusive gates* : if gates are considered exclusive, any cell that has been filtered in more than one gate will be recorded as "unknown", same as those that were filtered by no gates. <br>
# If gates are not exclusive, cells filtered by more than on gate will be recorded with all the relevant gates, concatenated as a single string separate with ",". 

# %%
assign_supervised <- function(
  spe,
  data,
  keep_gate = NULL,
  col_to_select = NULL,
  patient_id,
  patient_id_col = "patient_id",
  colData_name,
  exclusive = TRUE,
  filter_name = "OmiqFilter",
  quiet = FALSE
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  #Exclusive transforms double positive cells into "unknown"
  #If false, both name are made into a single string, separated by ","
  data <- data %>%
    dplyr::select(grep(filter_name, colnames(data)))

  colnames(data) <- colnames(data) %>%
    stringr::str_split_i(filter_name, 2)

  if (is.null(col_to_select)) {
    col_to_select <- rlang::expr(
      dplyr::all_of(keep_gate)
    )
  }

  if (!quiet) message("Preparing data for integration")
  celltypes <- data %>%
    dplyr::select(!!col_to_select) %>%
    dplyr::rowwise() %>%
    dplyr::summarize(
      celltype = paste(
        names(.)[which( #get the name
          dplyr::c_across(tidyselect::everything()) == 1
        )], #For which value = 1
        collapse = "," #collapse them as a single string
      )
    ) %>%
    dplyr::pull(celltype)

  #handle "" and double values
  if(exclusive){
    celltypes[grepl(",", celltypes) | celltypes == ""] <- "unknown"
  } else {
    celltypes[celltypes == ""] <- "unknown"
  }

  #add to coldata
  spe <- modify_subset_colData(
    spe,
    colData_name = colData_name,
    subset_by = patient_id_col,
    subset_value = patient_id,
    vector = celltypes,
    quiet = quiet
  )

  return(spe)
}

# %% [markdown]
# #### modify_subset_colData

# %%
modify_subset_colData <- function(
  spe,
  colData_name,
  subset_by,
  subset_value,
  vector,
  quiet
) {
  #Helper to integrate colData based on a subset
  #looks for colData(spe)[[subset_by]] %in% subset_value
  #The vector needs to be already for the subset only
  #return the spe with the colData changed
  spe_subset <- rlang::expr(
    SummarizedExperiment::colData(spe)[[!!subset_by]] %in% !!subset_value
  )
  if (!quiet) message("Integrating data into the object")

  SummarizedExperiment::colData(spe[, eval(spe_subset)])[[colData_name]] <- spe %>% #nolint
    SummarizedExperiment::colData() %>%
    as.data.frame %>%
    dplyr::filter(!!rlang::sym(subset_by) %in% subset_value) %>%
    dplyr::mutate(!!colData_name := vector) %>%
    dplyr::pull(!!colData_name)

  return(spe)
}

# %% [markdown]
# # Classifier

# %% [markdown]
# ## train_model

# %%
train_model <- function(
  spe,
  colData_name,
  model_name = colData_name,
  assay_name = "counts",
  marker_to_use =  rownames(spe)[SummarizedExperiment::rowData(spe)$use_channel], #nolint
  use_reducedDim = FALSE,
  reducedDim_name = "HARMONY",
  unknown_type = c("unknown"),
  training_partition = 0.75,
  ntree = 1000,
  tuneLength = 5,
  color_pal = RColorBrewer::brewer.pal(12, "Paired"),
  ncores = 1,
  quiet = FALSE,
  ...
) {
  #Create a predictive model (like random forest) using caret::train
  #To predict unknow values in a colData, based on a assay or a reducedDim
  #The model is then stored in the metadata(spe)$predictive_model[[model_name]]
  #... is passed to caret::train

  validate_train_model(spe, colData_name, model_name,
                       assay_name, marker_to_use, use_reducedDim,
                       reducedDim_name, ncores)

  if (!quiet) message("Creating dataframe")
  df <- extract_col_data(
    spe,
    colData_name = colData_name,
    use_reducedDim = use_reducedDim,
    reducedDim_name = reducedDim_name,
    assay_name = assay_name,
    marker_to_use = marker_to_use,
    unknown_type = unknown_type,
    assign_marker_to_use = TRUE,
    env = parent.frame()
  )

  #train_index will be stored at the end
  #to be able to find the cells used for training
  train_index <- df %>%
    dplyr::filter(!(!!sym(colData_name) %in% unknown_type)) %>%
    dplyr::pull(!!colData_name) %>%
    as.factor() %>%
    caret::createDataPartition(p = training_partition) %>%
    .[[1]] #Extract the vector from the list


  df_train <- df %>%
    dplyr::filter(!(!!sym(colData_name) %in% unknown_type)) %>%
    dplyr::slice(train_index)

  if (!quiet) message("Starting training")

  cluster <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cluster)

  train_control <- trainControl(
    method = "cv",
    number = 5,
    allowParallel = TRUE,
    verboseIter = !quiet
  )

  model <- df_train %>%
    dplyr::select(dplyr::all_of(marker_to_use))  %>%
    caret::train(
      y = factor(df_train[[colData_name]]),
      method = "parRF",
      ntree = ntree,
      tuneLength = tuneLength,
      trControl = train_control
    )

  parallel::stopCluster(cluster)
  foreach::registerDoSEQ()

  #Integrate the metadata
  spe <- store_multiple_metadata(
    spe,
    data = list(model, train_index),
    category = "predictive_model",
    name = model_name,
    type = c("model", "train_index")
  ) 

  spe <- evaluate_model(
    spe,
    model_name,
    df_train,
    marker_to_use,
    colData_name,
    color_pal
  )

  return(spe)
}

# %% [markdown]
# ### evaluate_model

# %%
evaluate_model <- function(
  spe,
  model_name,
  training_df,
  marker_to_use,
  colData_name,
  color_pal
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))

  model <- spe %>%
    get_metadata(
      category = "predictive_model",
      name = model_name,
      type = "model"
    )
  #Curve of accuaracy depending on number of predictor
  accuracy_curve <- ggplot2::ggplot(model) +
    ggplot2::geom_errorbar(
      data = model$results,
      aes(
        ymin = Accuracy - AccuracySD,
        ymax = Accuracy + AccuracySD
      ),
      width = 0.4
    ) +
    ggplot2::theme_classic(base_size = 15)
  #Plot of importance of variance per feature
  variance_plot <- caret::varImp(model)

  # Predict the cell phenotype labels of the test data
  # to use for confusion matrix
  train_data <- training_df %>%
    dplyr::select(dplyr::all_of(marker_to_use))

  predicted_type <- model %>%
    predict(newdata = train_data)
  
  prediction <- model %>%
    predict(newdata = train_data, type = "prob")

  confusion_matrix <- predicted_type %>%
    caret::confusionMatrix(
      reference = factor(training_df[[colData_name]]),
      mode = "everything"
    )

  sensitivity_plot <- confusion_matrix %>%
    .$byClass %>%
    data.frame() %>%
    dplyr::mutate(
      class = sub("Class: ", "", rownames(confusion_matrix$byClass))
    ) %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(
      aes(
        1 - Specificity,
        Sensitivity,
        size = Detection.Rate,
        fill = class
      ),
      shape = 21
    ) +
    ggplot2::scale_fill_manual(values = color_pal) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::ylab("Sensitivity (TPR)") +
    ggplot2::xlab("1 - Specificity (FPR)")

  prediction$truth <- factor(training_df[[colData_name]])

  class_prediction_plot <- prediction %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(unique(training_df[[colData_name]]))
    ) %>%
    ggplot2::ggplot() +
    ggplot2::geom_boxplot(
      aes(x = name, y = value, fill = name),
      outlier.size = 0.5
    ) +
    ggplot2::facet_wrap(. ~ truth, ncol = 1) +
    ggplot2::scale_fill_manual(values = color_pal)  +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )


  spe <- spe %>%
    store_multiple_metadata(
      data = list(accuracy_curve, variance_plot, confusion_matrix,
                  sensitivity_plot, class_prediction_plot),
      category = "predictive_model",
      name = model_name,
      type = list("accuracy_curve", "variance_plot", "confusion_matrix",
                  "sensitivity_plot", "class_prediction_plot")
    )

  return(spe)
}

# %% [markdown]
# ## classify_new_cells

# %%
classify_new_cells <- function(
  spe,
  model_name,
  colData_name,
  model = NULL,
  probability_threshold = 0.6,
  predicted_name = paste0("predicted_", colData_name),
  assay_name = "counts",
  marker_to_use =  rownames(spe)[SummarizedExperiment::rowData(spe)$use_channel], #nolint
  use_reducedDim = FALSE,
  reducedDim_name = "HARMONY",
  unknown_type = c("unknown"),
  color_pal = RColorBrewer::brewer.pal(12, "Paired")
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  validate_classify_cell(
    spe, model_name, colData_name,
    model, marker_to_use, use_reducedDim,
    reducedDim_name, unknown_type
  )

  data <- extract_col_data(
    spe,
    colData_name = colData_name,
    use_reducedDim = use_reducedDim,
    reducedDim_name = reducedDim_name,
    assay_name = assay_name,
    marker_to_use = marker_to_use,
    unknown_type = unknown_type,
    assign_marker_to_use = TRUE,
    env = parent.frame()
  )

  unknown_data <- data %>%
    dplyr::filter(!!rlang::sym(colData_name) %in% unknown_type)

  prediction <- predict(
    model,
    newdata = unknown_data,
    type = "prob"
  ) %>%
    dplyr::mutate(
      max_prob = do.call(pmax, .),
      predicted_class = colnames(.)[max.col(., ties.method = "first")]
    ) %>%
    dplyr::mutate( #add object_id as the reference
      object_id = unknown_data$object_id,
      .before = 1
    )

  #Merge the prediction, using object_id as a reference
  data <- data %>%
    dplyr::left_join(prediction, by = "object_id") %>%
    dplyr::mutate( #replaces NA (non predicted) cells with previous knowledge
      predicted_class = dplyr::coalesce(
        predicted_class,
        !!rlang::sym(colData_name)
      )
    ) %>%
    dplyr::mutate( #Handles below probability threshold values
      predicted_class = dplyr::case_when(
        max_prob < !!probability_threshold ~ "unknown",
        .default = predicted_class
      )
    )
  prediction_plot <- plot_prediction_ridgeplot(prediction, color_pal)
  # Integrate into spe using object_id (rownames of colData(spe))
  # To avoid shifting cells around
  SummarizedExperiment::colData(spe)[[predicted_name]] <- spe %>%
    SummarizedExperiment::colData() %>%
    tibble::as_tibble(rownames = "object_id") %>%
    dplyr::select(object_id, !!colData_name) %>%
    dplyr::left_join(
      dplyr::select(data, object_id, predicted_class), 
      by = "object_id"
    ) %>%
    dplyr::pull(predicted_class)


  spe <- spe %>%
    store_multiple_metadata(
      list(prediction_plot, prediction),
      category = "predictive_model",
      name = model_name,
      type = c("prediction_plot", "prediction_df")
    )

  return(spe)
}

# %% [markdown]
# ### plot_prediction_ridgeplot

# %%
plot_prediction_ridgeplot <- function(prediction, color_pal) {
  #Helper function to do one simple plot
  prediction |>
    ggplot2::ggplot() +
    ggridges::geom_density_ridges(
      ggplot2::aes(x = max_prob, y = predicted_class, fill = predicted_class)
    ) +
    ggplot2::scale_fill_manual(values = color_pal) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::xlab("Maximum probability") +
    ggplot2::ylab("Cell type") +
    ggplot2::xlim(c(0, 1.2))
}

# %% [markdown]
# # Spatial

# %% [markdown]
# ## neighborhood_analysis
# 
# Wrapper aroung aggregate_neighborhood and kmeans clustering. K_means is vectorized to allow to only perform the aggregate_neighborhood only once. 
# 
# **ARGUMENTS**
# - `colPair` : character string, name of the colPair(spe) that contains the neighbor graph
# - `name` : name of the colData that will contain the DataFrame of the aggregated Neighbor
# - `aggregate_by` : either "metadata" to aggregate based on a colData (like celltype), "expression" or "reducedDim". See aggregate_integrated_neighbors for more information
# - `count_by` : for "metadata", name of the colData to count
# - `assay_name` : for "expression", name of the assay to use for aggregation
# - `subset_row` : for "expression" : passed to aggregateNeighbors
# - `reduced_dim_name` : for "reducedDIm", name of the reducedDim(spe) to use
# - `statistic` : character value. Possible values are "mean", "median", "sd" and "var". Passed to aggregate_integrate_neighbors to choose the type of summary statistic to aggregate
# - `k_value` : integer, number of centers for kmeans clustering. Supports integer vectors, and one kmeans will be performed and recorded for each value of k_value.
# - `kmean_name` : character, name of the kmean colData. If k_value is a vector, you can provide a vector of equal length for the names. If only one name is provided, will attempt to recycle it in the format "{kmean_name}_{k_value}"
# - `...` : additional arguments passed to kmeans
# 
# **VALUE**
# - the aggregated neighbors in colData(spe)[[neighbor_name]]
# - a clustering of neighboorhoods in colData(spe)[[kmean_name]]
# 
# #TODO : argument checker

# %%
perform_neighborhood_analysis <- function(
  spe,
  colPairName = "knn",
  neighbor_name = "neigborhood",
  aggregate_by = c("metadata", "expression", "reducedDim"),
  count_by = "celltype",
  assay_name = "counts",
  subset_row = NULL,
  reduced_dim_name = NULL,
  statistic = c("mean", "median", "sd", "var"),
  k_value = c(6),
  kmean_name = c("neighborhood_clusters"),
  ...
) {
  validate_perform_neighborhood(
    spe, colPairName, neighbor_name, aggregate_by, count_by, assay_name,
    subset_row, reduced_dim_name, statistic, k_value, kmean_name
  )
  #Using the slightly modify aggregateNeighbors, to also allow aggregation over reducedDim
  spe <- aggregate_integrated_neighbors(
    spe,
    colPairName = colPairName,
    aggregate_by = aggregate_by,
    count_by = count_by,
    assay_name = assay_name,
    subset_row = subset_row,
    reduced_dim_name = reduced_dim_name,
    statistic = statistic,
    name = neighbor_name
  )
  #attempt recycling kmean_name
  if (length(k_value) > 1 &&
        length(kmean_name) == 1) {
    kmean_name <- paste0(
      kmean_name, "_", k_value
    )
  }

  spe <- purrr::reduce2(
    k_value,
    kmean_name,
    function(spe, x, y) {
      #handle NA values from the knn (cells without neighbors)
      mat <- as.matrix(SummarizedExperiment::colData(spe)[[neighbor_name]])
      mat[is.na(mat)] <- 0
      
      cn <- kmeans(
        mat,
        centers = x,
        ...
      )
      #integrate into spe
      SummarizedExperiment::colData(spe)[[y]] <- as.factor(cn$cluster)

      spe
    },
    .init = spe
  )

  return(spe)
}

# %% [markdown]
# ### aggregate_integrated_neighbors
# 
# This function is a slightly modified version of imcRtools::aggregateNeighbors that also provide the possibility to aggregate by reducedDim. <br>
# This is usefull to aggregate integrated data such as after batch effect correction (by HARMONY or other). <br>
# This function is only a marginal modification of the original code of imcRtools::aggregateNeighbors, simply changing the intial input to accomodate reducedDim. All credits to the authors of that function.

# %%
aggregate_integrated_neighbors <- function(
  object,
  colPairName,
  aggregate_by = c("metadata", "expression", "reducedDim"), 
  count_by = NULL,
  proportions = TRUE,
  assay_name = NULL,
  subset_row = NULL,
  reduced_dim_name = NULL,
  statistic = c("mean", "median", "sd", "var"),
  name = NULL)
{
  aggregate_by <- match.arg(aggregate_by)
  
  if(aggregate_by %in% c("metadata", "expression")){
    #Delegates to aggregateNeighbors what can be done with that function
    object <- imcRtools::aggregateNeighbors(
      object = object,
      colPairName = colPairName,
      aggregate_by = aggregate_by,
      count_by = count_by,
      proportions = proportions,
      assay_type = assay_name,
      subset_row = subset_row,
      statistic = statistic,
      name = name
    )
    return(object)
  } else {
    object |> verify_reducedDim(reduced_dim_name)
    
    cur_dat <- data.table::as.data.table(
      SingleCellExperiment::colPair(object, colPairName)
    )
    #Only change compared to aggregateNeighbors : making it compatible with dimension reduction such as HARMONY
    #Add the dimension reduction to the data.table
    #Unlike the assay, it does not need transposition (natively row = observation)
    cur_dat <- cbind(
      cur_dat,
      SingleCellExperiment::reducedDim(object, reduced_dim_name)[cur_dat$to, ]
    )

    cur_dat <- melt(cur_dat, id.vars = c("from", "to"))
    #Performs the statistic (such as mean(value)), 
    #grouped by from (cell) and variable (value in the space of dimension reduction)
    cur_dat <- cur_dat[ ,
                       eval(parse(text = paste0(statistic, "(value)"))),
                       by = c("from", "variable")]

    cur_dat <- data.table::dcast(
      cur_dat,
      formula = "from ~ variable",
      value.var = "V1"
    )

    name <- ifelse(is.null(name),
                   paste0(statistic, "_aggregatedExpression"), name)

    out_dat <- DataFrame(matrix(data = NA,
                                nrow = ncol(object),
                                ncol = ncol(cur_dat) - 1))
    names(out_dat) <- names(cur_dat)[-1]
    out_dat[cur_dat$from, ] <- cur_dat[,-1]

    SummarizedExperiment::colData(object)[[name]] <- out_dat

    return(object)
  }
}

# %% [markdown]
# ## Patch detection

# %% [markdown]
# ### plot_fraction_patch :
# 
# Used to create a plot that represent the proportion of given celltype in patches defined by Patch detection (see imcRtools)
# 
# **ARGUMENTS**
# - `spe` : a SpatialExperiment object
# - `cell_type` : character string, name of the colData (cluster) to represent
# - `chosen_type` : vector (usually character or integer) indicating the types within cell_type to represent
# - `chosen_name` : character string, used for naming the label of the plot
# - `patch_name` : character string, name of the colData containing the patch ids
# - `patient_id` : character string, name of the colData containing the patient/sample ids
# - `color_by` : character string, name of colData that will guide coloring
# 
# **VALUE**
# Return a ggplot object

# %%
plot_fraction_patch <- function(
  spe,
  cell_type,
  chosen_type,
  chosen_name = chosen_type[[1]],
  patch_name = "immune_patch",
  patient_id = "patient_id",
  color_by = "indication"
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  suppressPackageStartupMessages(stopifnot(require("rlang")))
  verify_multiple_colData(spe, c(cell_type, patch_name, patient_id, color_by))

  #Separate to avoid duplicating a column :
  if (color_by != patient_id && color_by != patch_name) {
    df <- SummarizedExperiment::colData(spe) %>%
      tibble::as_tibble(rownames = "object_id") %>%
      dplyr::group_by(
        !!sym(patch_name),
        !!sym(patient_id),
        !!sym(color_by)
      )
  } else {
    df <- SummarizedExperiment::colData(spe) %>% 
      tibble::as_tibble(rownames = "object_id") %>%
      dplyr::group_by(
        !!sym(patch_name),
        !!sym(patient_id)
      )
  }

  count_name <- paste0(chosen_name, "_count")
  freq_name <- paste0(chosen_name, "_freq")

  df %>%
    dplyr::summarize(
      !!sym(count_name) := sum(!!sym(cell_type) %in% chosen_type),
      patch_size = dplyr::n(),
      !!sym(freq_name) :=  !!sym(count_name)/patch_size,
      .groups = "drop"
    ) %>%
    ggplot2::ggplot() +
      ggplot2::geom_point(
        ggplot2::aes(
          x = log10(patch_size),
          y = !!sym(freq_name),
          color = !!sym(color_by)
        )
      ) +
      ggplot2::theme_classic()
}

# %% [markdown]
# ### plot_expression_patch
# 
# Used to create a plot that represent the mean expression of markers in patches defined by Patch detection (see imcRtools)
# 
# **ARGUMENTS**
# - `spe` : a SpatialExperiment object
# - `marker` : character string(NB : not vectorized) indicating which marker expression to plot
# - `assay_name` : character string, indicationg which assay to retrieve the values from
# - `patch_name` : character string, name of the colData containing the patch ids
# - `img_id` : character string, name of the colData containing the patient/sample ids
# - `color_by` : character string, name of colData that will guide coloring
# 
# **VALUE**
# Return a ggplot object

# %%
plot_expression_patch <- function(
  spe,
  marker,
  assay_name = "counts",
  patch_name = "immune_patch",
  img_id = "patient_id",
  color_by = "indication"
){
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  suppressPackageStartupMessages(stopifnot(require("rlang")))
  spe %>%
    verify_assay(assay_name) %>%
    verify_markers(marker) %>%
    verify_multiple_colData(c(patch_name, img_id, color_by))

  if (color_by != img_id && color_by != patch_name){
    df <- extract_col_data(
      spe,
      colData_name = patch_name,
      assay_name = assay_name,
      additional_col = c(img_id, color_by)
    )
  } else {
    df <- extract_col_data(
      spe,
      colData_name = patch_name,
      assay_name = assay_name,
      additional_col = c(img_id)
    )
  }

  mean_name <- paste0("mean_", marker)

  df %>%
    dplyr::group_by(
      !!rlang::sym(patch_name),
      !!rlang::sym(img_id),
      !!rlang::sym(color_by)
    ) %>%
    dplyr::summarize(
      !!rlang::sym(mean_name) := mean(!!rlang::sym(marker)),
      patch_size = dplyr::n(),
      .groups = "drop"
    ) %>%
    ggplot2::ggplot() +
      ggplot2::geom_point(
        ggplot2::aes(
          x = log10(patch_size),
          y = !!rlang::sym(mean_name),
          color = !!rlang::sym(color_by)
        )
      ) +
      ggplot2::theme_classic()
}

# %% [markdown]
# ## Area and density

# %% [markdown]
# ### compute_patch_area

# %%
compute_patch_area <- function(
  spe,
  region_id,
  img_id = "patient_id",
  colPairName = "expansion",
  spatial_coord = c("Pos_X", "Pos_Y"),
  min_patch_size = 3,
  concavity = 1,
  BPPARAM = BiocParallel::SerialParam(),
  quiet = FALSE
) {
  #returns a dataframe containing the region name and the computed area

  #Get continuous patch of each region, see defin_region_patch
  if(!quiet) message("Defining region patches")
  spe <- define_region_patch(
    spe,
    region_id = region_id,
    colPairName = colPairName,
    img_id = img_id,
    min_patch_size = min_patch_size,
    name = "patch_region"
  )
  
  #We extract the required information in a dataframe :
  if (!quiet) message("Building dataframe")
  df <- SpatialExperiment::spatialCoords(spe) %>%
    tibble::as_tibble()

  df <- SummarizedExperiment::colData(spe) %>%
    tibble::as_tibble(rownames = "object_id") %>%
    dplyr::select(object_id, !!img_id, patch_region) %>%
    cbind(df, .) %>%
    dplyr::filter(!is.na(patch_region))

  # Objective : compute the area covered by the regions in region_id
  # We will iterate over the region_id,
  # over the img_id and over each patch in that image
  # For each patch, we will compute the concave area using concaveman()
  # We will add patches of the same regiontype together across all sample
  # To estimate the total area covered by the region
  region_type <- df$patch_region %>%
    stringr::str_split_i("___", 1) %>%
    unique()

  area_vector <- purrr::map_dbl(
    region_type,
    function(region) { #region will be a region type
      region_area <- 0

      region_area <- purrr::reduce(
        unique(df[[img_id]]),
        function(region_area, sample) {
          # Extract coordinates of each patches of the corresponding region
          # For one sample
          sample_df <- df %>%
            dplyr::filter(
              !!rlang::sym(img_id) == sample,
              stringr::str_split_i(patch_region, "___", 1) == region
            )

          # Escape cases with no cell
          if (nrow(sample_df) != 0) {
            # We will extract each individual patch, 
            # Estimate a concave polygon of that patch
            # And compute the area of this polygon

            patch_areas <- BiocParallel::bplapply(
              unique(sample_df$patch_region),
              compute_polygon_area,
              BPPARAM = BPPARAM,
              spatial_coord = spatial_coord,
              sample_df = sample_df,
              concavity = concavity
            )

            sample_area <- sum(unlist(patch_areas))
          } else {
            #No cells means area = 0
            sample_area <- 0
          }

          #Add the area of each patient to the final results
          region_area <- region_area + sample_area

          return(region_area)
        },
        .init = region_area
      )
      return(region_area)
    }
  )
  #Finally prepare a dataframe with the results and return it
  returning_df <- data.frame(
    region_type = stringr::str_split_i(region_type, "patch_", 2),
    area = area_vector
  )
  return(returning_df)
}

# %% [markdown]
# #### define_region_patch

# %%
define_region_patch <- function(
  spe,
  region_id,
  colPairName,
  img_id,
  min_patch_size = 3,
  name = "patch_region"
) {
  spe <- purrr::reduce(
    unique(spe[[region_id]]),
    function(spe, x) {
      spe <- imcRtools::patchDetection(
        spe,
        patch_cells = spe[[region_id]] == x,
        colPairName = colPairName,
        min_patch_size = min_patch_size,
        name = paste0("patch_", x),
        img_id = img_id
      )
    },
    .init = spe
  )

  # Now we create a vector that combine each patch identification
  # And we record both the region type and the patch id
  # Separated by a special"___" separator that we will use
  # To extract each region_type separatly
  spe[[name]] <-
    SummarizedExperiment::colData(spe) %>%
    tibble::as_tibble(rownames = "object_id") %>%
    dplyr::select(
      dplyr::all_of(
        paste0(
          "patch_",
          unique(spe[[region_id]])
        )
      )
    ) %>%
    dplyr::mutate(
      patch_region = purrr::pmap_chr(
        .,
        function(...) {
          vals <- list(...)
          if(any(!is.na(unlist(vals)))) {
            paste0(
              names(vals)[!is.na(unlist(vals))], #colname = region
              "___", 
              vals[!is.na(unlist(vals))] #value = number of the patch
            )
          } else {
            NA_character_
          }
        }
      )
    ) %>%
    dplyr::pull(patch_region)

  #Free the unwanted colData
  #Returning only the coalesced column
  spe <- purrr::reduce(
    paste0("patch_", unique(spe[[region_id]])),
    function(spe, x){
      spe[[x]] <- NULL

      return(spe)
    },
    .init = spe
  )

  return(spe)

}


# %% [markdown]
# #### compute_polygon_area

# %%
compute_polygon_area <- function(
  patch,
  sample_df,
  spatial_coord,
  concavity
) {
  #Function designed to be used within bplapply
  # Using base pipe here to avoid having to attach dplyr inside each cluster
  # Honnestly should always use bas pipe as I don't care about not being compatible with older R version
  if (sample_df |>
        dplyr::filter(patch_region == patch) |>
        nrow() != 0) {
    sample_df |>
      dplyr::filter(patch_region == patch) |>
      sf::st_as_sf(coords = spatial_coord) |>
      concaveman::concaveman(concavity = concavity) |>
      sf::st_area()
  } else {
    return(0)
  }
}

# %% [markdown]
# # Troubleshooting functions

# %% [markdown]
# ## review_tag_names
# 
# Function used to see review the exact name in the panel in the files to be prepared (in prepare_csv). Will only review the first file.
# 
# **ARGUMENTS**
# - `in_dir` : Character string, path to the folder containing the files to review
# - `file_type` : Character string, file extension. Supported ".csv" and ".tsv"
# - `cel_part` : Character string, used to subset the columns for the panel. Expected to be cellpart considered (like "Entire.cell") if expression of multiple cell part exists in the input files
# - `exclude` : Character string, used to exclude some column
# 
# **VALUE** <br>
# Returns the colnames that match cell_part and are not in exclude, usually the colnames corresponding to the panel.

# %%
review_tag_names <- function(in_dir,
                             file_type,
                             cell_part,
                             exclude) {
  files <- list.files(path = in_dir[[1]], pattern = file_type)
  #Read the files in df
  if (file_type == ".csv") {
    df <- read.delim(paste0(in_dir[[1]], files[1]),
                     sep = ",")
  } else if(file_type == ".tsv") {
    df <- read.delim(paste0(in_dir[[1]], files[1]),
                     sep = "\t")
  } else {#returns a error if unsupported extension
    stop("Error = File extension not supported : Supported pattern '.csv' ou '.tsv'") # nolint: line_length_linter.
  }

  chosen <- grepl(cell_part,
                  colnames(df))
  if (!is.null(exclude)) {
    chosen <- chosen & !grepl(pattern = exclude, colnames(df))
  }

  colnames(df[, chosen])
}

# %% [markdown]
# ## review_intensity_column
# 
# Function used to verify if the provided tag vector and panel vector match.
# 
# **ARGUMENTS**
# - `in_dir` : character string, path to the folder containing the files
# - `file_type` : character string, extension of the files. Supported extension are ".tsv" and ".csv"
# - `tag` : character vector, containing the tag (metal, fluorescent, ec.) names
# - `cell_part` : character vector,  used to subset the columns for the panel. Expected to be cellpart considered (like "Entire.cell") if expression of multiple cell part exists in the input files
# - `panel` : character vector, used to indentify the panel (the markers recognized by the antibodies)
# - `exclude` : character string, used to exclude some column
# 
# **VALUE** <br>
# Returns the colnames found based on the tag vector. <br>
# Also prints the number of column identified respectively by tag and panel (which should match), and the colnames of those that don't perfectly match one to one.

# %%
review_intensity_column <- function(in_dir,
                                    file_type,
                                    tag,
                                    cell_part,
                                    panel,
                                    exclude) {
  files <- list.files(path = in_dir[[1]], pattern = file_type)
  #Read the files in df
  if (file_type == ".csv") {
    df <- read.delim(paste0(in_dir[[1]], files[1]),
                     sep = ",")
  } else if (file_type == ".tsv") {
    df <- read.delim(paste0(in_dir[[1]], files[1]),
                     sep = "\t")
  } else {#returns a error if unsupported extension
    stop("Error = File extension not supported : Supported pattern '.csv' ou '.tsv'") # nolint: line_length_linter.
  }
  #Finding the column
  col_tag <- map(tag, function(x) {

    correct <- grepl(pattern = x, colnames(df)) &
      grepl(pattern = cell_part, colnames(df))

    if (!is.null(exclude)) {
      correct <- correct & !grepl(pattern = exclude, colnames(df))
    }
    correct_col <- colnames(df)[correct]
    return(correct_col)
  })

  #Checking the number of column
  col_length <- map_int(col_tag, function(x) length(x))
  if (any(col_length != 1)) {
    cat("The following element of tag did not return a singular column:\n")
    cat(unlist(col_tag[col_length != 1]), sep = "\n")
  }

  cat("Length of tag : ", length(tag), "\n")
  cat("Length of panel : ", length(panel), "\n")
  cat("Length of columns : ", length(unlist(col_tag)), "\n")

  unlist(col_tag)
}

# %% [markdown]
# ## review_tag_column
# 
# Function used to review the name of the column that do not contain the tag, to identify those which contains the regionprops (x, y, area, etc.) data.
# 
# **ARGUMENTS**
# - `in_dir` : character string, path to the folder containing the files
# - `file_type` : character string, extension of the files. Supported extension are ".tsv" and ".csv"
# - `tag` : character vector, containing the tag (metal, fluorescent, ec.) names
# 
# **VALUE**<br>
# Returns all the colnames that don't match the tag.

# %%
review_nontag_column <- function(in_dir,
                                file_type,
                                tag) {
  #We assum the names of the columns to be the same across all files
  files <- list.files(path = in_dir[[1]], pattern = file_type)
  #Read the files in df
  if (file_type == ".csv") {
    df <- read.delim(paste0(in_dir[[1]], files[1]),
                     sep = ",")
  } else if (file_type == ".tsv") {
    df <- read.delim(paste0(in_dir[[1]], files[1]),
                     sep = "\t")
  } else {#returns a error if unsupported extension
    stop("Error = File extension not supported : Supported pattern '.csv' ou '.tsv'") # nolint: line_length_linter.
  }

  correct_tag <- lapply(tag, function(x) {
    #We keep the columns that :
    #1 correspond to the element x of the tag
    #2 correspond to the selected cellpart subset
    correct <- grepl(pattern = x, colnames(df))
    return(correct)
  })

  correct <- rep(FALSE, length(correct_tag[[1]]))
  for (i in seq_along(correct_tag)){
    correct <- correct | correct_tag[[i]]
  }

  colnames(df[, !correct])
}

# %% [markdown]
# ## review_regionprops_column
# 
# Function used to check the match between the keep_prop vector and the associated colnames.
# 
# **ARGUMENTS**
# - `in_dir` : character string, path to the folder containing the files
# - `file_type` : character string, extension of the files. Supported extension are ".tsv" and ".csv"
# - `keep_prop` : character vector, containing the names of the regionprops feature to extract.
# 
# **VALUE**
# Returns the colnames found based on the keep_prop vector. <br>
# Also prints the number of column identified by keep_prop and which value did not identify one single column.

# %%
review_regionprops_column <- function(in_dir,
                                      file_type,
                                      keep_prop) {
  files <- list.files(path = in_dir[[1]], pattern = file_type)
  #Read the files in df
  if (file_type == ".csv") {
    df <- read.delim(paste0(in_dir[[1]], files[1]),
                     sep = ",")
  } else if (file_type == ".tsv") {
    df <- read.delim(paste0(in_dir[[1]], files[1]),
                     sep = "\t")
  } else {#returns a error if unsupported extension
    stop("Error = File extension not supported : Supported pattern '.csv' ou '.tsv'") # nolint: line_length_linter.
  }
  #Finding the column
  col_props <- lapply(keep_prop, function(x) {
    #We keep the columns that :
    #1 correspond to the element x of the tag
    #2 correspond to the selected cellpart subset
    correct <- grepl(pattern = x, colnames(df))
    correct_col <- colnames(df)[correct]
    return(correct_col)
  })

  #Checking the number of column
  col_length <- unlist(lapply(col_props, function(x) length(x)))
  if (any(col_length != 1)) {
    cat("The following element of tag did not return a singular column:\n")
    cat(col_props[col_length != 1], sep = "\n")
  }

  cat("Length of keep_props : ", length(keep_prop), "\n")
  cat("Length of columns : ", length(unlist(col_props)), "\n")

  unlist(col_props)
}

# %% [markdown]
# # Utils

# %% [markdown]
# ## metadata

# %% [markdown]
# ### store_metadata
# 
# The metadata slot will be organized in  3 layers of lists
# - category
# - name
# - type
# 
# Is included a builder for this organization, and an acessor. <br>
# The builder, store_metadata is pipe compatible, and also exists in a vectorized form, store_multiple_metadata,which can store multiple types of the same category/name at once, which is a common pattern. <br>
# tje accessor get_metadata can also returns category level or name level if the other argument are empty. <br>
# 
# Not all metadata will follow this structure, for example clustering_methods will likely 

# %%
store_metadata <- function(spe, data, category, name, type = NULL) {
  #Simple helper to encapsulate storing metadata
  #WIth a 3 level list :
  #category such as "statistic" or "normality test"
  #name : name of the cluster associated with the metadata
  #Type : the data that is stored
  if (!category %in% names(S4Vectors::metadata(spe))){
    S4Vectors::metadata(spe)[[category]] <- list()
  }

  if (!name %in% names(S4Vectors::metadata(spe)[[category]])) {
    S4Vectors::metadata(spe)[[category]][[name]] <- list()
  }
  #Allow storing at the name level also
  #Can be shorter than using store_multiple_metadata
  if (is.null(type)) {
    S4Vectors::metadata(spe)[[category]][[name]] <- data
  }else {
    S4Vectors::metadata(spe)[[category]][[name]][[type]] <- data
  }

  return(spe)
}

# %%
store_multiple_metadata <- function(spe, data, category, name, type) {
  # Simple helper to encapsulate storing metadata
  # this version is vectorized and can assign a list of object, provided a ,
  # which can be easier when storing multiple metadata
  # Only work for lists of object to store, else will store the entire object
  # A little more overhead initially, but should be faster for multiple metadata
  #Of note, it will store elements of the list, 
  #so if your data are litteral lists, please provide data = list(my_list1, my_list2)
  if (length(data) != length(type)) {
    cli::cli_abort(c(
      "Incorrect length of data and type when storing multiple metadata",
      "x" = "Length of data and type argument must be the same",
      "i" = "Length of data {length(data)}, length of type {length(type)}"
    ))
  }

  if (!category %in% names(S4Vectors::metadata(spe))) {
    S4Vectors::metadata(spe)[[category]] <- list()
  }

  if (!name %in% names(S4Vectors::metadata(spe)[[category]])) {
    S4Vectors::metadata(spe)[[category]][[name]] <- list()
  }
  if (is.list(data) && length(data) > 1) {
    spe <- purrr::reduce2(
      data,
      type,
      function(spe, x, y) {
        S4Vectors::metadata(spe)[[category]][[name]][[y]] <- x
        return(spe)
      },
      .init = spe
    )
  } else {
    S4Vectors::metadata(spe)[[category]][[name]][[type]] <- data
  }

  return(spe)
}


# %% [markdown]
# ### get_metadata
# 
# Simple helper to get the specified metadata using category/name/type structure

# %%
get_metadata <- function(spe, category, name = NULL, type = NULL) {
  if (is.null(name)) {
    S4Vectors::metadata(spe)[[category]]
  } else if (is.null(type)) {
    S4Vectors::metadata(spe)[[category]][[name]]
  } else {
    #Works both for lists and tibble
    S4Vectors::metadata(spe)[[category]][[name]][[type]]
  }
}

# %% [markdown]
# ### mutate_metadata

# %%
mutate_metadata <- function(
  spe,
  category,
  var,
  value,
  name = NULL,
  type = NULL
){
  #Wrapper aroung dplyr::mutate to apply
  #mutate to the one metadata and store it after
  #Will coerce metadata as tibble
  suppressPackageStartupMessages(stopifnot(require(c("dplyr", "rlang"))))
  spe %>%
    get_metadata(
      category = category,
      name = name,
      type = type
    ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      !!rlang::sym(var) := value
    ) %>%
    store_metadata(
      spe,
      data = .,
      category = category,
      name = name,
      type = type
    )
}


# %% [markdown]
# ### choose k

# %%
choose_k <- function(spe, cluster_name, k, replace = FALSE){
  #k : integer vector, corresponding to the chosen metacluster k
  # replace : logical, whether to replace (true) or add to the existing k if any
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  if(replace) {
    spe %>%
      store_metadata(
        k,
        category  = "clustering_data",
        name = cluster_name,
        type = "metacluster_k"
      )
  } else {
    spe %>%
      get_metadata(
        category  = "clustering_data",
        name = cluster_name,
        type = "metacluster_k"
      ) %>%
      c(k) %>%
      store_metadata(
        spe,
        data = .,
        category  = "clustering_data",
        name = cluster_name,
        type = "metacluster_k"
      )
  }
}

# %% [markdown]
# ## address

# %% [markdown]
# ### address_to_metacluster

# %%
address_to_metacluster <- function(
  cluster_name,
  k,
  spe
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  #Will skip verification of cluster_name
  # Assuming it is already checked in the calling function
  use_meta <- get_metadata(
    spe,
    category = "clustering_methods",
    name = cluster_name
  ) %>%
    .[[1]] %>%
    .$metaclustering


  if (!use_meta) {#short cirtcuit
    return(cluster_name)
  } else {
    #Retrieve the k values stored in parameters
    k_values <- get_metadata(
      spe,
      category = "clustering_data",
      name = cluster_name,
      type = "metacluster_k"
    )
    #prepare the name vector
    meta_names <- cluster_name %>%
      paste0("_k", k_values)
    #Return the vector
    meta_names
  }
}

# %% [markdown]
# ### address_to_umap

# %%
address_to_umap <- function(
  cluster_name,
  spe
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  #For clarity
  dim_red <- spe %>%
    get_metadata(
      category = "clustering_methods",
      name = cluster_name,
      type = "reducedDim"
    )

  #Address to the corresponding UMAP(s)
  if(is.null(dim_red)){
    # Return a default if NULL (which indicates metadata not found)
    # To avoid failure
    return("UMAP")
  }
  dplyr::case_when(
    grepl("UMAP", dim_red, ignore.case = TRUE) ~ dim_red,
    grepl("harmony", dim_red, ignore.case = TRUE) ~ "UMAP_harmony_corrected",
    grepl("fastmMNN", dim_red, ignore.case = TRUE) ~ "UMAP_mnn_corrected",
    .default = "UMAP"
  )
}

# %% [markdown]
# ## prepare_device

# %%
prepare_device <- function(device_type = c("pdf", "postscript", "svg",
                  "bitmap", "png", "jpeg", "bmp", "tiff")) {
  # Function factory which prepares the device function to use
  # All the functions will have the same file argument and ...
  # Raster graphic functions will all use inches
  # Position matching beyond the file argument is strongly discouraged,
    #Results will vary depending on the device called  
  switch(
    device_type,
    "pdf" = function(file, ...) grDevices::pdf(file, ...),
    "postscript" = function(file, ...) grDevices::postscript(file, ...),
    "svg" = function(file, ...) grDevices::svg(filename = file, ...),
    "bitmap" = function(file, ...) grDevices::bitmap(file, ...),
    "png" = function(file, ...)  grDevices::png(..., filename = file, units = "in"), #nolint
    "jpeg" = function(file, ...)  grDevices::jpeg( ..., filename = file, units = "in"), #nolint
    "bmp" = function(file, ...)  grDevices::bmp(..., filename = file, units = "in"), #nolint
    "tiff" = function(file, ...)  grDevices::tiff(..., filename = file, units = "in"), #nolint
    .default = stop("Error, unrecognized graphics device")
  )
}

# %% [markdown]
# ## prepare_mapper
# 
# Simple function, prepared to do great things, but for now only returns each marker

# %%
prepare_mapper <- function(
  spe,
  variable,
  subset = NULL
) {
  mapper <- switch(
    variable,
    "marker" = rownames(spe),
    "colData" = names(colData(spe)),
    "reducedDim" = SingleCellExperiment::reducedDimNames(spe)
  )

  if (!is.null(subset)) {
    if(is.logical(subset)){
      mapper <- mapper[subset]
    } else {
      mapper <- mapper[mapper %in% subset]
    }
  }
}

# %% [markdown]
# ## Extract_roi
# 
# Function used in the preprocessing steps to extract only the corresponding roi

# %%
extract_roi <- function(df, file_name) {
  stopifnot(#Base pipe not compatible here
    suppressPackageStartupMessages(
      require("dplyr")
    )
  )

  roi <- file_name %>%
    stringr::str_split_i("_", 1) %>%
    stringr::str_split_i("roi", 2) %>%
    paste0("ROI 00", .)

  col_roi <- grep("ROIName", colnames(df), value = TRUE)

  df <- df %>%
    dplyr::filter(
      !!rlang::sym(col_roi) == roi
    )
}

# %% [markdown]
# ## organize_and_store
# 
# Function that takes a list of plot and stores them organized either into one file (organized as a grid) or into separate files
# 
# **ARGUMENTS**
# - `plot_list` : list of objects to plot
# - `plot_name` : character vector containing the name of each plot
# - `plot_dir` : character vector, path to the folder where the plot are to be stored
# - `separate` : logical value, whether to plot in separate files or in one file, organized with grid.arrange
# - `device_type` : character value, indicating the graphic device to use. <br> Supported graphic devices are pdf, postscript, svg, bitmap, png, jpeg, bmp and tiff, and only the version from {grDevices}
# - `file_extension` : character value, indicating the extension (like ".pdf") of the file(s). <br> If left to default value NULL, will attribute an extension based on device_type
# - `base_height` and `base_width` : numeric values, base size of the files **in inches**.
# - `...` additional argument passed to the respective graphics device (pdf, svg, etc.)
# 
# **VALUE**<br>
# Creates one or multiple files with the plot in the plot_dir folder. <br>
# Returns plot_list invisibly.
# 
# **DETAILS**
# - ... : due to how the device is called, position and partial argument matching for ... is strongly discouraged, and results may vary depending on the device called
# - height and width : for separate files, the size of each file will be the base dimensions. For grid-arranged files, the size of the file will be arranged, multiplying the base dimension by the number of plot in each row and column respectively.

# %%
organize_and_store <- function(
  plot_list,
  plot_name = names(plot_list),
  plot_dir = "./",
  separate = FALSE,
  device_type = c("pdf", "postscript", "svg",
                  "bitmap", "png", "jpeg", "bmp", "tiff"),
  file_extension = NULL,
  base_height = 5, #Always in inches
  base_width = 5,
  ...
) {

  #stop device on exit
  if (dev.cur() > 0) {
    on.exit(graphics.off(), add = TRUE)
  }

  #Creates the graphic device function depending on device_type
  device_type <- match.arg(device_type)
  graphic_device <- prepare_device(device_type)
  # Prepares file_extension with default values depending on device_type
  if (is.null(file_extension)) {
    file_extension <- switch(
      device_type,
      "pdf" = ".pdf",
      "postscript"  = ".ps",
      "svg" = ".svg",
      "bitmap" = ,
      "png" = ".png",
      "jpeg" = ".jpg",
      "bmp" = ".bmp",
      "tiff" = ".tif",
      .default = stop("Error, unsupported graphics device")
    )
  }

  if (separate) {
    purrr::walk2(
      plot_list,
      plot_name,
      function(x, y) {
        file_name <- paste0(plot_dir, y, file_extension)
        graphic_device(
          file = file_name,
          height = base_height,
          width = base_width,
          ...
        )
        plot(x)
        graphics.off()
      }
    )
  } else {
    file_name <- paste0(plot_dir, plot_name, file_extension)

    #prepare col_layout for grid.arrange
    col_layout <- case_when(
      length(plot_list) <= 3 ~ length(plot_list),
      length(plot_list) <= 6  ~ 2,
      length(plot_list) <= 11 ~ 3,
      length(plot_list) <= 19 ~ 4,
      .default = 5
    )

    graphic_device(
      file = file_name,
      width = col_layout * base_width,
      height = ceiling(length(plot_list) / col_layout) * base_height,
      ...
    )
    gridExtra::grid.arrange(
      grobs = plot_list,
      ncol = col_layout
    )
    graphics.off()
  }
  invisible(plot_list)
}

# %% [markdown]
# ## translate_metacluster

# %%
translate_metaclusters <- function(
  cluster,
  spe,
  k
) {
  # take a name of a cluster in spe and returns a VECTOR (factor)
  # With the value of the metacluster for each cell
  # I am not verifying the presence of k in the metaclusters,
  # But if any error, should be within dplyr::pull, so usually helpfull error
  suppressPackageStartupMessages(stopifnot(require("dplyr")))

  spe %>%
    verify_colData(cluster) %>%
    verify_metadata(
      "meta_code",
      category = "metaclustering",
      name = cluster
    ) %>%
    get_metadata(
      category = "metaclustering",
      name = cluster,
      type = "meta_code"
    ) %>%
    dplyr::pull(paste0("meta", k)) %>%
    .[SummarizedExperiment::colData(spe)[[cluster]]] %>%
    as.factor()
}

# %% [markdown]
# ### store_translation

# %%
store_translation <- function(spe, cluster, k, name = paste0(cluster, "_", k)) {
  #SImple wrapper aroung translate_metacluster that can store it into spe
  SummarizedExperiment::colData(spe)[[name]] <- 
    translate_metaclusters(
      cluster,
      spe,
      k
    )

  return(spe)
}

# %% [markdown]
# ## Extract_col_data
# 
# Utility that gets the data associated with colData, returning either the assay or the reducedDim. <br>
# These will be returned as a tibble with rownames in "object_id" which will can be used for merging <br>
# additional colData can be added to that tibble using the additional_col argument. <br>
# For reducedDim, if assign_marker_to_use is true, we can change marker_to_use to be the colnames of the reducedDim (PC1, etc.)

# %%
extract_col_data <- function(
  spe,
  colData_name,
  use_reducedDim = FALSE,
  reducedDim_name = NULL,
  assay_name = "counts",
  unknown_type = "unknown",
  marker_to_use = rownames(spe)[rowData(spe)$use_channel],
  assign_marker_to_use = FALSE,
  env = parent.frame(),
  additional_col = NULL
) {
  #Helper function to extract the colData and assay or reducedDim
  #And return these data as a tibble
  
  #Assign_marker_to_use controls wether to modify marker_to_use
  #in the env environment (typically the caller environment)
  # It is usefull to handle the dimension names (PC1, HARMONY1 etc.)
  # of the reducedDim
    # env is lazily evaluated from the function, 
    # so parent.frame will always return the caller environment
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  data <- spe %>%
    SummarizedExperiment::colData() %>%
    tibble::as_tibble(rownames = "object_id") %>%
    dplyr::select(
      dplyr::all_of(
        c(
          "object_id",
          !!colData_name,
          !!additional_col
        )
      )
    )

  if (!use_reducedDim) {
    assay <- spe %>%
      SummarizedExperiment::assay(assay_name) %>%
      t() %>%
      tibble::as_tibble(rownames = "object_id") %>%
      dplyr::select(dplyr::all_of(c("object_id", marker_to_use)))
  } else {
    #No transposition needed for reducedDim
    assay <- spe %>%
      SingleCellExperiment::reducedDim(reducedDim_name) %>%
      tibble::as_tibble(rownames = "object_id")

    marker_to_use <- assay %>%
      dplyr::select(-object_id) %>%
      colnames()

    if (assign_marker_to_use) {
      assign(
        "marker_to_use",
        value = marker_to_use,
        envir = env
      )
    }
  }

  df <- merge.data.frame(data, assay, by.y = "object_id")

  return(df)
}

# %% [markdown]
# ## Filter_spe
# 
# Helper to subset spe objects based on a colData
# Also makes sure that a object_id exists to keep track of which cells where filtered

# %%
filter_spe <- function(spe, filter_by, keep, object_id_col = "object_id") {
  # Helper to subset spe objects based on a colData
  # Also makes sure that a object_id exists to keep track
  # of which cells where filtered
  if (!object_id_col %in% names(colData(spe))) {
    SummarizedExperiment::colData(spe)[[object_id_col]] <- 
      SummarizedExperiment::colData(spe) %>%
      tibble::as_tibble(rownames = object_id_col) %>%
      dplyr::pull(!!sym(object_id_col))
  }

  spe_subset <- rlang::expr(
    SummarizedExperiment::colData(spe)[[!!filter_by]] %in% !!keep
  )

  spe[, eval(spe_subset)]
}

# %% [markdown]
# ## reintegrate_subset
# 
# Function used to reintegrate a clustering that has been done only on a subset of spe. <br>
# Reintegration is being done via matching the object_id_col

# %%
reintegrate_subset <- function(
  spe,
  subset,
  cluster_col,
  subcluster_col,
  new_col = subcluster_col,
  object_id_col = "object_id"
) {
  df <- SummarizedExperiment::colData(subset) %>%
    tibble::as_tibble() %>%
    dplyr::select(
      dplyr::all_of(c(
        !!object_id_col,
        !!subcluster_col
      ))
    )

  SummarizedExperiment::colData(spe)[[new_col]] <-
    SummarizedExperiment::colData(spe) %>%
    tibble::as_tibble() %>%
    dplyr::left_join(df, by = object_id_col) %>%
    dplyr::mutate(
      !!rlang::sym(subcluster_col) := dplyr::coalesce(
        !!sym(subcluster_col),
        !!sym(cluster_col)
      )
    ) %>%
    dplyr::pull(
      !!rlang::sym(subcluster_col)
    )

  return(spe)
}

# %% [markdown]
# ## Annotation

# %% [markdown]
# ### Add_annontation_category
# 
# Essentially a wrapper around dplyr::case_match to add automatic annotation recoded based on previous annotation. <br>
# Will return the list annotation with new categorisation
# 
# **Arguments**
# - `annotation` : list containing named list of cluster category
# - `from` : character string, name of the list inside annotation to recode
# - `to` : character string, name of the new list to add to annotation. <br>
# By default creates unique names if necessary, see `refuse_name`
# - `matching` : character vector. Each value of the from list will be checked if present in this vector.
# - `true_value` and `false_value` : scalar value, usually character string or logical, value that will be recoded depending on the match.
# - `verbose` : logical value, whether to display more message
# - `refuse_name` : logical value,  whether to throw an error if the name `to` already exists instead of making a unique one with make.unique().
# 
# **Value**
# Returns the annotation list with a new element appended
# 
# **Details**
# The type of true_value and false_value must be the same, so that it can be made into a single atomic vector. <br>
# While character strings ("tumor", "immune_cells", "epithelial_cells") or logical values (tumor: TRUE/FALSE) are likely to be the most usefull, this function does technically support any atomic vector of length one as true_value and false_value, including S3 type vector like date, date-time (POSIXct), duration and factor
# 

# %%
add_annotation_category <- function(
  annotation,
  from,
  to = NULL,
  matching,
  true_value = TRUE,
  false_value = FALSE,
  verbose = FALSE,
  refuse_name = FALSE
) {
  #### Checking arguments ##########################################
  if (!from %in% names(annotation)) {
    abort_not_found(
      arg_name = "from",
      valid = "Starting annotation",
      choices = names(annotation),
      not = from
    )
  }

  if (refuse_name && to %in% names(annotation)) {
    msg <- c(
      "Name 'to' already present in 'annotation'",
      "x" = "Name 'to' has to be unique",
      "i" = "You can set refuse_name to FALSE to automatically handle creation of unique name" #nolint
    )
    cli::cli_abort(msg)
  }

  #matching and from
  if (canCoerce(matching, typeof(from))) {
    if (verbose) message("Coercing matching as typeof(from)")
    matching <- as(matching, typeof(from))
  }
  if (typeof(matching) != typeof(from)) {
    msg <- c(
      "Values of 'matching' must be comparable to the value of 'from'",
      "x" = "Type of 'matching' and 'from' arguments are different",
      "i" = "Type of 'matching' : {typeof(matching)}",
      "i" = "Type of 'from' : {typeof(from)}"
    )
    cli::cli_abort(msg)
  }

  #Values arguments
  if (canCoerce(false_value, typeof(true_value))) {
    if (verbose) message("Coercing false_value as typeof(true_value)")
    matching <- as(false_value, typeof(true_value))
  }
  if (typeof(true_value) != typeof(false_value)) {
    msg <- c(
      "x" = "Type of 'true_value' and 'false_value' must be the same",
      "i" = "Type of 'true_value' : {typeof(true_value)}",
      "i" = "Type of 'false_value' : {typeof(false_value)}"
    )
    cli::cli_abort(msg)
  }
  #prepare a name if none was provided
  #since name will be required when adding to colData
  if (is.null(to)) {
    warning(
      "No name provided in 'to' argument, \n
      true_value will be used as column name.\n"
    )
    to <- true_value
  }
  ### Adding the argument
  annotation <- annotation %>%
    .[[from]] %>%
    dplyr::case_match(
      matching ~ true_value,
      .default = false_value
    ) %>%
    list() %>%
    setNames(to) %>%
    append(annotation, values = .)
  #Make sure the names are unique
  names(annotation) <- make.unique(names(annotation))
  annotation
}

# %% [markdown]
# ### add_annotation

# %%
add_annotation <- function(
  spe,
  annotation,
  init_name
) {
  ### Conditions ##################################################
  #Package
  suppressPackageStartupMessages(stopifnot(require("dplyr")))

  #check spe_object
  if (!inherits(spe, "SpatialExperiment")) {
    cli::cli_abort(message = c(
      "Incorrect spe_object argument",
      "x" = "spe must be a SpatialExperiment",
      "i" = "Provided spe is : {class(spe)}"
    ))
  }
  spe |> verify_colData(init_name)

  # Check  annotation :
  #Prepare cluster_number
  cluster_number <- SummarizedExperiment::colData(spe) %>%
    .[[init_name]] %>%
    unique() %>%
    length()
  #prepare must_list that will check the conditions in must_multiple
  annotation_must <- list(
    "be_list" = \(x) is.list(x),
    "have_a_value_per_cluster" = function(x) {
      #We check for each list inside annotation
      purrr::every(x, function(y) length(y) == cluster_number)
    }
  )
  #Verify the conditions
  annotation <- must_multiple(
    arg_name = "annotation",
    must_list = annotation_must,
    arg = annotation
  )
  ### Execute #############################################
  spe <- purrr::reduce2(
    annotation,
    names(annotation),
    function(spe, ann, name) {
      SummarizedExperiment::colData(spe)[[name]] <-
        ann %>%
        .[SummarizedExperiment::colData(spe)[[init_name]]]

      return(spe)
    },
    .init = spe
  )
}

# %% [markdown]
# # Parameters functions

# %% [markdown]
# ## set_parameter
# 
# Accessor : 
# - parameters
# - global_parameters
# - param_names() gives the list of param_name already stored (including global_parameters)
# 
# Setter :
# - user friendly parameters<- and the pipe-friendly set_parameters
# - user friendly global_parameters and pipe-friendly set_global_parameters

# %%
parameters <- function(spe, record_name) {
  spe |>
    get_metadata(category = "parameters", name = record_name)
}

# %%
global_parameters <- function(spe) {
  spe |>
    get_metadata(category = "parameters", name = "global_parameters")
}

# %%
`parameters<-` <- function(x, record_name, value) {
  if(record_name == "global_parameters") {
    cli::cli_abort(c(
      "Error while setting parameters",
      "x" = "'global_parameters' is reserved",
      "i" = "To modify or access global parameters, please use global_parameters functions" #nolint
    ))
  }

  if(length(value) == 1) {
    name <- names(value)
    value <- value[[1]]
    x |>
    store_metadata(
      unname(value),
      category = "parameters",
      name = record_name,
      type = name
    )
  } else {
    x |>
      store_multiple_metadata(
        value,
        category = "parameters",
        name = record_name,
        type = names(value)
      )
  }
}

# %%
set_parameters <- function(x, record_name, value) {
  parameters(x, record_name) <- value

  return(x)
}

# %%
`global_parameters<-` <- function(x, value) {
  if (length(value) == 1) {
    name <- names(value)
    value <- value[[1]]
    x |>
      store_metadata(
        unname(value),
        category = "parameters",
        name = "global_parameters",
        type = name
      )
  } else {
    x |>
      store_multiple_metadata(
        value,
        category = "parameters",
        name = "global_parameters",
        type = names(value)
      )
  }
}

# %%
set_global_parameters <- function(x, value) {
  global_parameters(x) <- value

  return(x)
}

# %%
param_names <- function(spe) {
  spe |>
    get_metadata(
      category = "parameters"
    ) |>
    names()
}

# %%
map_parameters <- function(spe, param_names, param_list) {

  param_list <- purrr::transpose(param_list)

  if (length(param_list) != length(param_names)) {
    cli::cli_abort(c(
      "Error trying to store multiple parameters",
      "x" = "Number of parameter sets does not match number of names",
      "i" = "Number of parameter sets : {length(param_list)}",
      "i" = "Number of names : {length(param_names)}"
    ))
  }

  spe <- purrr::reduce2(
    param_names,
    param_list,
    function(spe, name, param) {
      parameters(spe, name) <- param
    },
    .init = spe
  )
}

# %% [markdown]
# ## Address_parameters
# 
# This function is designed to be called within other function to access parameters in the metadata.
# 
# **ARGUMENTS :**
# - `spe` : SpatialExperiment object
# - `name` : character string name of the metadata where the parameters are/will be stored
# - `parameters` : **named** list (name-value pairs) of the parameters to address
# - `envir` : environment, the environment of calling function, where the parameters must be evaluated and assigned. See details
# 
# **VALUE**
# 
# Return the spe object, with modified metadata, available in metadata(spe)$parameters[[name]]
# 
# **DETAILS**
# 
# The motivation for this function is to provide a tool that both records in the metadata any parameter used for the functions (for reproductibility and for transparency) AND provide a convenient short-hand to avoid long function calls, by reusing parameters that can be reused.
# 
# Here is the general idea :
# - If an argument is provided in the function call, it will be used in priority and recorded in the metadata
# - If an argument is missing from the function call but is present in the metadata, the value from the metadata will be used
# - We attempt to use specific metadata first if available, then global parameter
# - If an argument is missing from the function call and the metadata, it will use the default value and that default value will be recorded in the metadata
# 
# To do this, this function will look into the environment provided by the envir argument to evaluate whether the parameters are missing from the call are what the default value are. <br>
# 
# Because we try to know if the argument are missing from the call, it is primordial to call this function before any evaluation is made of the parameters. In particular, trying to call match.arg() before this function will make this function not attempt to look within the metadata, and only the value return by match.arg() will be used.
# 
# The defaut value of address_parameters tries to get the correct environment, but may not be applicable to all usage. It is therefore recommanded to explicitly provide an environment like this : <br>
# <code>myfunction <- function(spe, arg1, arg2) {<br>
#   env <- environment() <br>
#   address_parameters(spe, name "my_function", parameters = list(arg1 = arg1, arg2 = arg2), envir = env) <br>
# }</code>
# 
# 

# %%
address_parameters <- function(
  spe,
  param_name,
  parameters,
  envir = parent.frame()
) {
  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  suppressPackageStartupMessages(stopifnot(require("rlang")))
  # evaluate before passing anything to pmap()
  envir <- eval(envir)
  #prepare dataframe for pmap
  df <- prepare_missing(param = parameters, envir = envir)

  #initialize parameter if absent
  if (!"parameters" %in% names(S4Vectors::metadata(spe))) {
    S4Vectors::metadata(spe)[[param_name]] <- list()
  } else if (!param_name %in% names(S4Vectors::metadata(spe)$parameters)) {
    S4Vectors::metadata(spe)$parameters[[param_name]] <- list()
  }
  
  param <- df %>% 
    purrr::pmap(
      prepare_parameters,
      record_name = param_name,
      envir = envir,
      spe = spe
    ) %>%
    setNames(df$parameters)

  #Then we can record these parameter in the metadata
  parameters(spe, param_name) <- param
  # And assign it to the envir
  purrr::walk2(
    param,
    names(param),
    function(param, param_type) {
      assign(param_type, param, envir = envir)
    }
  )

  return(spe)
}

# %% [markdown]
# ### prepare_missing
# 
# Helper that takes the parameters and evaluates in envir wether they are missing from the call and what is the default value (from formals())

# %%
prepare_missing <- function(param, envir) {
  #Prepare a list-column containing the provided values
  # To account for the fact that different types will exist in the parameters
  # First eval the parameters symbol in the provided environment
  parameters <- param %>% 
    purrr::map(\(x) eval(x, envir = envir)) %>%
    tibble::as_tibble_col(column_name = "provided") %>%
    dplyr::mutate(parameters = names(param), .before = 1)

  missing <- purrr::map_lgl(
    names(param),
    function(arg) {
      eval(rlang::expr(missing(!!arg)), envir = envir)
    }
  ) %>%
    setNames(names(param)) %>%
    as.list() %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      cols = colnames(.),
      names_to = "parameters",
      values_to = "missing"
    )

  formals <- eval(rlang::expr(formals()), envir = envir)

  formal <- formals %>%
    as.list() %>%
    tibble::as_tibble_col(column_name = "formals") %>%
    dplyr::mutate(parameters = names(formals), .before = 1)

  parameters %>%
    dplyr::left_join(missing, by = "parameters") %>%
    dplyr::left_join(formal, by = "parameters")
}

# %% [markdown]
# ### prepare_parameters
# 
# called within address to parameters via pmap, will address one parameter to either the provided value, the recorded value or the default value depending on the situation, returning one value.
# 
# Spe will be found in the upper layers of pmap (in the execution environment of address_parameters). Not passed explicitly to avoid including it in the df for pmap (risk of copy on modify)

# %%
prepare_parameters <- function(
  parameters,
  provided,
  missing,
  formals,
  record_name,
  envir,
  spe
) {
  if (!missing) {
    #if provided, we use the input
    return(provided)
  } else if (require_metadata(
    spe, arg = parameters, category = "parameters", name = record_name
  )) {
    #If not provided but available in the metadata, we use the metadata
    get_metadata(
      spe,
      category = "parameters",
      name = record_name,
      type = parameters
    )
  } else if (require_metadata(
    spe, arg = parameters, category = "parameters", name = "global_parameters"
  )) {# else we look into the global parameters
    get_metadata(
      spe,
      category = "parameters",
      name = "global_parameters",
      type = parameters
    )
  } else {
    # if missing from call and the metadata, we use the default value
    eval(formals, envir = envir)
  }
}

# %% [markdown]
# ## get_called_param
# 
# Helper function designed to prepare the called arguments to pass to adress parameters. <br>
# By default, spe and param_name are excluded, which is likely the most common behavior. <br>
# 
# Due to the functionment of parameters, arguments passed by ... will always be excluded.

# %%
get_called_param <- function(env, excluded = c("spe", "param_name")) {
  # returns a list in the form of list(arg1 = arg1)
  # containing all the arguments in the formals except for 
  # ... (ellipsis) and the excluded arguments
  # designed to pass the list to address_parameter
  formals <- eval(expr(formals()), envir = env)

  formals <- formals[!names(formals) %in% c(excluded, "...")]
  
  names <- names(formals)
  purrr::map(
    names,
    \(x) expr(!!sym(x))
  ) |> setNames(names)
}

# %% [markdown]
# ## Default_parameter
# 
# Used as default parameter for parameters that don't have a default value. <br>
# It will try to find it in the metadata of spe, but if absent will return an error telling the user that the parameter has to be specified (in the call or in the metadata)

# %%
default_parameter <- function(
  spe,
  param_name,
  param
) {
  if (
    ! param %in% names(parameters(spe, param_name)) &&
      ! param %in% names(parameters(spe, "global_parameters"))
  ) {
    cli::cli_abort(
      c(
        "Incorrect {param} argument",
        "x" = "No {param} argument provided without default value"
      )
    )
  } else if (
    param %in% names(parameters(spe, param_name))
  ) {
    return(parameters(spe, param_name)[[param]])
  } else if(
    param %in% names(parameters(spe, "global_parameters"))
  ) {
    return(parameters(spe, "global_parameters")[[param]])
  }
}

# %% [markdown]
# # Condition functions

# %%
abort_bad_argment <- function(arg_name, must, not = NULL) {
  #Type error handler
  msg <- c("x" = "`{arg_name}` must {must}")
  if (!is.null(not)) {
    not <- typeof(not)
    msg <- c(msg, "i" = "Provided `{arg_name}` is of type {not}")
  }
  cli::cli_abort("error_bad_argument",
            message = msg,
            arg = arg_name,
            must = must,
            not = not)
}

# %%
abort_not_found <- function(arg_name, valid, choices, not = NULL) {
  #Not found error handler
  if(!is.null(not)){
    msg <- c(
      "`{arg_name}` not found in `{valid}`",
      "x" = "`{arg_name}` provided : {not}",
      "i" = "Valid `{valid}` are : {choices}"
    )
  } else {
    msg <- c(
      "`{arg_name}` not found in `{valid}`",
      "i" = "Valid `{valid}` are : {choices}"
    )
  }
  

  cli::cli_abort("error_not_found",
            message = msg,
            arg = arg_name,
            valid = valid,
            choices = choices,
            not = not)
}

# %%
warn_drop <- function(arg_name, valid, choices, not = NULL) {
  #Condition handler to warn that one invalid argument will be dropped
  #Mostly used within other condition handler
  dropped <- not[!not %in% choices]
  msg <- c("!" = "Some invalid `{valid}`  found in `{arg_name}`",
           "!" = "`{arg_name}` provided {not}",
           "i" = "Only valid`{valid}` will be used",
           "i" = "`{valid}` dropped : {dropped}",
           "i" = "Valid `{valid}` are {choices}")

  cli::cli_warn("warn_drop",
           message = msg,
           arg = arg_name,
           valid = valid,
           choices = choices,
           not = not)
}

# %%
match_arg_multiple <- function(arg, valid, choices, arg_name) {
  # Equivalent to match.arg that accept multiple valid values
  if (all(!arg %in% choices)) {
    abort_not_found(arg_name, valid, choices, not = arg)
  } else if (any(!arg %in% choices)) {
    warn_drop(arg_name, valid, choices, not = arg)
  }
  # Return the matched arguments (i.e., the ones that are in choices)
  return(arg[arg %in% choices])
}

# %%
must_multiple <- function(arg_name, must_list, arg) {
  #Condition handler, needs to be given a predicate list (must_list)
  #Will return an error with only the condition that are not met
  cond_pass <- purrr::map_lgl(must_list, ~.x(arg))

  if (all(cond_pass)) {
    return(arg)
  } else {
    must <- names(cond_pass)[!cond_pass & !is.na(cond_pass)]
    abort_bad_argment(arg_name, must = must, not = arg)
  }
}

# %% [markdown]
# ## Argument checker

# %% [markdown]
# ### validate_perform_cluster

# %%
validate_perform_cluster <- function(
  spe,
  assay_name,
  dim_red_for_clusters,
  cluster_method = c("flowsom", "louvain", "blusparam"),
  weight_type_parameter = c("jaccard", "rank"),
  BLUSPARAM = list(),
  BPPARAM = BiocParallel::SerialParam()
) {
  if (dim_red_for_clusters == "none") {
    verify_assay(spe, assay_name)
  } else {
    verify_reducedDim(spe, dim_red_for_clusters)
  }
  cluster_method <- tolower(cluster_method)
  assign(
    "cluster_method",
    match.arg(cluster_method),
    envir = parent.frame()
  )
  assign(
    "weight_type_parameter",
    match.arg(weight_type_parameter),
    envir = parent.frame()
  )
  if (cluster_method == "BLUSPARAM") {
    stopifnot(inherits(BLUSPARAM, "BlusterParam"))
  }
  stopifnot(inherits(BPPARAM, "BiocParallelParam"))
}

# %% [markdown]
# ### validate_build_metacluster

# %%
validate_build_metacluster <- function(
  spe,
  cluster_name,
  assay_name = NULL,
  dim_red = NULL,
  max_k,
  clustering_metadata,
  quiet,
  BPPARAM
) {
  if (!quiet) {
    message("Checking arguments. \n")
  }
  #Check arguments :
  #spe
  stopifnot(inherits(spe, "SingleCellExperiment"))
  spe |>
    verify_metadata(
      arg = cluster_name,
      category = clustering_metadata
    )
  #assay (only used if dim_red is NULL)
  if (is.null(dim_red)) {
    verify_assay(spe, assay_name)
  } else {
    verify_reducedDim(spe, dim_red)
  }

  stopifnot(canCoerce(max_k, "integer"))
  max_k <- as(max_k, "integer")
  if (max_k <= 2 ||
        !is.finite(max_k) ||
        !length(max_k) == 1) {
    msg <- c(
      "Incorrect max_k argument",
      "x" = "max_k argument must be a finite integer scalar strictly superior to one", #nolint
      "x" = "max_k provided : {max_k}",
      "i" = "max_k correspond to the maximum number of cluster to test,",
      "i" =  "with value between 2 and maxK been tested"
    )
    cli::cli_abort(msg)
  }

  stopifnot(inherits(BPPARAM, "BiocParallelParam"))

  invisible(TRUE)
}

# %% [markdown]
# ### check_prepare_col

# %%
check_prepare_col <- function(
  spe,
  remove_clustering,
  only_keep,
  metaclusters
) {

  suppressPackageStartupMessages(stopifnot(require("dplyr")))
  stopifnot(inherits(spe, "SingleCellExperiment"))

  if(!"clustering_methods" %in% names(metadata(spe))){
    cli::cli_abort(c(
      "Error, no clustering_methods found in metadata of spe",
      "x" = "spe must have clustering_methods recording the clustering done",
      "i" = "this metadata is usually created during multi_cluster function",
      "i" = "available metadata : {names(metadata(spe))}"
    ))
  }
  
  clustering_methods <- metadata(spe)$clustering_method

  if (!is.null(remove_clustering)) {
    if (!all(remove_clustering %in% clustering_methods$name)) {
      abort_not_found(
        "remove_clustering",
        "Clusters to remove",
        choices = clustering_methods$name,
        not = remove_clustering
      )
    }
  }

  if (!is.null(only_keep)) {
    if(!all(only_keep %in% clustering_methods$name)){
      abort_not_found(
        "only_keep",
        "Clusters to keep",
        choices = clustering_methods$name,
        not = only_keep
      )
    }
  }

  if (!is.null(metaclusters)) {
    #Number of clustering for which the metaclustering has been done
    n_clustering <- sum(clustering_methods$metaclustering)

    
    
    if (length(metaclusters) == 1){
      metaclusters <- must_multiple(
        "metaclusters",
        list("be_numerical" = \(x) is.numeric(x)),
        metaclusters
      )
      #Vectorize for each clustering
      metaclusters <- rep_len(metaclusters, n_clustering)
    } else if (length(metaclusters) > 1) {
      metaclusters <- must_multiple(
        "metaclusters",
        list(
          "have_correct_length" = \(x) length(x) == n_clustering,
          "be_list_or_numeric" = \(x) is.list(x) || is.numeric(x)
        ),
        metaclusters
      )
    } else {
      cli::cli_abort(c(
        "Error incorrect argument length for 'metaclusters'",
        "x" = "Metaclusters must have length 1 or equal to the number of clustering", #nolint
        "i" = "Length of provided metaclusters : {length(metaclusters)}",
        "i" = "Number of clustering : {n_clustering}"
      ))
    }

    if (!all(purrr::map_lgl(metaclusters, \(x) all(x >= 2)))) {
      cli::cli_abort(
        message = c(
          "Incorrect values in metaclusters.",
          "x" = "Values in metaclusters must be higher than 2",
          "i" = "Provided metaclusters values : {metaclusters}"
        )
      )
    }
  }  
  metaclusters <<- metaclusters

  invisible(spe)
}

# %% [markdown]
# ### validate_compare_indication

# %%
validate_compare_indication <- function(
  spe,
  col_to_compare,
  indication_col,
  patient_id_col,
  paired,
  sample_id_col,
  clever_test,
  force_test = c("wilcox", "t.test", "welch"),
  count_type = c("count", "proportion", "density"),
  tissue_area_name,
  quiet,
  adjust = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
  signif_threshold = 0.05,
  skip_normality
) {
  if(!inherits(spe, "SummarizedExperiment")) {
    cli::cli_abort(c(
      "Error, spe object not a SingleCellExperiment",
      "i" = "Class of provided spe object : {class(spe)}"
    ))
  }
  purrr::walk(
    c(paired, clever_test, quiet, skip_normality),
    \(x) stopifnot(is.logical(x))
  )
  if (skip_normality && clever_test){
    cli::cli_abort(c(
      "Unable to use clever_test without normality assessment",
      "x" = "Clever_test is not compatible with skip_normality"
    ))
  }

  spe |>
    verify_metadata(var = tissue_area_name) |>
    verify_multiple_colData(
      c(col_to_compare, indication_col, patient_id_col, sample_id_col)
    )

  if (!clever_test) {
    force_test <<- match.arg(force_test)
  }
  count_type <<- match.arg(count_type)
  adjust <<- match.arg(adjust)

  if (count_type == "density") {
    verify_metadata(spe, tissue_area_name, "tissue_area_name")

    tissue_area <- metadata(spe)[[tissue_area_name]]
    if(paired){
      tissue_id <- unique(colData(spe)[[sample_id_col]])
    } else {
      tissue_id <- unique(colData(spe)[[patient_id_col]])
    }

    if (any(!names(tissue_area) %in% tissue_id)) {
      cli::cli_abort(
        c(
          "Error in tissue_are_name",
          "x" = "For density counts, please provide a named vector with tissue area of each sample", #nolint
          "i" = "Provided sample_id : {tissue_id}",
          "i" = "Provided names of tissue_area : {names(tissue_area)}"
        )
      )
    }
  }
}

# %% [markdown]
# ### validate_call_test

# %%
validate_call_test <- function(
  formula,
  x,
  y,
  paired,
  test = c("wilcox", "t_test", "welch")) {

  if (paired) {
    if (is.null(x) || is.null(y)) {
      cli::cli_abort(c(
        "Unsupported statistical test",
        "x" = "x and y data are necessary for paired test",
        "i" = "Provided x : {x}",
        "i" = "Provided y : {y}"
      ))
    }
  }
  if(is.null(x) && is.null(formula)){
    cli::cli_abort(c(
        "Unsupported statistical test",
        "x" = "No data recognized in x or in formula",
        "i" = "Provided x = {x}",
        "i" = "Provided formula = {formula}"
    ))
  }
  #Handle test depending on the input, keeping the lowest common test
  if(all(test == "t_test")) {
    test <- "t_test"
  } else if (all(test %in% c("welch", "t_test"))){
    test <- "welch"
  } else {
    test <- "wilcox"
  }

  return(test)
}

# %% [markdown]
# ### validate_train_model

# %%
validate_train_model <- function(
  spe,
  colData_name,
  model_name = colData_name,
  assay_name = "counts",
  marker_to_use =  rownames(spe)[SummarizedExperiment::rowData(spe)$use_channel], #nolint
  use_reducedDim = FALSE,
  reducedDim_name = "HARMONY",
  ncores
) {
  if (!inherits(spe, "SummarizedExperiment")) {
    cli::cli_abort(c(
      "Error, incorrect spe input",
      "x" = "spe must be a SpatialExperiment, SingleCellExperiment or SummarizedExperiment", #nolint
      "i" = "Provided spe is a : {class(spe)}"
    ))
  }

  if (!is.numeric(ncores) ||
        ncores > parallel::detectCores()) {
    cli::cli_abort(
      c(
        "Incorrect ncores argument",
        "x" = "ncores must be a number lower than the number of available CPU cores",
        "i" = "Number of available CPU cores : {parallel::detectCores()}"
       )
    )
  }

  verify_colData(spe, colData_name)
  if (!use_reducedDim) verify_markers(spe, marker_to_use)
  if (use_reducedDim) verify_reducedDim(spe, reducedDim_name)
}

# %% [markdown]
# ### validate_classify_cell

# %%
validate_classify_cell <- function(
  spe,
  model_name,
  colData_name,
  model = NULL,
  marker_to_use =  rownames(spe)[SummarizedExperiment::rowData(spe)$use_channel],
  use_reducedDim = FALSE,
  reducedDim_name = "HARMONY",
  unknown_type = c("unknown")
) {

  verify_colData(spe, colData_name)
  if (!use_reducedDim) verify_markers(spe, marker_to_use)
  if (use_reducedDim) verify_reducedDim(spe, reducedDim_name)

  if (is.null(model)) {
    model <- spe |>
      verify_metadata(
        model_name,
        category = "predictive_model"
      ) |>
      get_metadata(
        category = "predictive_model",
        name = model_name,
        type = "model"
      )

    assign("model", model, envir = parent.frame())
  }
}

# %% [markdown]
# ### validate_perform_neighborhood

# %%
validate_perform_neighborhood <- function(
  spe,
  colPairName = "knn",
  neighbor_name = "neigborhood",
  aggregate_by = c("metadata", "expression", "reducedDim"),
  count_by = "celltype",
  assay_name = "counts",
  subset_row = NULL,
  reduced_dim_name = NULL,
  statistic = c("mean", "median", "sd", "var"),
  k_value = c(6),
  kmean_name = c("neighborhood_clusters")
) {
  stopifnot(inherits(spe, "SummariedExperiment"))
  aggregate_by <- match.arg(aggregate_by)

  spe |>
    verify_colPair(colPairName) |>
    verify_colData(neighbor_name)

  dplyr::case_match(
    aggregate_by,
    "metadata" ~ verify_colData(spe, count_by),
    "expression" ~ verify_assay(spe, assay_name),
    "reducedDim" ~ verify_reducedDim(spe, reduced_dim_name)
  )

  assign(aggregate_by, envir = parent.frame())
  statistic |> match.arg() |> assign(envir = parent.frame())
  invisible(spe)
}

# %% [markdown]
# ## Verify series
# 
# Verify_* functions are argument checker that assert if one argument is within the colData, rowData, or metadata of spe respectively. <br>
# If absent, they send an error via cli_abort, if present, they return spe invisibly to make them compatible with pipes.

# %%
verify_colData <- function(spe, arg) {
  arg_name <- rlang::enexpr(arg) |>
    as.character()

  if (!arg %in% names(SummarizedExperiment::colData(spe))) {
    cli::cli_abort(
      c(
        "Error, incorrect {arg_name} colData",
        "x" = "{arg_name} not found in spe",
        "x" = "Provided {arg_name} : {arg}",
        "i" = "Available colData : {names(SummarizedExperiment::colData(spe))}"
      )
    )
  }

  invisible(spe)
}

# %%
verify_multiple_colData <- function(spe, args) {
  #verifying multiple colData at once is a common pattern
  # It is recommended to call this function like this :
    # verify_multiple_colData(spe, args = c(col1, col2))
  # rather than :
    # arg_list <- list((col1, col2))
    # verify_multiple_colData(spe, args = arg_list)
  # The first call will allow this function to return
  # the names of the missing arguments
  arg_name <- rlang::enexpr(args) %>%
    as.character() %>%
    .[-1]
  
  missing <- !args %in% names(colData(spe))
  if (any(missing)) {
    cli::cli_abort(
      c(
        "Error, incorrect {arg_name[missing]} colData",
        "x" = "{arg_name[missing]} not found in spe",
        "x" = "Provided {arg_name[missing]} : {args[missing]}",
        "i" = "Available colData : {names(SummarizedExperiment::colData(spe))}"
      )
    )
  }
  invisible(spe)
}

# %% [markdown]
# ### verify_metadata
# 
# The most convoluted of the verify series, to look for presence of all the specified category/name/type specified and the presence of arg in the first non specified slot. <br>
# For example if category is NULL, but not name and type, will check for the presence of arg in the category slot, and then for the presence of name and type in their respective slot within that "arg" category. <br>

# %%
verify_metadata <- function(
  spe,
  arg,
  category = NULL,
  name = NULL,
  type = NULL
) {
  arg_name <- rlang::enexpr(arg) |>
    as.character()

  categories <- names(S4Vectors::metadata(spe))
  if (!is.null(category)) {
    meta_names <- names(S4Vectors::metadata(spe)[[category]])
    if (!is.null(name)) {
      meta_types <- names(S4Vectors::metadata(spe)[[category]][[name]])
    } else {
      meta_types <- names(S4Vectors::metadata(spe)[[category]][[arg]])
    }
  } else {
    meta_names <- names(S4Vectors::metadata(spe)[[arg]])
    if (!is.null(name)) {
      meta_types <- names(S4Vectors::metadata(spe)[[arg]][[name]])
    }
  }

  if (is.null(category)) {
    if (!arg %in% categories) {
      cli::cli_abort(c(
        "Error, incorrect {arg_name} metadata",
        "x" = "{arg} not found in metadata(spe)",
        "i" = "Provided {arg_name} : {arg}",
        "i" = "Available metadata : {categories}"
      ))
    } else if (!is.null(name) &&
                 !name %in% meta_names) {
      cli::cli_abort(c(
        "Error, incorrect {arg_name} metadata",
        "x" = "{name} not found in metadata(spe)[[{arg}]]",
        "i" = "Provided {arg_name} : {arg}",
        "i" = "Available metadata : {meta_names}"
      ))
    } else if (!is.null(name) && !is.null(type) &&
                 !type %in% meta_types){ 
      cli::cli_abort(c(
        "Error, incorrect {arg_name} metadata",
        "x" = "{type} not found in metadata(spe)[[{arg}]][[{name}]]",
        "i" = "Provided {arg_name} : {arg}",
        "i" = "Available metadata : {meta_types}"
      ))
    }
  } else {
    if (!category %in% categories) {
      cli::cli_abort(c(
        "Error, incorrect {arg_name} metadata",
        "x" = "{category} not found in metadata(spe)",
        "i" = "Provided category : {category}",
        "i" = "Available metadata : {categories}"
      ))
    }
    if (is.null(name)) {
      if (!arg %in% meta_names) {
        cli::cli_abort(c(
          "Error, incorrect {arg_name} metadata",
          "x" = "{arg} not found in metadata(spe)[[{category}]]",
          "i" = "Provided {arg_name} : {arg}",
          "i" = "Available metadata : {meta_names}"
        ))
      } else if (!is.null(type) &&
                   !type %in% meta_types) {
        cli::cli_abort(c(
          "Error, incorrect {arg_name} metadata",
          "x" = "{type} not found in metadata(spe)[[{category}]][[{arg}]]",
          "i" = "Provided {arg_name} : {arg}",
          "i" = "Available metadata : {meta_types}"
        ))
      }
    } else {
      if (!name %in% meta_names) {
        cli::cli_abort(c(
          "Error, incorrect {arg_name} metadata",
          "x" = "{name} not found in metadata(spe)[[{category}]]",
          "i" = "Provided {arg_name} : {arg}",
          "i" = "Available metadata : {meta_names}"
        ))
      } else if (!arg %in% meta_types) {
        cli::cli_abort(c(
          "Error, incorrect {arg_name} metadata",
          "x" = "{arg} not found in metadata(spe)[[{category}]][[{name}]]",
          "i" = "Provided {arg_name} : {arg}",
          "i" = "Available metadata : {meta_types}"
        ))
      }
    }
  }
  invisible(spe)
}

# %%
require_metadata <- function(
  spe,
  arg,
  category = NULL,
  name = NULL,
  type = NULL
) {
  # require_metadata works simillarily from verify_metadata,
  # Excepts it does not trigger a error, instead
  # returs FALSE if missing, and TRUE if present
  arg_name <- rlang::enexpr(arg) |>
    as.character()

  categories <- names(S4Vectors::metadata(spe))
  if (!is.null(category)) {
    meta_names <- names(S4Vectors::metadata(spe)[[category]])
    if (!is.null(name)) {
      meta_types <- names(S4Vectors::metadata(spe)[[category]][[name]])
    } else {
      meta_types <- names(S4Vectors::metadata(spe)[[category]][[arg]])
    }
  } else {
    meta_names <- names(S4Vectors::metadata(spe)[[arg]])
    if (!is.null(name)) {
      meta_types <- names(S4Vectors::metadata(spe)[[arg]][[name]])
    }
  }

  if (is.null(category)) {
    if (!arg %in% categories) {
      return(FALSE)
    } else if (!is.null(name) &&
                 !name %in% meta_names) {
      return(FALSE)
    } else if (!is.null(name) && !is.null(type) &&
                 !type %in% meta_types){ 
      return(FALSE)
    }
  } else {
    if (!category %in% categories) {
      return(FALSE)
    }
    if (is.null(name)) {
      if (!arg %in% meta_names) {
        return(FALSE)
      } else if (!is.null(type) &&
                   !type %in% meta_types) {
        return(FALSE)
      }
    } else {
      if (!name %in% meta_names) {
        return(FALSE)
      } else if (!arg %in% meta_types) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

# %%
verify_rowData <- function(spe, arg) {
  arg_name <- rlang::enexpr(arg) |>
    as.character()

  if (any(!arg %in% names(SummarizedExperiment::rowData(spe)))) {
    cli::cli_abort(c(
      "Error, incorrect {arg_name} rowData",
      "x" = "{arg_name} not found in spe",
      "i" = "Provided {arg_name} : {arg}",
      "i" = "Available rowData : {names(SummarizedExperiment::rowData(spe))}"
    ))
  }

  invisible(spe)
}

# %%
verify_colPair <- function(spe, arg) {
  arg_name <- rlang::enexpr(arg) |>
    as.character()

  if (any(!arg %in% SingleCellExperiment::colPairNames(spe))) {
    cli::cli_abort(c(
      "Error, incorrect {arg_name} colPair",
      "x" = "{arg_name} not found in spe",
      "i" = "Provided {arg_name} : {arg}",
      "i" = "Available colPair : {SingleCellExperiment::colPairNames(spe)}"
    ))
  }

  invisible(spe)
}

# %%
verify_reducedDim <- function(spe, arg) {
  arg_name <- rlang::enexpr(arg) |>
    as.character()

  if (any(!arg %in% SingleCellExperiment::reducedDimNames(spe))) {
    cli::cli_abort(c(
      "Error, incorrect {arg_name} reducedDim",
      "x" = "{arg_name} not found in spe",
      "i" = "Provided {arg_name} : {arg}",
      "i" = "Available reducedDim : {SingleCellExperiment::reducedDimNames(spe)}"
    ))
  }

  invisible(spe)
}

# %%
verify_markers <- function(spe, arg) {
  arg_name <- rlang::enexpr(arg) |>
    as.character()

  if(is.logical(arg)) arg <- rownames(spe)[arg]

  if (any(!arg %in% rownames(spe))) {
    cli::cli_abort(c(
      "Error, incorrect {arg_name} colPair",
      "x" = "{arg_name} not found in spe",
      "i" = "Provided {arg_name} : {arg}",
      "i" = "Available markers : {rownames(spe)}"
    ))
  }

  invisible(spe)
}

# %%
verify_assay <- function(spe, arg) {
  arg_name <- rlang::enexpr(arg) |>
    as.character()

  if (!arg %in% SummarizedExperiment::assayNames(spe)) {
    cli::cli_abort(c(
      "Error, incorrect {arg_name} rowData",
      "x" = "{arg_name} not found in spe",
      "i" = "Provided {arg_name} : {arg}",
      "i" = "Available assay : {SummarizedExperiment::assayNames(spe)}"
    ))
  }

  invisible(spe)
}

# %%
verify_spatialCoords <- function(spe, arg) {
  arg_name <- rlang::enexpr(arg) |>
    as.character()

  if (!arg %in% SpatialExperiment::spatialCoordsNames(spe)) {
    cli::cli_abort(c(
      "Error, incorrect {arg_name} spatialCoordsNames",
      "x" = "{arg_name} not found in spe",
      "i" = "Provided {arg_name} : {arg}",
      "i" = "Available colData : {SpatialExperiment::spatialCoordsNames(spe)}"
    ))
  }

  invisible(spe)
}


