

#' This function combines the files in the given directory and returns as a data table or saves it to the out_put_file_path.
#' Option 1 : Set header = T and no col_names. In this case the function will read the file and use the header of the first 
#' file in the directory. 
#' Option 2 : Set header = F and specify col_names manually. In this case the function will read the file and will use the given array as column names
#' The assumption is that it does not have a header (default)
#' 
#' @param dir_path Absolute path of the directory containing the files to be combined. 
#' @param header Boolean value of if the header needs to be considered
#' @param col_names Character vector to replace the header with a custom column name.
#' @param Outputl File path - absolute path. By default its null and can be omitted. Produces a tsv file.
dir_rbind <- function(dir_path, header = F, col_names = NULL, out_file_path = NULL) {
  assert_that(isDirectory(dir_path))
  filename_list <- list.files(dir_path, full.names = T)
  
  assert_that(length(filename_list) >= 1)
  cat("Reading Init",filename_list[[1]])
  combine_data <- fread(filename_list[[1]], sep="\t", header = header, na.strings="NA",
                        stringsAsFactors = FALSE)
  if (is.null(header) & is.null(col_names)) {
    warning("Header and Col Names are both NULL")
    col_names <- seq(1:dim(combine_data)[1], by=1)
  }
  
  flog.info(paste0("reading "))
  
  number_of_cpus <- detectCores(all.tests = T, logical = T)
  cl <- makeCluster(number_of_cpus - 1)  
  registerDoParallel(cl)
  
  combine_data <- foreach(l = filename_list, .packages='data.table', .combine = rbind)  %dopar% {
    data <- fread(l, sep="\t", header = header, na.strings="NA",
                  stringsAsFactors = FALSE, verbose =T)
    data
  }
  
  stopCluster(cl)
  

  combine_data <- unique(combine_data)
  flog.info(paste0("reading over",filename_list))
  flog.info(paste0("Out Dataframe : rows ",nrow(combine_data), " and cols : ", ncol(combine_data)))
  
  
  if (is.null(header) & is.null(col_names)) {
    warning("Header and Col Names are both NULL")
    setnames(combine_data, names(combine_data), col_names)
  } else {
    setnames(combine_data, names(combine_data), names(combine_data))
  }
  if (!is.null(out_file_path)) {
    write.table(x = combine_data, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(combine_data)
}

