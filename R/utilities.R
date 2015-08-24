

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
  flog.info(paste0("Reading Init ",filename_list[[1]]))
  head_data <- fread(filename_list[[1]], sep="\t", header = header, na.strings="NA",
                        stringsAsFactors = FALSE, nrows = 10)
  if (is.null(header) & is.null(col_names)) {
    warning("Header and Col Names are both NULL")
    col_names <- seq(1:dim(head_data)[1], by=1)
  }
  
  flog.info(paste0("reading "))
  
  number_of_cpus <- detectCores(all.tests = T, logical = T)
  cl <- makeCluster(2)  
  registerDoParallel(cl)
  
  combine_data <- foreach(l = filename_list, .packages='data.table', .combine = rbind)  %dopar% {
    data <- fread(l, sep="\t", header = header, na.strings="NA",
                  stringsAsFactors = FALSE, verbose =T)
    data
  }
  
  stopCluster(cl)
  

  combine_data <- unique(combine_data)
  flog.info(paste0("reading over ",filename_list))
  flog.info(paste0("Out Dataframe : rows ",nrow(combine_data), " and cols : ", ncol(combine_data)))
  
  
  if (is.null(header) & is.null(col_names)) {
    warning("Header and Col Names are both NULL")
    setnames(combine_data, names(combine_data), norm_var_names(col_names))
  } else {
    setnames(combine_data, names(combine_data), norm_var_names(names(combine_data)))
  }
  if (!is.null(out_file_path)) {
    assert_that(is.writeable(out_file_path))
    write.table(x = combine_data, file=out_file_path, quote = F, sep = "\t", row.names = F)
  }
  return(combine_data)
}




#'
#'
#'
#'example call : t <- dir_merge("/home/data/tmp", col_names = c("#CHROM","POS","ID","REF","ALT"),sep = ",", out_file_path="/home/data/combine-tmp.csv")
#'
#'

dir_merge <- function(dir_path, col_names, sep = "\t", out_file_path = NULL) {
  assert_that(isDirectory(dir_path))
  filename_list <- list.files(dir_path, full.names = T)
  
  assert_that(length(filename_list) >= 1)
  custom_merge <- function(x, y) {
    flog.debug(paste0("x names : ", names(x)))
    flog.debug(paste0("y names : ", names(y)))
    merge(x, y, by=norm_var_names(col_names), all = T)
  }
  
  
  flog.info(paste0("reading "))
  
  number_of_cpus <- detectCores(all.tests = T, logical = T)
  cl <- makeCluster(2)  
  registerDoParallel(cl)
  
  combine_data <- foreach(l = filename_list, .packages='data.table', .combine = custom_merge)  %do% {
    df <- fread(l, sep = sep, header = T, na.strings="NA",
                  stringsAsFactors = FALSE)
    data.table::setnames(df, names(df), norm_var_names(names(df)))
    suffix_list <- setdiff(names(df), norm_var_names(col_names))
    data.table::setnames(df, suffix_list, norm_var_names(paste0(suffix_list, "_", basename(l))))
    
    df
  }
  
  stopCluster(cl)
  
  
  combine_data <- unique(combine_data)
  flog.info(paste0("reading over",filename_list))
  flog.info(paste0("Out Dataframe : rows ",nrow(combine_data), " and cols : ", ncol(combine_data)))
  
  if (!is.null(out_file_path)) {
    write.table(x = combine_data, file=out_file_path, quote = F, sep = sep, row.names = F)
  }
  return(combine_data)
}






#' This function normalises the input string vector
#' _ and small case output, no spaces, no . etc
#' @param a character vector of names
#' @return a normalised character vector of names
norm_var_names <- function(vars, sep="_") {
  
  assert_that(is.character(vars))
  
  if (sep == ".") sep <- "\\."
  
  vars <- make.names(vars)
  
  # Replace all _ and . and ' ' with the nominated separator.
  
  pat  <- '_|\\.| |,'
  rep  <- sep
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  
  # Remove any resulting initial or trailing underscore or multiples:
  #
  # _2level -> 2level
  
  vars <- stringr::str_replace(vars, "^_+", "")
  vars <- stringr::str_replace(vars, "_+$", "")
  vars <- stringr::str_replace(vars, "__+", "_")
  
  # Convert to lowercase
  
  vars <- tolower(vars)
  
  # Remove repeated separators.
  
  pat  <- paste0(sep, "+")
  rep  <- sep
  vars <- stringr::str_replace_all(vars, pat, rep)
  
  return(vars)
}


















