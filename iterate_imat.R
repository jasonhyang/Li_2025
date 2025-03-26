library(gembox)
library(tools)
library(parallel)
library(reticulate)

load_RData <- function(file_name) {
  load(file_name)
  get(ls()[ls() != "file_name"])
}

struct_make <- function(result_par) {
  if (typeof(result_par) == "S4") {
    return(result_par)
  } else if (typeof(result_par) == "character" && length(result_par) == 1) {
    return(result_par)
  } else {
    return(as.matrix(result_par))
  }
}

individual_imat2mat <- function(rs_data, gem_model, q_low, q_high, exp_dir, sp_io, condition) {
  samp_name <- colnames(rs_data)
  dir.create(paste(exp_dir, "/", samp_name, sep=""))
  
  exprs <- exprs2int(x=rs_data, model=gem_model, q.lo=q_low, q.hi=q_high, na2zero=TRUE)
  imat_result <- imat(expr=exprs, model=gem_model, samp.pars=NULL)
  imat_result[["result.model"]][["grRules"]][is.na(imat_result[["result.model"]][["grRules"]])] <- ""
  
  gem_struct <- list(lapply(imat_result[["result.model"]], struct_make))
  names(gem_struct) <- samp_name
  sp_io$savemat(paste(exp_dir, "/", samp_name, "/", samp_name, ".mat", sep=""), gem_struct)
  
  exprs_mat <- as.matrix(exprs)
  colnames(exprs_mat) <- samp_name
  rownames(exprs_mat) <- gem_model[["genes"]]
  write.table(x=exprs_mat, file=paste(exp_dir, "/", condition, "_", samp_name, " ", "GD.tsv", sep=""), sep="\t",
              col.names=NA)
}

iterate_imat2mat <- function(rs_data, i, gem_model, q_low, q_high, exp_dir, sp_io) {
  samp_name <- colnames(rs_data)[i]
  dir.create(paste(exp_dir, "/", samp_name, sep=""))
  
  rs_mat <- as.matrix(rs_data[, i])
  colnames(rs_mat) <- samp_name
  
  exprs <- exprs2int(x=rs_mat, model=gem_model, q.lo=q_low, q.hi=q_high,
                     na2zero=TRUE)
  imat_result <- imat(expr=exprs, model=gem_model, samp.pars=NULL)
  imat_result[["result.model"]][["grRules"]][is.na(imat_result[["result.model"]][["grRules"]])] <- ""
  
  gem_struct <- list(lapply(imat_result[["result.model"]], struct_make))
  names(gem_struct) <- samp_name
  sp_io$savemat(paste(exp_dir, "/", samp_name, "/", samp_name, ".mat", sep=""), gem_struct)
  
  exprs_mat <- as.matrix(exprs)
  colnames(exprs_mat) <- samp_name
  return(exprs_mat)
}

iterate_imat <- function(data_dir, gem_dir, exp_dir, q_low, q_high, run_samp, num_samp,
                         conda_env, condition) {
  use_condaenv(conda_env)
  sp_io <- import("scipy.io")
  
  if (file_ext(gem_dir) == "mat") {
    gem_model <- read.matlab.model(gem_dir)
  } else if (file_ext(gem_dir) == "RData") {
    gem_model <- load_RData(gem_dir)
  } else if (file_ext(gem_dir) == "xml") {
    gem_model <- read.sbml.model(gem_dir)
  }
  
  rs_data <- as.matrix(read.table(data_dir, header=TRUE, quote="", row.names=1, sep="\t"))
  
  if (length(colnames(rs_data)) == 1) {
    individual_imat2mat(rs_data=rs_data, gem_model=gem_model, q_low=q_low, q_high=q_high,
                        exp_dir=exp_dir, sp_io=sp_io, condition=condition)
    
    if (run_samp == TRUE) {
      source_python("iterate_pysampling.py")
      individual_pysampling(exp_dir=exp_dir, num_samp=as.integer(num_samp), condition=condition)
      print("iMAT + Sampling is completed")
    } else if (run_samp == FALSE) {
      print("iMAT is completed")
    }
    
  } else {
    exprs_list <- mclapply(seq_along(colnames(rs_data)), iterate_imat2mat, rs_data=rs_data,
                           gem_model=gem_model, q_low=q_low, q_high=q_high, exp_dir=exp_dir, sp_io=sp_io,
                           mc.cores=detectCores()-1)
    exprs_summary <- do.call(cbind, exprs_list)
    rownames(exprs_summary) <- gem_model[["genes"]]
    
    if (run_samp == TRUE) {
      order_list <- as.list(colnames(rs_data))
      source_python("iterate_pysampling.py")
      iterate_pysampling(gem_dir=gem_dir, exp_dir=exp_dir, num_samp=as.integer(num_samp), order_list=order_list,
                         condition=condition)
      write.table(x=exprs_summary, file=paste(exp_dir, "/", condition, " ", "GD.tsv", sep=""), sep="\t",
                  col.names=NA)
      print("iMAT + Sampling is completed")
    } else if (run_samp == FALSE) {
      write.table(x=exprs_summary, file=paste(exp_dir, "/", condition, " ", "GD.tsv", sep=""), sep="\t",
                  col.names=NA)
      print("iMAT is completed")
    }
  }
}