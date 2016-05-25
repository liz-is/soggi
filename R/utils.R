check_params <- function(param, values){
  if (!(param %in% values)){
    stop(paste0(param, " not one of: ", paste(values, collapse = ", ")))
  }
}

check_param_type <- function(param, type){
  if(!is(param, type)){
    stop(paste0(param, " should be class ", type, ", but is ", class(param)))
  }
}
