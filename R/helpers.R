#' @import stats stringr stringi
decompose_formula <- function(formulas, colnamesX) {
  res <- vector("list", length(formulas))
  for (i in 1:length(formulas)) {
    trms <- attributes(terms(formulas[[i]]))
    terms <- trms$term.labels
    brack_terms <- terms[str_detect(terms, stringr::fixed('('))]
    unbrack_terms <- terms[!str_detect(terms, stringr::fixed('('))]
    s_ <- brack_terms[str_detect(brack_terms, stringr::fixed('s('))]
    te_ <- brack_terms[str_detect(brack_terms, stringr::fixed('te('))]
    structured <- c("1", s_, te_, unbrack_terms)
    sform <- paste(structured, collapse = "+")
    sform <- paste0("Y ~ ", sform)
    sformula <- as.formula(sform)
    u <- terms[!(terms %in% structured)]
    re <- "\\(([^()]+)\\)"
    uu <- vector("list", length(u))
    if (length(uu) > 0) {
      for (j in 1:length(u)) {
        uu[[j]] <- gsub(" ", "", gsub(re, "\\1", str_extract_all(u[[j]], re)[[1]]))
      }
      deep_vars <- sapply(uu, strsplit, ",")
      names(deep_vars) <- stringr::str_sub(u, 1, stringr::str_locate(u, stringr::fixed('('))[, 1] - 1)
      use <- rep(F, length(deep_vars))
      for (j in 1:length(deep_vars)) {
        use[j] <- all(deep_vars[[j]] %in% c("time", colnamesX))
      }
      deep_vars <- deep_vars[use]
    } else {
      deep_vars <- NULL
    }
    res[[i]] <- list(structured = sformula, deep_vars = deep_vars)
  }
  res
}

scale_ <- function(scaler, data) {
  for (i in 1:length(scaler[[1]])) {
    data[,names(scaler$means)[i]] <- 
      (data[,names(scaler$means)[i]] - scaler$means[[i]]) / scaler$sds[[i]]
  }
  data
}

#' @import keras tensorflow dplyr
s_layer <- function(x, P, lambda_factor, trainable = TRUE, pretrained_weights = NULL) {
  P = P * lambda_factor
  if (is.null(pretrained_weights)) {
    x %>% layer_dense(1, use_bias = FALSE, trainable = trainable,
                      kernel_initializer = initializer_he_normal(7L),
                      kernel_regularizer = function(x)
                        k_mean(
                          tf$multiply(
                            tf$multiply(k_transpose(x),
                                        tf$constant(P)), x)))  
  } else {
    x %>% layer_dense(1, use_bias = FALSE, trainable = trainable,
                      kernel_initializer = initializer_constant(pretrained_weights),
                      kernel_regularizer = function(x)
                        k_mean(
                          tf$multiply(
                            tf$multiply(k_transpose(x),
                                        tf$constant(P)), x)))  
  }
  
}

data_unlist <- function(x) {
  res <- list()
  if (length(x) > 1) {
    for (i in 1:length(x)) {
      for (j in 1:length(x[[i]])) {
        res <- append(res, x[[i]][j])
      }
    }
  } else {
    res <- x[[1]]
  }
  res
}

smart_append <- function(structured, deep = NULL, unstructured = NULL) {
  is_null_deep <- is.null(deep)
  is_null_unstructured <- is.null(unstructured)
  if (!is_null_unstructured) unstructured <- list(unstructured)
  if (is_null_deep & is.null(unstructured)) {
    structured
  } else if (is_null_deep & !(is_null_unstructured)) {
    append(structured, unstructured)
  } else if (!(is_null_deep) & is_null_unstructured) {
    append(structured, deep) 
  } else {
    append(append(structured, deep), unstructured)
  }
}
#' @import dplyr
recompute_intlen <- function(data) {
  lagged_times <- as.numeric(dplyr::lag(data$times))
  lagged_times[1] <- 0
  intlen <- data$times - lagged_times
  intlen[intlen < 0] <- intlen[1]
  as.numeric(intlen)
}

#' Create list of structured parameters
#' 
#' This function creates a list of parameters used for the structured part
#' of the deeppamm object.
#' The list is passed towards the $make_model method.
#' 
#' @param trainable_struct a logical indicating whether the structured
#' network part should be trained or not. Only use FALSE if model is pretrained
#' (pretain = TRUE).
#' @param pretrain a logical indicating whether the structured network should be 
#' pretrained using a PAMM.
#' @param lambdas a list of lambda values for the structured subnetwork.
#' Each specified effect can receive a lambda value that results in a 
#' L2 penalty on structured coefficient and a prefactor for splines.
#' Defaults to 0 for linear effects and 1 for splines meaning no regularization
#' for linear effects and the default penalization suggested by mgcv for splines.
#' For linear effects the lambda value is simply the regularization value
#' while for splines a value >1 means a stronger regularization than the
#' default while <1 results in smaller regularization.
#' @return a list that is complaint with the <param> argument of the
#' <deeppamm> method <make_net>.
#' 
#' @export
make_params <- function(trainable_struct = FALSE, pretrain = TRUE, lambdas = NULL) {
  list(trainable_struct = trainable_struct, pretrain = pretrain, lambdas = lambdas)
}

subset_mm <- function(data, frac = 0.8) {
  res <- vector("list", 2L)
  train_id <- sample(1:nrow(data[[1]]), round(frac * nrow(data[[1]])))
  test_id <- setdiff(1:nrow(data[[1]]), train_id)
  res[[1]] <- res[[2]] <- vector("list", length(data))
  names(res[[1]]) <- names(res[[2]]) <- names(data)
  for (j in 1:length(res[[1]])) {
    if (length(dim(data[[j]])) == 2L) {
      res[[1]][[j]] <- data[[j]][train_id, , drop = FALSE]
      res[[2]][[j]] <- data[[j]][test_id, , drop = FALSE]
    } else if (length(dim(data[[j]])) == 3L) {
      res[[1]][[j]] <- data[[j]][train_id, , , drop = FALSE]
      res[[2]][[j]] <- data[[j]][test_id, , , drop = FALSE]
    } else if (length(dim(data[[j]])) == 4L) {
      res[[1]][[j]] <- data[[j]][train_id, , , , drop = FALSE]
      res[[2]][[j]] <- data[[j]][test_id, , , , drop = FALSE]
    } else if (length(dim(data[[j]])) == 5L) {
      res[[1]][[j]] <- data[[j]][train_id, , , , , drop = FALSE]
      res[[2]][[j]] <- data[[j]][test_id, , , , , drop = FALSE]
    }
  }
  names(res) <- c("train", "test")
  res
}

