make_ped <- function(data, trafo_fct, cut, cr, n_cr, multimodal) {
  if (multimodal) {
    make_multimodal_ped(data, trafo_fct, cut, cr, n_cr)
  } else {
    make_ped_(data, trafo_fct, cut, cr, n_cr)
  }
}

make_ped <- function(data, formulas, trafo_fct, cut, cr, n_cr, multimodal) {
  if (multimodal) {
    data_unstruct <- data[2:length(data)]
    data <- data[[1]]
  }
  if (cr) {
    crs <- 1:n_cr
    data_ <- ped_data <- X <- X2 <- P <- sp <- pam <- vector("list", n_cr)
  } else {
    crs <- 1
    data_ <- ped_data <- X <- X2 <- P <- sp <- pam <- vector("list", 1L)
  }
  if (!is.list(formulas)) {
    formulas <- rep(list(formulas), length(crs))
  }
  formulas <- decompose_formula(formulas, colnames(data))
  for (i in 1:length(data_)) {
    data_[[i]] <- data %>% mutate(status = ifelse(status == crs[i], 1, 0))
    ped_data[[i]] <- as_ped(data_[[i]], trafo_fct, cut = cut)
    colnames(ped_data[[i]])[colnames(ped_data[[i]]) %in% c("tend", "ped_status")] <- c("time", "Y")
    pam[[i]] <- gam(formulas[[i]]$structured, data = ped_data[[i]], offset = offset, family = "poisson")
    form2 <- as.formula(paste0(
      "Y~",
      paste(
        colnames(ped_data[[i]][,!(colnames(ped_data[[i]]) %in% c("id", "tstart", "interval", "offset", "Y"))]),
        collapse = "+")))
    pam2 <- gam(form2, data = ped_data[[i]], offset = offset, family = "poisson")
    mm <- model.matrix(pam[[i]])
    mm <- cbind(mm, offset = ped_data[[i]]$offset)
    mm2 <- model.matrix(pam2)
    processed_terms <- colnames(mm)[!colnames(mm) %in% c("(Intercept)", colnames(ped_data[[i]]))]
    unprocessed_terms <- c(colnames(ped_data[[i]])[colnames(ped_data[[i]]) %in% colnames(mm)], "(Intercept)", "offset")
    processed_terms_ <- stringr::str_sub(processed_terms, 1,
                                         stringr::str_locate(processed_terms, stringr::fixed('.'))[, 1] - 1)
    processed_terms_t <- split(processed_terms, f = processed_terms_)
    X[[i]] <- vector("list", length(unique(processed_terms_)) + length(unprocessed_terms))
    P[[i]] <- vector("list", length(unique(processed_terms_)))
    names(X[[i]]) <- c(unique(processed_terms_), unprocessed_terms)
    names(P[[i]]) <- unique(processed_terms_)
    for (j in 1:length(processed_terms_t)) {
      X[[i]][[j]] <- mm[, processed_terms_t[[j]], drop = FALSE]
      P[[i]][[j]] <- pam[[i]]$smooth[[j]]$S[[1]] * pam[[i]]$sp[j]
    }
    k = j
    for (j in 1:length(unprocessed_terms)) {
      X[[i]][[k + j]] <- mm[, unprocessed_terms[[j]], drop = FALSE]
    }
    X2[[i]] <- vector("list", 2)#length(formulas[[i]]$deep_vars))
    for (j in 1:length(X2)) {
      use <- mm2[, colnames(mm2) %in% formulas[[i]]$deep_vars[[j]]]
      X2[[i]][[j]] <- matrix(use, nrow = nrow(use), ncol = ncol(use))
    }
  }
  X <- reshape(X, ped_data)
  X2 <- reshape(X2, ped_data)
  if (multimodal) {
    self$data = append(X, X2, data_unstruct)
    self$P = P
    self$weights = weights
  } else {
    self$data = append(X, X2)
    self$P = P
    self$weights = weights
  }
  self$related_pamm = pam
  self$processing_pam = pam2
}
