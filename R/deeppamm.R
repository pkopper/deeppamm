#' DEEPPAMM
#' 
#' R6 object to create a deeppamm instance resembling the DEEPPAMM method 
#' proposed in Kopper et al. (2022), https://arxiv.org/abs/2202.07423
#' The instance can be used to create a DEEPPAMM model, fit it, and predict from it.
#' DEEPPAMM is a semi-structured deep learning method for survival analysis.
#' It can model both, unstructured and structured data.
#' Structured data can be modeled through a GAM subnetwork (shallow) or deep subnetwork.
#' DEEPPAMM is capable of fitting competing risks.
#' The initialization does not build the network yet or majorly transforms data 
#' apart from scaling if applicable.
#' The workflow is as follows: $initialize(...), $make_model(...), $train(...), $preidict
#' After making the model a Keras model will be contained in $model; 
#' it can be worked with in the same fashion as all Keras models.
#' Keras and Tensorflow are used in R which are powered through the reticulate
#' package.
#' @import stats Formula reticulate pammtools mgcv keras tensorflow dplyr stringr stringi purrr R6 checkmate survival
#' @export
deeppamm <- R6::R6Class(
  "deeppamm",
  public = list(
    model = NULL,
    formulas = NULL,
    raw_data = NULL,
    data = NULL,
    X = NULL,
    Y = NULL,
    P = NULL,
    weights = NULL,
    deeppamm_data = NULL,
    trafo_fct = NULL,
    cut = NULL,
    cr = FALSE,
    n_cr = 1,
    multimodal = FALSE,
    lr = NULL,
    status = NULL,
    time = NULL,
    scale = TRUE,
    scaler = NULL,
    scaled_data = NULL,
    related_pamm = NULL,
    pam_coeffs = NULL,
    processing_pam = NULL,
    precision = "float32",
    subnetworknames = NULL,
    latest_test_data = NULL,
    deep_architectures = NULL,
    deep = TRUE,
    built = FALSE,
    ped = FALSE,
    subnets = NULL,
    used_vars = NULL,
    Nout = NULL,
    tabular_terms = NULL,
    partial_domain = NULL,
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param formulas
    #' Either a list of formulas for all competing risks or a single formula 
    #' that is applied to all competing risks. The formula MUST contain
    #' Y as the target (lhs) and consider the time variable (rhs; this is to 
    #' account for the baseline hazard). Next to that rules for formulas from
    #' mgcv::gam apply. Additionally, you can include deep architectures in
    #' formulas that need to be supplied in deep_architectures
    #' @param deep_architectures (`list()`)\cr
    #' A named list of all deep architectures that are used in formulas.
    #' See examples to see the construction convention.
    #' @param data 
    #' Either a data.frame (if there is no unstructured data) or a list that 
    #' contains the survival data set in the first entry and all other modalities
    #' in the following ones. All unstructured modalities must be named so that 
    #' they can be matched in formulas
    #' @param trafo_fct (`formula`)\cr
    #' Transformation formula that creates a PED object (pammtools package).
    #' LHS needs to be changed in accordance with variable names for event time
    #' (time) and survival status (status), RHS can be changed if only a subset
    #' of features of the data is neeeded.
    #' @param cut
    #' A vector of discretization points used for the piecewise exponential 
    #' model. Does not need to be provided, in this case all distinct event 
    #' times are taken (not suggested).
    #' @param lr
    #' A numeric scalar to determine the learning rate of the Adam optimizer
    #' used for the NN.
    #' @param scale
    #' A logical that indicates whether or not the structured data should be 
    #' scaled. The event time (even though it is used as feature) is not scaled.
    #' Data augmentation and scaling for the unstructured data needs to be
    #' applied before supplying them. 
    initialize = function(formulas,
                          deep_architectures = list(),
                          data,
                          trafo_fct = Surv(time, status) ~.,
                          cut = NULL,
                          lr = 0.001,
                          scale = TRUE)  {
      check_formulas(formulas)
      check_architectures(deep_architectures)
      check_data(data)
      check_trafo_fct(trafo_fct)
      if (!is.null(cut)) assert_numeric(cut, min.len = 2L, lower = 0, sorted = T)
      assert_number(lr, lower = 0, upper = 1)
      assert_logical(scale, max.len = 1L)
      Y = gsub(" ", "", as.character(trafo_fct)[2])
      time = c(6, unlist(gregexpr(',', Y)[1]) - 1) 
      status = c(time[2] + 2, unlist(gregexpr(')', Y)[1]) - 1)
      self$time = substr(Y, time[1], time[2])
      self$status = substr(Y, status[1], status[2])
      if (is.list(data) & !(is.data.frame(data))) {
        multimodal = TRUE
        data_ <- data[[1]]
      } else {
        multimodal = FALSE
        data_ <- data
      }
      if (length(unique(data_[[self$status]])[unique(data_[[self$status]]) != 0]) > 1) {
        self$n_cr = sum(unique(data_[[self$status]]) != 0)
        self$cr = TRUE
      } else {
        self$n_cr = 1
        self$cr = FALSE
      }
      self$deep_architectures = deep_architectures
      self$raw_data = data_
      self$trafo_fct = trafo_fct
      self$formulas = formulas
      self$cut = cut
      self$lr = lr
      self$scale = scale
      if (scale) {
        scale_data <- data_[, !(colnames(data_) %in% c(self$status, self$time)), drop = FALSE]
        is_numeric <- sapply(scale_data, is.numeric)
        means <- lapply(scale_data[, is_numeric, drop = FALSE], mean)
        sds <- lapply(scale_data[, is_numeric, drop = FALSE], sd)
        self$scaler <- list(means = means, sds = sds)
        data_ = scale_(self$scaler, data_)
        scaled_data <- self$scaled_data
      } 
      if (multimodal) {
        self$multimodal = TRUE
        d <- data
        d[[1]] <- data_
        self$data = d 
      } else {
        self$data <- data_
      }
    }, 
    #' @description 
    #' Changes TF precision
    #' 
    #' @param dtype
    #' The data type to which the precision in TF is to be changed.
    #' Default change is from 32bit (default) to 64bit.
    change_precision = function(dtype = "float64") {
      self$precision = dtype
    },
    #' @description 
    #' Create the model and reshape the data so that we can train a DEEPPAMM
    #' in the following.
    #' 
    #' @param params 
    #' optional parameter, if supplied needs to be a list with parameters 
    #' for the strucured part of the architecture. The list must contain:
    #' params$trainable_struct, params$pretrain, params$lambdas
    #' referring to whether or not the structured part should be trained (only
    #' meaningful if pretrained), the structured mode should be pretrained 
    #' with a GAM and a list of lambdas to apply L2 regularization for linear 
    #' effects or to scale the penalization applied to smooth terms.
    #' Use this interface when tuning the structured part; tuning deep parameters
    #' can only be done via an outer loop.
    make_model = function(params = NULL) {
      tf$keras$backend$set_floatx(self$precision)
      self$make_ped(self$data, self$formulas, 
                    self$trafo_fct, 
                    self$cut, cr = self$cr, 
                    n_cr = self$n_cr, multimodal = self$multimodal,
                    train = TRUE)
      if (is.null(params)) {
        self$make_net(self$deep_architectures, self$formulas)
      } else {
        self$make_net(self$deep_architectures, self$formulas, params$trainable_struct, params$pretrain, params$lambdas)
      }
    },
    #' @import dplyr
    make_ped = function(data, formulas, trafo_fct, cut, cr, n_cr, multimodal, 
                        train = TRUE, partial_covar = NULL, partial_effect = NULL, 
                        Nout = 1000L) {
      if (self$ped & train) stop("Train data already transformed to PED.")
      if (multimodal) {
        data_unstruct <- data[2:length(data)]
        data <- data[[1]]
      }
      if (cr) {
        crs <- 1:n_cr
        data_ <- ped_data <- X <- X2 <- P <- sp <- pam <- coeffs <-  vector("list", n_cr)
      } else {
        crs <- 1
        data_ <- ped_data <- X <- X2 <- P <- sp <- pam <- coeffs <- vector("list", 1L)
      }
      if (!is.list(formulas)) {
        formulas <- rep(list(formulas), length(crs))
        self$formulas <- formulas
      }
      formulas <- decompose_formula(formulas, colnames(data))
      tabular_terms <- vector("list", length(data_))
      for (i in 1:length(data_)) {
        data_[[i]] <- data %>% mutate(status = ifelse(status == crs[i], 1, 0)) %>%
          as.data.frame()
        ped_data[[i]] <- as_ped(data_[[i]], trafo_fct, cut = cut)
        colnames(ped_data[[i]])[colnames(ped_data[[i]]) %in% c("tend", "ped_status")] <- c("time", "Y")
        if (train) {
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
        } else {
          partial_ <- !is.null(partial_covar) | !is.null(partial_effect)
          partial <- c(partial_covar, partial_effect)
          if (partial_) {
            Nout <- min(Nout, nrow(ped_data[[i]]))
            partial_type <- ifelse(!is.null(partial_covar), "covar", "effect")
            Nout <- (Nout %/% length(cut)) * length(cut)
            self$Nout <- Nout
            ped_data[[i]] <- ped_data[[i]][1:Nout, ]
            is_structured <- get_partial_type(partial, self$tabular_terms[[i]])
            covars <- get_partial_vars(self$tabular_terms[[i]], partial, partial_type, is_structured)
            if (is.null(is_structured)) is_structured <- TRUE
            if (partial_type == "covar") {
              if (length(covars) != 1L) stop("Only one Covar at a time!")
              ped_data[[i]][, !(colnames(ped_data[[i]]) %in% c(covars, "id"))] <- 0
              min_ <- min(ped_data[[i]][, colnames(ped_data[[i]]) == covars])
              max_ <- max(ped_data[[i]][, colnames(ped_data[[i]]) == covars])
              ped_data[[i]][, colnames(ped_data[[i]]) == covars] <- 
                seq(min_, max_, length.out = Nout)
              self$partial_domain <- ped_data[[i]][, colnames(ped_data[[i]]) == covars, drop = FALSE] 
            } else {
              ped_data[[i]][, !(colnames(ped_data[[i]]) %in% c(covars, "id"))] <- 0
              if (length(covars) == 1L) {
                mins <- min(ped_data[[i]][, colnames(ped_data[[i]] == covars)])
                maxs <- max(ped_data[[i]][, colnames(ped_data[[i]] == covars)])
                ped_data[[i]][, colnames(ped_data[[i]]) %in% covars] <- seq(mins[1], maxs[1], length.out = Nout)
                self$partial_domain <- ped_data[[i]][, colnames(ped_data[[i]]) == covars, drop = F]
              } else {
                mins <- sapply(ped_data[[i]][, colnames(ped_data[[i]]) %in% covars], min)
                maxs <- sapply(ped_data[[i]][, colnames(ped_data[[i]]) %in% covars], max)
                filled <- fill(ped_data[[i]], covars, mins, maxs)
                ped_data[[i]] <- filled$filled
                self$partial_domain <- filled$partial_domain
                self$Nout <- filled$length.out
              }
            }
          }
          mm <- predict(self$related_pamm[[i]], ped_data[[i]], type = "lpmatrix")
          mm <- cbind(mm, offset = 0)
          mm2 <- predict(self$processing_pam, ped_data[[i]], type = "lpmatrix")
          if (partial_) {
            if (partial_type == "effect" & !is_structured) {
              mm <- mm * 0
            }
          }
        }
        processed_terms <- colnames(mm)[!colnames(mm) %in% c("(Intercept)", colnames(ped_data[[i]]))]
        unprocessed_terms <- c(colnames(ped_data[[i]])[colnames(ped_data[[i]]) %in% colnames(mm)], "(Intercept)")
        processed_terms_ <- stringr::str_sub(processed_terms, 1,
                                             stringr::str_locate(processed_terms, stringr::fixed('.'))[, 1] - 1)
        pt <- processed_terms[is.na(processed_terms_)]
        processed_terms <- processed_terms[!is.na(processed_terms_)]
        processed_terms_ <- processed_terms_[!is.na(processed_terms_)]
        unprocessed_terms <- c(pt, unprocessed_terms)
        processed_terms_t <- split(processed_terms, f = processed_terms_)
        ut <- as.list(unprocessed_terms)
        names(ut) <- unprocessed_terms
        all_terms <- append(processed_terms_t, ut[names(ut) != "offset"])
        coeffs[[i]] <- vector("list", length(all_terms))
        names(coeffs[[i]]) <- names(all_terms)
        for (j in 1:length(all_terms)) {
          coeffs[[i]][[j]] <- pam[[i]]$coefficients[all_terms[[j]]]
        }
        X[[i]] <- vector("list", length(unique(processed_terms_)) + length(unprocessed_terms))
        P[[i]] <- vector("list", length(unique(processed_terms_)))
        names(X[[i]]) <- c(unique(processed_terms_), unprocessed_terms)
        names(P[[i]]) <- unique(processed_terms_)
        for (j in 1:length(processed_terms_t)) {
          if (length(processed_terms_t) > 0) {
            X[[i]][[j]] <- mm[, processed_terms_t[[j]], drop = FALSE]
            P[[i]][[j]] <- pam[[i]]$smooth[[j]]$S[[1]] * (pam[[i]]$sp[j] / nrow(ped_data[[i]]))
          }
        }
        k = j
        for (j in 1:length(unprocessed_terms)) {
          X[[i]][[k + j]] <- mm[, unprocessed_terms[[j]], drop = FALSE]
        }
        X2[[i]] <- vector("list", length(formulas[[i]]$deep_vars))
        no_deep <- TRUE
        if (length(X2[[i]]) > 0L) {
          no_deep <- FALSE
          for (j in 1:length(X2[[i]])) {
            use <- mm2[, colnames(mm2) %in% formulas[[i]]$deep_vars[[j]], drop = FALSE]
            X2[[i]][[j]] <- matrix(use, nrow = nrow(use), ncol = ncol(use))
            names(X2[[i]])[j] <- names(formulas[[i]]$deep_vars)[j]
          }
        }
        structured_covars <- lapply(names(processed_terms_t), get_structured_covars, ped = ped_data[[i]])
        names(structured_covars) <- names(processed_terms_t)
        tabular_terms[[i]] <- list(structured = structured_covars,
                                   deep = formulas[[i]]$deep_vars)
        if (!train) {
          if (partial_) {
            if (partial_type == "effect") {
              if (!is_structured) {
                for (mi in 1:length(X2)) {
                  if (names(X2[[i]])[[mi]] != partial) {
                    X2[[i]][[mi]] <- X2[[i]][[mi]] * 0
                  }
                }
              }
            }
          }
        }
      }
      X <- reshape(X, ped_data, cuts = self$cut)
      X2 <- reshape(X2, ped_data, cuts = self$cut)
      Y <- make_Y(ped_data, cuts = self$cut)
      if (train) {
        weights <- reshape_weights(ped_data)
        if (multimodal & !no_deep) {
          self$data = list(structured = X, deep = X2, unstructured = data_unstruct)
        } else if (!multimodal & !no_deep) {
          self$data = list(structured = X, deep = X2)
        } else {
          self$data = list(structured = X)
        }
        self$Y = tf$constant(Y, dtype = self$precision)
        self$P = P
        self$pam_coeffs = coeffs
        self$weights = tf$constant(weights, dtype = self$precision)
        self$related_pamm = pam
        self$processing_pam = pam2
      } else {
        no_deep <- !self$deep
        if (multimodal & !no_deep) {
          self$latest_test_data = list(structured = X, deep = X2, unstructured = data_unstruct)
        } else if (!multimodal & !no_deep) {
          self$latest_test_data = list(structured = X, deep = X2)
        } else {
          self$latest_test_data = list(structured = X) 
        }
      }
      if (train) {
        self$tabular_terms <- tabular_terms
      }
      self$ped <- TRUE
    },
    #' @import keras tensorflow dplyr
    make_net = function(deep_architectures, formulas, trainable_struct = TRUE, 
                        pretrain = FALSE, lambdas = NULL) {
      if (self$built) stop("Net already built.")
      P <- self$P
      data <- self$data
      t_ <- dim(data$structured[[1]][[1]])[2]
      structured <- structured_input <- vector("list", length(data[[1]])) 
      all_subnetnames <- all_covars <- c()
      for (i in 1:length(structured)) {
        d_a <- deep_architectures
        d_a <- d_a[sapply(paste0(names(d_a), stringr::fixed("(")), grepl, 
                          paste(deparse(formulas[[i]]), collapse = ""), fixed = T)]
        structured[[i]] <- structured_input[[i]] <- vector("list", length(data[[1]][[i]]))
        for (j in 1:length(data[[1]][[i]])) {
          subnetname <- names(data[[1]][[i]])[[j]]
          if (is.null(lambdas[[subnetname]])) {
            lf <- 0
          } else {
            lf <- lambdas[[subnetname]]
          }
          all_subnetnames <- c(all_subnetnames, paste0(subnetname, i, j))
          all_covars <- c(all_covars, names(data[[1]][[1]][[j]]))
          structured_input[[i]][[j]] <- layer_input(c(t_, dim(data[[1]][[1]][[j]])[3]), name = paste0(subnetname, i, j), dtype = self$precision)
          if (dim(data[[1]][[1]][[j]])[3] > 1L) {
            lf[lf == 0] <- 1
            if (!pretrain) {
              structured[[i]][[j]] <- structured_input[[i]][[j]] %>% 
                s_layer(P = tf$constant(P[[i]][[subnetname]], dtype = self$precision), 
                        trainable = trainable_struct, 
                        lambda_factor = lf) 
            } else {
              structured[[i]][[j]] <- structured_input[[i]][[j]] %>% 
                s_layer(P = tf$constant(P[[i]][[subnetname]], dtype = self$precision), 
                        lambda_factor = lf,
                        trainable = trainable_struct,
                        pretrained_weights = self$pam_coeffs[[i]][[subnetname]]) 
            }
          } else { 
            if (names(data[[1]][[i]])[j] != "offset") {
              if (!pretrain) {
                structured[[i]][[j]] <- structured_input[[i]][[j]] %>% 
                  layer_dense(1, use_bias = FALSE, 
                              trainable = trainable_struct,
                              kernel_regularizer = regularizer_l2(lf))
              } else {
                structured[[i]][[j]] <- structured_input[[i]][[j]] %>% 
                  layer_dense(1, use_bias = FALSE, 
                              trainable = trainable_struct,
                              kernel_regularizer = regularizer_l2(lf),
                              kernel_initializer = initializer_constant(self$pam_coeffs[[i]][[subnetname]]))
              }
            } else {
              structured[[i]][[j]] <- structured_input[[i]][[j]] %>% 
                layer_dense(1, trainable = FALSE, use_bias = FALSE,
                            kernel_initializer = initializer_constant(1))
            }
          }
        }
        structured[[i]] <- structured[[i]] %>% layer_add()
      }
      structured <- structured %>% layer_concatenate()
      structured_input <- unlist(structured_input)
      structured_data <- unname(data_unlist(data[[1]]))
      is_deep <- "deep" %in% names(data)
      is_unstructured <- "unstructured" %in% names(data)
      if (is_deep) {
        deep <- deep_input <- vector("list", length(data[["deep"]])) 
        for (i in 1:length(deep)) {
          deep[[i]] <- deep_input[[i]] <- vector("list", length(data[["deep"]][[i]]))
          for (j in 1:length(deep[[i]])) {
            subnetname <- names(data[["deep"]][[i]])[[j]]
            if (!(subnetname %in% names(d_a))) {
              break
            }
            all_subnetnames <- c(all_subnetnames, paste0(subnetname, i, j))
            all_covars <- c(all_covars, names(data[["deep"]][[i]][[j]]))
            deep_input[[i]][[j]] <- layer_input(c(t_, dim(data[["deep"]][[i]][[j]])[3]), dtype = self$precision, name = paste0(subnetname, i, j))
            arch <- d_a[[subnetname]]
            deep[[i]][[j]] <- deep_input[[i]][[j]] %>% arch
            deep[[i]][[j]] <- deep[[i]][[j]] %>% layer_add()
          }
          deep[[i]] <- deep[[i]][!sapply(deep[[i]], is.null)]
          deep_input[[i]] <- deep_input[[i]][!sapply(deep_input[[i]], is.null)]
          deep[[i]] <- deep[[i]] %>% layer_add()
        }
        deep <- deep %>% layer_concatenate()
        deep_input <- unlist(deep_input)
        deep_data <- unname(data_unlist(data[["deep"]]))
      } else {
        deep_input = NULL
        deep = NULL
        deep_data = NULL
        self$deep = FALSE
      }
      if (is_unstructured) {
        unstructured <- unstructured_input <- vector("list", length(data[["unstructured"]])) 
        for (i in 1:length(unstructured)) {
          subnetname <- names(data[["unstructured"]])[[i]]
          if (!(subnetname %in% names(d_a))) {
            break
          }
          all_subnetnames <- c(all_subnetnames, paste0(subnetname, i, j))
          all_covars <- c(all_covars, names(data[["unstructured"]][[i]]))
          unstructured_input[[i]] <- layer_input(dim(data[["unstructured"]][[i]])[-1], dtype = self$precision, name = paste0(subnetname, i, j))
          arch <- d_a[[subnetname]]
          unstructured[[i]] <- unstructured_input[[i]] %>% arch
          unstructured[[i]] <- unstructured[[i]] %>% 
            layer_repeat_vector(t_) %>%
            layer_dense(self$n_cr)
        }
        unstructured <- unstructured[!sapply(unstructured, is.null)]
        unstructured_input <- unstructured_input[!sapply(unstructured_input, is.null)]
        if (length(unstructured) > 0L) {
          unstructured <- unstructured %>% layer_add()
          unstructured_data <- unname(data_unlist(data[["unstructured"]]))
        } else { 
          unstructured <- unstructured_input <- unstructured_data <- NULL
          self$multimodal <- FALSE
        }
      } else {
        unstructured <- unstructured_input <- unstructured_data <- NULL
        self$multimodal <- FALSE
      }
      inputs <- smart_append(unlist(structured_input), unlist(deep_input), unstructured_input)
      XX <-lapply(smart_append(structured_data, deep_data, unstructured_data), 
                  tf$constant, dtype = self$precision)
      names(XX) <- all_subnetnames
      self$X <- XX
      output <- list(structured, deep, unstructured) 
      output <- output[!sapply(output, is.null)]
      output <- output %>% layer_add() %>% layer_activation("exponential") # fix
      self$subnets <- list(names = all_subnetnames, vars = all_covars)
      self$model <- keras_model(inputs, output)
      compile(self$model, loss = loss_poisson, optimizer = optimizer_adam(lr = self$lr), sample_weight_mode = "temporal", weighted_metrics = list())
      self$built <- TRUE
    },
    #' @description 
    #' 
    #' train the DEEPPAMM after constructing the model and data with $make_model
    #' @param epochs
    #' number of epochs used for training (integer+)
    #' @param batch_size
    #' batch size used for training (integer+)
    #' @param callbacks
    #' callbacks passed down to keras::fit
    #' @param val_split 
    #' validation split used for training (early stopping)
    train = function(epochs, batch_size, callbacks, val_split = 0.1) {
      fit(self$model, x = self$X, y = self$Y,
          sample_weight = self$weights,
          epochs = epochs, batch_size = batch_size, 
          callbacks = callbacks,
          validation_split = val_split,
          view_metrics = FALSE)
    },
    #' @description 
    #' Hazard Prediction
    #' 
    #' Predict hazards based on the intervals used for training on new data.
    #' 
    #' @param new_data
    #' New data provided for prediciton. MUST have identical format to training 
    #' data, also names MUST match.
    #' @param full
    #' logical, should the full follow up be predicted or only until event?
    #' @param verbose
    #' logical, verbosity passed to keras::predict
    predictHaz = function(new_data, full = TRUE, verbose = FALSE, time = "time", partial_covar = NULL, partial_effect = NULL, Nout = 1000L) {
      partial_ <- !is.null(partial_covar) | !is.null(partial_effect)
      if ((!is.null(partial_covar)) & !(is.null(partial_effect))) {
        stop("You can either generate partial predictions covariate wise or effect wise. Not both at the same time.")
      }
      if (!is.data.frame(new_data) & is.list(new_data) & !self$multimodal) {
        new_data <- new_data[[1]]
      }
      if (is.data.frame(new_data)) {
        n <- nrow(new_data)
        maxt = max(new_data[[time]])
        if (self$scale) { 
          new_data <- scale_(self$scaler, new_data)
        }
        if (full) new_data[[time]] <- rep(maxt, n)
      } else {
        n <- nrow(new_data[[1]])
        maxt = max(new_data[[1]][[time]])
        if (self$scale) new_data[[1]] <- scale_(self$scaler, new_data[[1]])
        if (full) new_data[[1]][[time]] <- rep(maxt, n)
      }
      self$make_ped(new_data, self$formulas, self$trafo_fct, self$cut, self$cr, 
                    self$n_cr, self$multimodal, train = FALSE, 
                    partial_covar = partial_covar, 
                    partial_effect = partial_effect)
      structured_data <- unname(data_unlist(self$latest_test_data[[1]]))
      dd <- self$latest_test_data[["deep"]]
      if (!is.null(dd)) {
        deep_data <- unname(data_unlist(dd)) 
      } else {
        deep_data <- NULL
      }
      ud <- self$latest_test_data[["unstructured"]]
      if (!is.null(ud)) {
        unstructured_data <- unname(data_unlist(ud))
        if (partial_) {
          for (mi in 1:length(unstructured_data)) {
            unstructured_data[[mi]] <- unstructured_data[[mi]] * 0
          }
        }
      } else {
        unstructured_data <- NULL
      }
      self$latest_test_data <- lapply(smart_append(structured_data, deep_data, unstructured_data), 
                                      tf$constant, dtype = self$precision)
      predict(self$model, self$latest_test_data, verbose = verbose)
    },
    #' @description 
    #' Standardized Cause Specific Hazard Prediction
    #' 
    #' Predict hazards based on new intervals on new data for a single cause.
    #' The whole follow-up is predicted.
    #' 
    #' @param new_data
    #' New data provided for prediction. MUST have identical format to training 
    #' data, also names MUST match.
    #' @param intervals
    #' integer+, intervals used for prediction, if not supplied: use the training
    #' intervals
    #' @param cause
    #' integer+, cause to predict cause-specific hazards for; in single risk case:
    #' 1L
    #' @import dplyr
    predictStdHaz = function(new_data, intervals = NULL, cause = 1L) {
      assert_numeric(intervals, lower = 0)
      if (is.null(intervals)) intervals <- self$cut
      id_var <- "id"
      time = self$time
      brks  <- self$cut
      if ( max(intervals) > max(brks) ) {
        stop("Cannot predict beyond the last time point used during model estimation.
        Check the 'times' argument.")
      }
      predicted_hazards <- self$predictHaz(new_data, full = TRUE)
      if (is.list(new_data) & !is.data.frame(new_data)) {
        newdata <- new_data[[1]]
      } else {
        newdata <- new_data
      }
      newdata$id = 1:nrow(newdata)
      ni <- dim(predicted_hazards)[1]
      pred_frame <- data.frame(id = rep(newdata[[id_var]], each = length(brks)), 
                               time = rep(brks, ni), pred = as.numeric(t(predicted_hazards[, , cause])))
      ped_times <- sort(unique(union(c(0, brks), intervals)))
      # extract relevant intervals only, keeps data small
      ped_times <- ped_times[ped_times <= max(intervals)]
      # obtain interval information
      ped_info <- get_intervals(brks, ped_times[-1])
      # add adjusted offset such that cumulative hazard and survival probability
      # can be calculated correctly
      ped_info[["intlen"]] <- c(ped_info[["times"]][1], diff(ped_info[["times"]]))
      newdata <- combine_df(ped_info, newdata)
      env_times <- intervals
      newdata <- inner_join(newdata, pred_frame, by = c("id" = "id", "tend" = "time"))
      newdata[["intlen"]] <- recompute_intlen(newdata)
      newdata <- newdata %>%
        arrange(.data$id, .data$times) %>%
        group_by(.data$id) %>%
        mutate(stdhaz = (.data$pred * .data$intlen)) %>%
        ungroup() %>%
        dplyr::filter(.data[["times"]] %in% env_times)
      stdhaz <- matrix(newdata$stdhaz, nrow = ni, 
                       ncol = length(intervals), byrow = TRUE)
      stdhaz
    },
    #' @description 
    #' Cause Specific Cumulative Hazard Prediction
    #' 
    #' Predict cumulative hazards based on new intervals on new data for a 
    #' single cause. The whole follow-up is predicted.
    #' 
    #' @param new_data
    #' New data provided for prediction. MUST have identical format to training 
    #' data, also names MUST match.
    #' @param intervals
    #' integer+, intervals used for prediction, if not supplied: use the training
    #' intervals
    #' @param cause
    #' integer+, cause to predict cause-specific hazards for; in single risk case:
    #' 1L
    predictCumHaz = function(new_data, intervals = NULL, cause = 1L) {
      if (is.null(intervals)) intervals = self$cut
      stdhaz <- self$predictStdHaz(new_data, intervals, cause = cause)
      tf$cumsum(stdhaz, axis = 1L)$numpy()
    },
    #' @description 
    #' (Cause Specific) Survival Probability Prediction
    #' 
    #' Predict survival probabilities based on new intervals on new data for a 
    #' single cause. The whole follow-up is predicted.
    #' In most cases, you rather want to use predictCIFs in competing risks 
    #' settings. Hence, we recommend to use this function just for single risk.
    #' @param new_data
    #' New data provided for prediction. MUST have identical format to training 
    #' data, also names MUST match.
    #' @param intervals
    #' integer+, intervals used for prediction, if not supplied: use the training
    #' intervals
    #' @param cause
    #' integer+, cause to predict cause-specific hazards for; in single risk case:
    #' 1L
    predictSurvProb = function(new_data, intervals = NULL, cause = 1L) {
      if (is.null(intervals)) intervals = self$cut
      exp(-self$predictCumHaz(new_data, intervals, cause = cause))
    },
    #' @description 
    #' Cumulative Incidence Functions (CIFs) Predictions
    #' 
    #' Predict CIFs based on new intervals on new data for a all causes.
    #' The whole follow-up is predicted.
    #' 
    #' @param new_data
    #' New data provided for prediction. MUST have identical format to training 
    #' data, also names MUST match.
    #' @param intervals
    #' integer+, intervals used for prediction, if not supplied: use the training
    #' intervals
    predictCIFs = function(new_data, intervals) {
      stdhaz <- CIFs <- cumhaz <- vector("list", self$n_cr)
      for (j in 1:self$n_cr) {
        stdhaz[[j]] <- self$predictStdHaz(new_data, intervals, cause = j)
        cumhaz[[j]] <- self$predictCumHaz(new_data, intervals, cause = j)
      }
      overall_survival <- exp(-Reduce("+", cumhaz))
      lagged_os <- cbind(1, overall_survival[, 1:(ncol(overall_survival) - 1)])
      for (j in 1:self$n_cr) {
        CIFs[[j]] <- tf$cumsum(stdhaz[[j]] * lagged_os, axis = 1L)$numpy()
      }
      
      CIFs
    }, 
    predict_partial = function(partial_covar = NULL, partial_effect = NULL, time = "time", Nout = 1000L) {
      if (is.null(partial_covar) & is.null(partial_effect)) {
        stop("You must supply either the covariable or effect you are interested in.")
      } 
      new_data <- self$raw_data
      hazards <- self$predictHaz(new_data, full = TRUE, verbose = FALSE,
                                 time = time, partial_covar = partial_covar, 
                                 partial_effect = partial_effect) 
      if (!self$cr) {
        haz <- data.frame(haz = as.numeric(t(hazards[, , 1])),
                          self$partial_domain)
      }
      else {
        haz <- vector("list", length(dim(hazards)[3L]))
        for (i in 1:length(haz)) {
          haz[[i]] <- data.frame(haz = as.numeric(t(hazards[, , i])),
                                 self$partial_domain)
        }
      }
      haz[1:self$Nout, ]
    }
  )
)
