context("Competing risks example")
library(deeppamm)
sim_pexp_cr <- function(formula, data, cut) {
  
  Form <- Formula(formula)
  F_rhs <- attr(Form, "rhs")
  l_rhs <- length(F_rhs)
  seq_rhs <- seq_len(l_rhs)
  
  data <- data %>%
    mutate(
      id     = row_number(),
      time   = max(cut),
      status = 1)
  
  # construct eta for time-constant part
  ped  <- split_data(
    formula = Surv(time, status)~.,
    data    = select_if (data, is_atomic),
    cut     = cut,
    id      = "id") %>%
    rename("t" = "tstart")
  
  # calculate cause specific hazards
  for(i in seq_rhs) {
    ped[[paste0("hazard", i)]] <-  exp(eval(F_rhs[[i]], ped))
  }
  ped[["rate"]] <- reduce(ped[paste0("hazard", seq_rhs)], `+`)
  
  # simulate survival times
  sim_df <- ped %>%
    group_by(id) %>%
    mutate(
      time   = msm::rpexp(rate = .data$rate, t = .data$t),
      status = 1L * (.data$time <= max(cut)),
      time   = pmin(.data$time, max(cut))) %>%
    dplyr::filter(.data$t < .data$time & .data$time <= .data$tend)
  sim_df$type <- apply(sim_df[paste0("hazard", seq_rhs)], 1,
                       function(probs)
                         sample(seq_rhs, 1, prob = probs))
  
  sim_df %>%
    mutate(type = ifelse(.data$status == 1, .data$type, 0)) %>%
    select(-one_of(c("t", "tend", "interval", "offset", "ped_status", 
                     "rate")))
}

predictEventProb_pamm <- function(pamms, times, newdata, ped, breaks) {
  trafo_args <- attr(ped, "trafo_args")
  id_var <- trafo_args[["id"]]
  brks       <- breaks#trafo_args[["cut"]]
  if ( max(times) > max(brks) ) {
    stop("Cannot predict beyond the last time point used during model estimation.
        Check the 'times' argument.")
  }
  n_row <- nrow(newdata)
  env_times <- times
  ped_times <- sort(unique(union(c(0, brks), times)))
  # extract relevant intervals only, keeps data small
  ped_times <- ped_times[ped_times <= max(times)]
  # obtain interval information
  ped_info <- get_intervals(brks, ped_times[-1])
  # add adjusted offset such that cumulative hazard and survival probability
  # can be calculated correctly
  ped_info[["intlen"]] <- c(ped_info[["times"]][1], diff(ped_info[["times"]]))
  newdata <- combine_df(ped_info, newdata)
  newdata[["intlen"]] <- recompute_intlen(newdata)
  hazards <- list(predict(pamms[[1]], newdata), predict(pamms[[2]], newdata))
  newdata$pred1 <- exp(hazards[[1]])
  newdata$pred2 <- exp(hazards[[2]])
  newdata <- newdata %>%
    arrange(.data$id, .data$times) %>%
    group_by(.data$id) %>%
    mutate(overall_survival = exp(-cumsum((.data$pred1 + .data$pred2) * .data$intlen))) %>%
    mutate(cif1 = cumsum(.data$pred1 * lag(.data$overall_survival, default = 1) * .data$intlen),
           cif2 = cumsum(.data$pred2 * lag(.data$overall_survival, default = 1) * .data$intlen)) %>%
    ungroup() %>%
    dplyr::filter(.data[["times"]] %in% env_times)
  cif1 <- matrix(newdata$cif1, nrow = n_row, 
                 ncol = length(times), byrow = TRUE)
  cif2 <- matrix(newdata$cif2, nrow = n_row, 
                 ncol = length(times), byrow = TRUE)
  list(cif1 = cif1, cif2 = cif2)
}

library(pammtools)
set.seed(18)
n <- 800
Te <- 90
Tn <- 30
data <- data.frame(x1 = runif(n, 0, 2), x2 = runif(n, 0, 2), x3 = runif(n, -1, 1))
f0_1 <- function(t) {
  dgamma(t, 8, 2) * 6
}
f0_2 <- function(t) {
  - t / 20
}
form <- ~ -4 + f0_1(t) -1 * x1 + 2 * sqrt(x2) - 1 * x1 * x3 | 
  -2 + f0_2(t) - sqrt(x1 + x2) + sin(3 * x3)
sim_data <- sim_pexp_cr(form, data = data, cut = seq(0.1, Tn, length.out = 250)) 
data <- sim_data %>% ungroup() %>%
  mutate(cens = runif(n, 0, 2 * Tn),
         status = type, type = NULL,
         status = ifelse(cens < time, 0, status),
         time = ifelse(cens < time, cens, time), 
         cens = NULL, hazard1 = NULL, hazard2 = NULL)

testthat::test_that("CR",{
  set.seed(1223)
  tr <- sample(1:nrow(data), round(0.8 * nrow(data)))
  te <- setdiff(1:nrow(data), tr)
  DEEPPAMM <- deeppamm$new(formulas = Y ~ s(time) + x1 + x2 + x3 + deep(x1, x2, x3),
                           deep_architectures = 
                             list(deep = function(x) x %>% 
                                    layer_dense(32, activation = "relu") %>% 
                                    layer_dropout(0.5) %>%
                                    layer_dense(32, activation = "relu") %>% 
                                    layer_dropout(0.5) %>%
                                    layer_dense(16, activation = "relu") %>% 
                                    layer_dense(1)),
                           data = data[tr, ] %>% select(-id),
                           trafo_fct = Surv(time, status) ~.,
                           cut = c(0.1, 0.5, 1:30),
                           lr = 0.001,
                           scale = FALSE)
  DEEPPAMM$make_model(make_params())
  DEEPPAMM$train(epochs = 1000, batch_size = 128, callbacks = list(
    callback_reduce_lr_on_plateau(min_lr = 0.0001, factor = 0.5, patience = 50L),
    callback_early_stopping(patience = 100L)))
  times <- seq(0.01, 20, length.out = 250L)
  y_pred <- DEEPPAMM$predictCIFs(data[te, ] %>% select(-id), intervals = times)
  
  ped <- as_ped(Surv(time, status) ~ ., data = data[tr, ] %>% select(-id) %>% as.data.frame(), 
                cut = c(0.1, 0.5, 1:30), combine = FALSE)
  pam_csh <- map(ped, ~ pamm(ped_status ~ s(tend) + x1 + x2 + x3, data = .x))
  pam_CIF <- predictEventProb_pamm(pamms = pam_csh, times = times, 
                                   newdata = data[te, ], ped = ped, 
                                   breaks = c(0.1, 0.5, 1:30))
  
  pec1 <- pec(
    list(PAM = pam_CIF[[1]], Ours = y_pred[[1]]),
    Hist(time, status) ~ 1, # formula for IPCW
    data = data[te, ], # new data not used for model fit
    times = times,
    start = 0.01,
    exact = FALSE
  )
  pec2 <- pec(
    list(PAM = pam_CIF[[2]], Ours = y_pred[[2]]),
    Hist(time, status) ~ 1, # formula for IPCW
    data = data[te, ], # new data not used for model fit
    times = times,
    start = 0.01,
    exact = FALSE,
    cause = 2L
  )
})
