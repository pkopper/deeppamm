context("Single Risk Example")
library(deeppamm)
library(pammtools)
set.seed(18)
n <- 600
Te <- 90
Tn <- 30
data <- data.frame(x1 = runif(n, 0, 2), x2 = runif(n, 0, 2), x3 = runif(n, -1, 1))
f0 <- function(t) {
  dgamma(t, 8, 2) *6
}
form <- ~ -3.5 + f0(t) -2 * x1 + 3 * sqrt(x2) - 0.1 * x1 * x3
sim_data <- sim_pexp(form, data = data, cut = seq(0.1, Tn, length.out = 250)) 
data <- sim_data %>% 
  mutate(cens = runif(n, 0, 2* Tn),
         status = ifelse(cens < time, 0, status),
         time = ifelse(cens < time, cens, time), 
         cens = NULL)


testthat::test_that("Single Risk", {
  expect_no_error({set.seed(1223)
  tr <- sample(1:nrow(data), round(0.8 * nrow(data)))
  te <- setdiff(1:nrow(data), tr)
  DEEPPAMM <- deeppamm$new(formulas = Y ~ s(time) + x1 + x2 + x3 + deep(x1, x2, x3),
                           deep_architectures = 
                             list(deep = function(x) x %>% 
                                    #layer_batch_normalization() %>%
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
  y_pred <- DEEPPAMM$predictSurvProb(data[te, ] %>% select(-id), intervals = times)
  
  ped <- as_ped(Surv(time, status) ~ ., data = data[tr, ] %>% select(-id), cut = c(0.1, 0.5, 1:30))
  pam <- pamm(ped_status ~ s(tend) + x1 + x2 + x3, data = ped)
  pam_o <- pamm(ped_status ~ s(tend) + x1 + s(x2) + x1:x3, data = ped)
  
  pec <- pec(
    list(PAM = pam, optimal = pam_o, Ours = y_pred),
    Surv(time, status) ~ 1, # formula for IPCW
    data = data[te, ], # new data not used for model fit
    times = times,
    start = 0.01,
    exact = FALSE
  )}, message = "Error")
})
