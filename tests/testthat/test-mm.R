context("Multimodal example")
library(deeppamm)
library(pammtools)
set.seed(18)
n <- 3200
Te <- 90
Tn <- 30
mnist <- keras::dataset_mnist()
mnist <- list(image = mnist$train$x[1:n,,], label = mnist$train$y[1:n])
mnist$image <- array_reshape(mnist$image, c(n, 28, 28, 1))
mnist$image <- mnist$image / 256
data <- data.frame(x1 = runif(n, 0, 2), x2 = runif(n, 0, 2), x3 = runif(n, -1, 1), label = mnist$label)
f0 <- function(t) {
  dgamma(t, 8, 2) * 6 + 0.75 * sin(t / 10)
}
form <- ~ -3.5 + f0(t) -2 * x1 + 3 * sqrt(x2) - 0.1 * x1 * x3 + 1.25 * as.integer(label > 5)
sim_data <- sim_pexp(form, data = data, cut = seq(0.1, Tn, length.out = 250)) 
data <- sim_data %>% 
  mutate(cens = runif(n, 0, 2 * Tn),
         status = ifelse(cens < time, 0, status),
         time = ifelse(cens < time, cens, time), 
         cens = NULL, id = NULL)

data <- list(data, image = mnist$image)


testthat::test_that("MM example", {
  set.seed(1223)
  data <- subset_mm(data)
  DEEPPAMM <- deeppamm$new(formulas = Y ~ s(time) + x1 + deep(x1, x2, x3),# + image(image),
                           deep_architectures = 
                             list(deep = function(x) x %>% 
                                    layer_dense(32, activation = "relu") %>% 
                                    layer_dropout(0.5) %>%
                                    layer_dense(32, activation = "relu") %>% 
                                    layer_dropout(0.5) %>%
                                    layer_dense(16, activation = "relu") %>% 
                                    layer_dense(1),
                                  image = function(x) x %>% 
                                    layer_conv_2d(32, kernel_size = c(3, 3), activation = "relu") %>%
                                    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
                                    layer_conv_2d(64, kernel_size = c(3, 3), activation = "relu") %>%
                                    layer_max_pooling_2d(pool_size = c(2, 2)) %>%
                                    layer_flatten() %>%
                                    layer_dropout(0.5)),
                           data = data[[1]],
                           trafo_fct = Surv(time, status) ~.,
                           cut = c(0.1, 0.5, 1:30),
                           lr = 0.01,
                           scale = FALSE)
  DEEPPAMM$make_model(make_params(TRUE, TRUE))
  DEEPPAMM$train(epochs = 1000, batch_size = 128, callbacks = list(
    callback_reduce_lr_on_plateau(min_lr = 0.0001, factor = 0.5, patience = 10L),
    callback_early_stopping(patience = 20L)))
  times <- seq(0.01, 20, length.out = 250L)
  y_pred <- DEEPPAMM$predictSurvProb(data[[2]], intervals = times)
  
  ped <- as_ped(Surv(time, status) ~ ., data = data[[1]][[1]], cut = c(0.1, 0.5, 1:30))
  pam <- pamm(ped_status ~ s(tend) + x1 + x2 + x3, data = ped)
  pam_o <- pamm(ped_status ~ s(tend) + x1 + s(x2) + x1:x3, data = ped)
  pam_oi <- pamm(ped_status ~ s(tend) + x1 + s(x2) + x1:x3 + label, data = ped)
  
  pec <- pec(
    list(PAM = pam, optimal_no_image = pam_o, optimal = pam_oi, Ours = y_pred),
    Surv(time, status) ~ 1, # formula for IPCW
    data = data[[2]][[1]], # new data not used for model fit
    times = times,
    start = 0.01,
    exact = FALSE
  )
})
