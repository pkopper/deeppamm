data_f <- data.frame(status = sample(c(0, 1, 2), 14, T), time = 1:14, 
                     x1 = runif(14), x2 = runif(14))

DEEPPAMM <- deeppamm$new(formulas = Y ~ s(time) + x1 + s(x2) + deep(time, x1) + deep(x2) + image(image),
                         deep_architectures = 
                           list(deep = function(x) x %>% layer_dense(10, activation = "relu") %>% layer_dense(1),
                                image = function(x) x %>% layer_dense(32, "relu") %>% layer_flatten %>% layer_dense(8)),
                         data = list(data_f, image = array(runif(14 * 5 * 5), c(14, 5, 5, 1))),
                         trafo_fct = Surv(time, status) ~.,
                         cut = 1:12,
                         lr = 0.001)
DEEPPAMM$make_model()
Y_pred <- DEEPPAMM$model$predict(lapply(DEEPPAMM$X, tf$constant))
loss_poisson(Y_pred,tf$constant(DEEPPAMM$Y))

DEEPPAMM$model$inputs
lapply(DEEPPAMM$X, dim)
DEEPPAMM$weights
DEEPPAMM$train(10, 4, callbacks = list(callback_early_stopping(patience = 10)))
ndata <- list(data_f, image = array(runif(14 * 5 * 5), c(14, 5, 5, 1)))
DEEPPAMM$predictHaz(new_data = ndata)
a <- DEEPPAMM$predictStdHaz(ndata, seq(0.01, 10, by = 0.1))
b <- DEEPPAMM$predictCumHaz(ndata, seq(0.01, 10, by = 0.1))
sp <- DEEPPAMM$predictSurvProb(ndata, seq(0.01, 10, by = 0.1))
CIFs <- DEEPPAMM$predictCIFs(ndata, seq(0.01, 10, by = 0.1))
plot(seq(0.01, 10, by = 0.1), CIFs[[2]][1,])
