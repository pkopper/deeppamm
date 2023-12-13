library(testthat)
library(pec)

if (reticulate::py_module_available("tensorflow") & 
    reticulate::py_module_available("keras")){
  test_check("deeppamm")
}