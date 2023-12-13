check_formulas <- function(formulas) {
  if (is.list(formulas)) {
    t <- all(sapply(formulas, is_formula))
  } else {
    t <- is_formula(formulas)
  }
  if (!t) stop("Argument <formulas> not compliant.")
}

check_data <- function(data) {
  if (is.list(data) & !is.data.frame(data)) {
    t <- TRUE
  } else {
    t <- all(names(data)[-1] != "")
  }
  if (!t) stop("<data> must be either a data.frame or a list of modalities, where all but the first entry must be named.")
}

check_trafo_fct <- function(trafo_fct) {
  assert_formula(trafo_fct)
}

check_architectures <- function(deep_architectures) {
  no_list <- !is.list(deep_architectures)
  if (!no_list) {
    no_names_at_all <- is.null(names(deep_architectures))
    no_names_partly <- any(names(deep_architectures) == "")
    no_names <- no_names_at_all | no_names_partly
  }
  if (no_list | no_names) {
    stop("<deep_architectures> must be supplied as a named list")
  }
}