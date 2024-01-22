reshape_weights <- function(ped) {
  ped_ <- ped[[1]]
  ids <- ped_$id
  tm <- unique(ped_[["time"]])
  res <- matrix(0, nrow = length(unique(ids)), ncol = length(tm))
  ped_ind <- 1
  for (i in 1:nrow(res)) {
    current_len <- sum(i == ids)
    res[i, 1:(current_len)] <- 1
    ped_ind <- ped_ind + current_len + 1
  }
  res
}

make_Y <- function(ped) {
  res <- array(0, dim = c(length(unique(ped[[1]]$id)), length(unique(ped[[1]]$time)), length(ped)))
  for (k in 1:length(ped)) {
    ped_ <- ped[[k]]
    ids <- ped_$id
    tm <- attributes(ped_)$trafo_args$cut
    ped_ind <- 1
    for (i in 1:dim(res)[1]) {
      current_len <- sum(i == ids)
      res[i, current_len, k] <- ped[[k]]$Y[ped_ind + current_len - 1]
      ped_ind <- ped_ind + current_len
    }
  }
  res
}

reshape <- function(X, ped, cuts = NULL) {
  if (is.null(cuts)) cuts <- length(unique(ped[[i]]$time))
  res <- vector("list", length(X))
  for (i in 1:length(X)) {
    ids <- ped[[i]]$id
    res[[i]] <- vector("list", length(X[[i]]))
    for (j in 1:length(X[[i]])) {
      res[[i]][[j]] <- array(0, dim = c(length(unique(ped[[i]]$id)), cuts, ncol(X[[i]][[j]])))
      ped_ind <- 1
      for (k in 1:nrow(res[[i]][[j]])) {
        current_len <- sum(k == ids)
        res[[i]][[j]][k, 1:(current_len),] <- X[[i]][[j]][ped_ind:(ped_ind + current_len - 1),]
        ped_ind <- ped_ind + current_len
      }
    }
  }
  names(res[[i]]) <- names(X[[i]])
  res
}