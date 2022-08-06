#' Simulating new multivariate datasets with shifted mean vector from existing empirical data
#'
#' @param data.input,id.vec,arm.vec,n.patient,n.simulation,seed Please refer to the function \link{copula.sim}.
#' @param validation.type,validation.sig.lvl,rmvnorm.matrix.decomp.method,verbose Please refer to the function \link{copula.sim}.
#' @param shift.vec.list A list of numeric vectors to specify the mean-shifted values for new arms.
#' @return Please refer to the function \link{copula.sim}.
#' @export
#' @author Pei-Shan Yen, Xuemin Gu, Jenny Jiao, Jane Zhang
#' @examples
#'
#' library(copulaSim)
#'
#' ## Generate Empirical Data
#'  # Assume that the single-arm, 3-dimensional empirical data follows multivariate normal data
#' library(mvtnorm)
#' arm1 <- rmvnorm(n = 80, mean = c(10,10.5,11), sigma = diag(3) + 0.5)
#' test_data <- as.data.frame(cbind(1:80, rep(1,80), arm1))
#' colnames(test_data) <- c("id", "arm", paste0("time_", 1:3))
#'
#' ## Generate 1 simulated datasets with one empirical arm and two new-arm.
#' ## The mean difference between empirical arm and
#'  # (i) the 1st new arm is assumed to be 2.5, 2.55, and 2.6 at each time point
#'  # (ii) the 2nd new arm is assumed to be 4.5, 4.55, and 4.6 at each time point
#' new.arm.copula.sim(data.input = test_data[,-c(1,2)],
#'   id.vec = test_data$id, arm.vec = test_data$arm,
#'   n.patient = 100 , n.simulation = 1, seed = 2022,
#'   shift.vec.list = list(c(2.5,2.55,2.6), c(4.5,4.55,4.6)))

new.arm.copula.sim <- function(data.input,
                               id.vec,
                               arm.vec,
                               shift.vec.list,
                               n.patient,
                               n.simulation,
                               seed = NULL,
                               validation.type = "none",
                               validation.sig.lvl = 0.05,
                               rmvnorm.matrix.decomp.method = "svd",
                               verbose = TRUE) {
  if (length(unique(arm.vec)) > 1){
    stop("Input data must be single ARM.")
  }

  if (!("list" %in% class(shift.vec.list)) || !all(sapply(shift.vec.list, is.numeric))) {
    stop("shift.vec.list must be a list of numeric vector.")
  }

  if (is.null(dim(data.input))) {
    stop("data.input must be a numeric data.frame or matrix.")
  }

  data.input <- as.matrix(data.input)
  if (any(sapply(shift.vec.list, length) == 1)) {
    shift.vec.list <- lapply(shift.vec.list, function(v){
      if (length(v) == 1) {
        rep(v, ncol(data.input))
      } else {
        v
      }
    })
  }

  if (any(sapply(shift.vec.list, length) != ncol(data.input))) {
    stop("The element of shift.vec.list must be a vector with length 1 or ncol(data.input).")
  }

  first_arm <-
    copula.sim(
      data.input,
      id.vec,
      arm.vec,
      n.patient,
      n.simulation,
      seed,
      validation.type,
      validation.sig.lvl,
      rmvnorm.matrix.decomp.method,
      verbose
    )
  data.new.arm <- mapply(function(i, s){
      mapply(function(simu.id, data.simu.mat){
        data.simu <- copula.sim(
          sweep(data.simu.mat[[1]], 2, s, `+`),
          seq.int(1L, nrow(data.simu.mat[[1]])),
          rep(i + 1L, nrow(data.simu.mat[[1]])),
          n.patient,
          1L,
          seed,
          validation.type,
          validation.sig.lvl,
          rmvnorm.matrix.decomp.method,
          verbose
        )$data.simul
        data.simu$sim.id <- simu.id
        return(data.simu)
    }, seq_along(extract.data.sim(first_arm)), extract.data.sim(first_arm), SIMPLIFY = FALSE) %>% bind_rows
  }, seq_along(shift.vec.list), shift.vec.list, SIMPLIFY = FALSE) %>%
    bind_rows

  res <- list(
    data.input = first_arm$data.input,
    data.transform = first_arm$data.transform,
    data.simul = bind_rows(first_arm$data.simul, data.new.arm)
  )
  class(res) <- c("copula.sim", "copula.sim.new.arm")
  return(res)
}
