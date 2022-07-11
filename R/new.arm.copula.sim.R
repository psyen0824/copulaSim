#' The function to simulate Copula datasets with new arms
#'
#' @param data.input,id.vec,arm.vec,n.patient,n.simulation,seed,std.norm.lb,std.norm.ub,validation.type,validation.sig.lvl,rmvnorm.matrix.decomp.method,verbosee Please refer to the function \link{copula.sim}.
#' @param shift.vec.list A list of numeric vectors to specify the mean-shifted values for new arms.
#' @return Please refer to the function \link{copula.sim}.
#' @export
new.arm.copula.sim <- function(data.input,
                               id.vec,
                               arm.vec,
                               shift.vec.list,
                               n.patient,
                               n.simulation,
                               seed = NULL,
                               validation.type = "none",
                               validation.sig.lvl = 0.05,
                               std.norm.lb = -3,
                               std.norm.ub = 3,
                               rmvnorm.matrix.decomp.method = "svd",
                               verbose = TRUE) {
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
      std.norm.lb,
      std.norm.ub,
      rmvnorm.matrix.decomp.method,
      verbose
    )
  return(first_arm)
}
