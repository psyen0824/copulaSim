#' Performing the comparison between empirical data and multiple simulated datasets.
#'
#' @param object A copula.sim object for the comparison.
#' @return Returned the comparison of marginal parameter and covariance.
#'  1. mean.comparison: comparison between empirical mean and average value of simulated mean.
#'     (1) simu.mean: average value of simulated mean
#'     (2) simu.sd: average value of simulated standard error
#'     (3) simu.mean.RB: relative bias for marginal mean
#'     (4) simu.mean.SB: standardized bias for marginal mean
#'     (5) simu.mean.RMSE: root mean square error for marginal mean
#'  2. cov.comparison: comparison between empirical covariance and average value of simulated covariance
#'
#' @export
#' @author Pei-Shan Yen, Xuemin Gu
#' @importFrom stats sd quantile
#' @importFrom dplyr left_join
compare.copula.sim <- function(object) {
  if (!("copula.sim" %in% class(object))) {
    stop("Only accept copula.sim object.")
  }

  empir.stat.df <- object$data.transform %>%
    group_by(.data$arm, .data$col.num) %>%
    summarise(
      empir.sample = n(),
      empir.mean = round(mean(.data$data.input), 4),
      empir.sd = round( sd(.data$data.input), 4),
    )

  size.stat.df <- object$data.simul %>%
    group_by(.data$arm, .data$col.num) %>%
    summarise(
      simu.sample = max(.data$id),
      n.simu = max(.data$sim.id),
    )

  simu.diff.df <- object$data.simul %>%
    group_split(.data$sim.id) %>%
    lapply(function(sim.df){
      sim.df %>%
        group_by(.data$arm, .data$col.num) %>%
        summarise(
          simu.mean = mean(.data$data.sim),
          simu.sd = sd(.data$data.sim),
        ) %>%
        left_join(empir.stat.df, c("arm", "col.num")) %>%
        mutate(
          simu.RB = (.data$simu.mean - .data$empir.mean)/.data$empir.mean*100,
          simu.SB = abs(.data$simu.mean - .data$empir.mean)/.data$simu.sd*100,
          simu.RMSE = (.data$simu.mean - .data$empir.mean)^2
        )
    }) %>%
    bind_rows %>%
    group_by(.data$arm, .data$col.num) %>%
    summarise(
      simu.mean = round(mean(.data$simu.mean),4),
      simu.sd = round(mean(.data$simu.sd),4),
      simu.mean.RB = round(mean(.data$simu.RB),4),
      simu.mean.SB = round(mean(.data$simu.SB),4),
      simu.mean.RMSE = round(sqrt(mean(.data$simu.RMSE)),4)
    )

  comparison.result.df <- empir.stat.df %>%
    left_join(size.stat.df, c("arm", "col.num")) %>%
    left_join(simu.diff.df, c("arm", "col.num")) %>%
    mutate(marginal.name = colnames(object$data.input)) %>%
    select(
      .data$marginal.name, .data$arm, .data$empir.sample, .data$simu.sample, .data$n.simu,
      .data$empir.mean, .data$simu.mean,
      .data$simu.mean.RB, .data$simu.mean.SB, .data$simu.mean.RMSE,
      .data$empir.sd, .data$simu.sd
    ) %>% ungroup

  data.transform.split <- object$data.transform %>%
    group_by(.data$arm) %>%
    group_split
  empir.cov <- data.transform.split %>%
    lapply(function(df){
      cov.mat <- df %>%
        arrange(.data$col.num, .data$id) %>%
        (function(temp.df) cov(matrix(temp.df$data.input, ncol = max(temp.df$col.num),
                                      dimnames = list(NULL, colnames(object$data.input)))))
    }) %>% set_names(sprintf("arm=%i", sapply(data.transform.split, function(df) df$arm[1])))

  sim.data.mat.list <- extract.data.sim(object)
  simu.cov <- sim.data.mat.list %>%
    lapply(function(arm.list) lapply(arm.list, cov)) %>%
    (function(zz){
      arm.name.list <- names(zz[[1]])
      lapply(arm.name.list, function(arm) {
        Reduce(function(x, y) x + y, lapply(zz, function(l) l[[arm]])) / length(zz)
      })
    }) %>% set_names(names(sim.data.mat.list[[1]]))

  return(list(
    mean.comparison = comparison.result.df,
    cov.comparison = list(empir.cov = empir.cov, simu.cov = simu.cov)
  ))
}
