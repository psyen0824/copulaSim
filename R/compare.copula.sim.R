#' The function to do the comparison between input dataset and simulated datasets.
#'
#' @param object A copula.sim object for the comparison.
#' @return Returned simulated data
#' @export
#' @importFrom stats sd
#' @importFrom dplyr left_join
compare.copula.sim <- function(object) {
  if (!("copula.sim" %in% class(object))) {
    stop("Only accept copula.sim object.")
  }

  empir.stat.df <- object$data.transform %>%
    group_by(.data$arm, .data$col.num) %>%
    summarise(
      empir.sample = n(),
      empir.mean = mean(.data$data.input),
      empir.sd = sd(.data$data.input),
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
      simu.RB = round(mean(.data$simu.RB),4),
      simu.SB = round(mean(.data$simu.SB),4),
      simu.RMSE = round(sqrt(mean(.data$simu.RMSE)),4)
    )

  simu.ci.df <- object$data.simul %>%
    group_split(.data$sim.id) %>%
    lapply(function(sim.df){
      sim.df %>%
        group_by(.data$arm, .data$col.num) %>%
        summarise(
          simu.mean = mean(.data$data.sim),
        )
    }) %>%
    bind_rows %>%
    group_by(.data$arm, .data$col.num) %>%
    summarise(
      simu.lower.ci = round(quantile(.data$simu.mean, probs = 0.025),4),
      simu.upper.ci = round(quantile(.data$simu.mean, probs = 0.975),4),
      simu.length.ci = simu.upper.ci - simu.lower.ci
    )


  comparison.result.df <- empir.stat.df %>%
    left_join(size.stat.df, c("arm", "col.num")) %>%
    left_join(simu.diff.df, c("arm", "col.num")) %>%
    left_join(simu.ci.df, c("arm", "col.num")) %>%
    mutate(marginal.name = colnames(object$data.input)) %>%
    select(
      .data$marginal.name, .data$arm, .data$empir.sample, .data$simu.sample,
      .data$n.simu, .data$empir.mean, .data$simu.mean, .data$empir.sd,
      .data$simu.sd, .data$simu.RB, .data$simu.SB, .data$simu.RMSE,
      .data$simu.lower.ci, .data$simu.upper.ci, .data$simu.length.ci
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
