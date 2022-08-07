# helper function to recover the simulated data to the range of empirical via quantile transformation
empirical.ecdf.inv <-
  function(sim.mvnorm.data,
           arm.id,
           data.transform,
           mean.vec,
           cov.mat) {
    data.df.split.arm <- data.transform %>% filter(.data$arm == arm.id) %>% group_split(.data$col.num)
    mapply(
      function(df, s.data, m, sd) {
        inv.data <-
          ecdf.inv(x = df$data.input.transformed, p = pnorm(s.data, mean = m, sd = sd), sort.flag = FALSE)
        if (df$data.is.integer[1]) { inv.data <- ceiling(inv.data)}
        return(inv.data)
      },
      data.df.split.arm,
      lapply(seq_len(ncol(sim.mvnorm.data)), function(j)
        sim.mvnorm.data[, j]),
      mean.vec,
      diag(cov.mat)
    )
  }

#' Performing the hypothesis test to compare the difference between the empirical data and the simulated data
#'
#' @param x A numeric matrix.
#' @param y A numeric matrix which is compared to \code{x}.
#' @param test.method A string to specify the hypothesis test used to detect the difference
#'   between input data and the simulated data. Default is "none". Possible methods are
#'   energy distance ("energy") and  ball divergence ("ball"). The R packages "energy" and "Ball" are needed.
#' @return A list with two elements.
#'   1. p.value: the p-value of the hypothesis test.
#'   2. test.result: the returned object of the hypothesis test.
data.diff.test <- function(x, y, test.method) {
  if (test.method == "energy") {
    if(!requireNamespace("energy")) {stop("You need to install R package energy.")}
     test.res <- energy::eqdist.etest(rbind(x, y), c(nrow(x), nrow(y)), R=999)
    return(list(p.value = test.res$p.value, test.result = test.res))
  } else if (test.method == "ball") {
    if(!requireNamespace("Ball")) {stop("You need to install R package Ball.")}
     test.res <- Ball::bd.test(x, y)
    return(list(p.value = test.res$p.value, test.result = test.res))
  } else {stop("No such test.method!")}
}

#' To generate simulated datasets from empirical data by utilizing the copula invariance property.
#'
#' Based on the empirical data, generating simulated datasets through the copula invariance property.
#' @param data.input The empirical patient-level data to be used to simulate new virtual patient data.
#' @param id.vec The ID for individual patient in the input data.
#' @param arm.vec The column to identify the arm in clinical trial.
#' @param n.patient The targeted number of patients in each simulated dataset.
#' @param n.simulation The number of simulated datasets.
#' @param seed The random seed. Default is NULL to use the current seed.
#' @param validation.type A string to specify the hypothesis test used to detect the difference
#'   between input data and the simulated data. Default is "none". Possible methods are
#'   energy distance ("energy") and  ball divergence ("ball"). The R packages "energy" and "Ball" are needed.
#' @param validation.sig.lvl The significant level (alpha) value for the hypothesis test.
#' @param rmvnorm.matrix.decomp.method The method to do the matrix decomposition used in the function \code{rmvnorm}. Default is "svd".
#' @param verbose A logical value to specify whether to print message for simulation process or not.
#' @return A copula.sim object with four elements.
#'   1. data.input: empirical data (wide-form)
#'   2. data.input.long: empirical data (long-form)
#'   3. data.transform: quantile transformation of data.input
#'   4. data.simul: simulated data
#' @export
#' @importFrom magrittr set_colnames
#' @importFrom stats pnorm qnorm runif cov
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr between bind_rows group_split if_else left_join
#' @importFrom dplyr select ungroup arrange group_by summarise filter mutate n
#' @importFrom mvtnorm rmvnorm
#' @references Sklar, A. (1959). Functions de repartition an dimensionset leursmarges., Paris: PublInst Stat.
#' @references Nelsen, R. B. (2007). An introduction to copulas. Springer Science & Business Media.
#' @references Ross, S. M. (2013). Simulation. Academic Press.
#' @author Pei-Shan Yen, Xuemin Gu
#' @examples
##'
##' library(copulaSim)
##'
##' ## Generate Empirical Data
##'  # Assume the 2-arm, 5-dimensional empirical data follows multivariate normal data.
##' library(mvtnorm)
##' arm1 <- rmvnorm(n = 40, mean  = rep(10, 5), sigma = diag(5) + 0.5)
##' arm2 <- rmvnorm(n = 40, mean  = rep(12, 5), sigma = diag(5) + 0.5)
##' test_data <- as.data.frame(cbind(1:80, rep(1:2, each = 40), rbind(arm1, arm2)))
##' colnames(test_data) <- c("id","arm",paste0("time_", 1:5))
##'
##' ## Generate 100 simulated datasets
##' copula.sim(data.input = test_data[,-c(1,2)], id.vec = test_data$id, arm.vec = test_data$arm,
##' n.patient = 100 , n.simulation = 100, seed = 2022)

copula.sim <- function(data.input,
                       id.vec,
                       arm.vec,
                       n.patient,
                       n.simulation,
                       seed = NULL,
                       validation.type = "none",
                       validation.sig.lvl = 0.05,
                       rmvnorm.matrix.decomp.method = "svd",
                       verbose = TRUE) {
  # check dimensions of data.input/id.vec/arm.vec are all valid.
  if (length(id.vec) != length(arm.vec)) {
    stop("The length of id.vec should be equal to the length of arm.vec.")
  }
  if (is.null(dim(data.input))) {
    stop("data.input must be a numeric data.frame or matrix.")
  }
  if (length(id.vec) != nrow(data.input)) {
    stop("The length of id.vec should be equal to the number of rows of data.input.")
  }

  # check data input are valid
  data.input <- as.matrix(data.input)
  if (is.null(colnames(data.input))) {
    colnames(data.input) <- sprintf("Visit%i", seq.int(1L, nrow(data.input)))
  }
  if (all(data.input %in% c("integer", "numeric"))) {
    stop("data.input should be all numeric.") # either numerical or integer is acceptable
  }
  if (any(is.na(data.input))) {
    stop("data.input must be non-NA.")
  }

  # check id and arm are valid
  if (any(is.na(id.vec))) {
    stop("id.vec must be non-NA.")
  }
  if (any(is.na(arm.vec))) {
    stop("arm.vec must be non-NA.")
  }

  # check options are valid
  if ((abs(n.patient - floor(n.patient)) > 1e-6) || (n.patient <= 0)) {
    stop("n.patient must be a positive integer.")
  }
  if ((abs(n.simulation - floor(n.simulation)) > 1e-6) || (n.simulation <= 0)) {
    stop("n.simulation must be a positive integer.")
  }
  if (!is.null(seed) && !is.numeric(seed)) {
    stop("seed must be a number or NULL.")
  }
  if (!(validation.type %in% c("none", "energy", "ball", "equalCopula"))) {
    stop(
      "wrong validation.type, it must be one of ['none', 'energy', 'ball', 'equalCopula']."
    )
  }
  if ((length(validation.sig.lvl) > 1) || is.na(validation.sig.lvl) ||
       !is.numeric(validation.sig.lvl) || !between(validation.sig.lvl, 0, 1)) {
    stop("validation.sig.lvl must be non-NA numeric between 0 and 1.")
  }
  if (!(rmvnorm.matrix.decomp.method %in% c("eigen", "svd", "chol"))) {
    stop("rmvnorm.matrix.decomp.method must be one of ['eigen', 'svd', 'chol']")
  }
  if (!is.logical(verbose)) {
    stop("verbose must be logical.")
  }

  # Step 0: set random seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Step 1: Input empirical data
  colname.df <- tibble(
    col.num = seq.int(1L, ncol(data.input)),
    col.name = colnames(data.input)
  )
  n.arm <- length(unique(arm.vec))
  data.df <- cbind(
    rep(id.vec, ncol(data.input)),
    rep(arm.vec, ncol(data.input)),
    rep(seq.int(1L, ncol(data.input)), each = nrow(data.input))
  ) %>%
    set_colnames(c("id", "arm", "col.num")) %>%
    as_tibble %>%
    left_join(colname.df, by = "col.num")
  data.df$data.input <- as.vector(as.matrix(data.input))

  # Step 2: calculate the ranks of empirical marginal for performing quantile transformation.
  #         and mapping marginal standard uniform to marginal standard normal Z
  data.transform <- data.df %>%
    group_by(.data$arm, .data$col.num) %>%
    # check the data.input (per dimension) is integers
    mutate(data.is.integer = all(mean(.data$data.input - floor(.data$data.input)) <= 1e-6)) %>%
    # add shift values if data.input is integers
    mutate(data.input.transformed = .data$data.input + if_else(.data$data.is.integer, runif(n()) - 1, 0)) %>%
    # quantile transformation: marginal dist CDF^(-1) follow Unif(0,1) with reasonable bounds
    mutate(data.unif = (rank(.data$data.input) - 0.5) / max(unique(.data$id)),
           data.norm = qnorm((rank(.data$data.input) - 0.5) / max(unique(.data$id)))) %>%
    ungroup

  # get mean vectors and covariance matrices
  data.split.arm <- data.transform %>%
    group_by(.data$arm) %>%
    group_split

  data.norm.mean.list <- data.split.arm %>%
    lapply(function(df){
      df %>% group_by(.data$col.num) %>%
        summarise(norm.mean = mean(.data$data.norm)) %>%
        arrange(.data$col.num) %>% ungroup %>%
        select(.data$norm.mean) %>% unlist
    })

  data.norm.cov.mat.list <- data.split.arm %>%
    lapply(function(df){
      cov.mat <- df %>%
        arrange(.data$col.num, .data$id) %>%
        (function(temp.df) cov(matrix(temp.df$data.norm, ncol = max(temp.df$col.num))))
    })

  arm.unique.vec <- data.split.arm %>% sapply(function(df) df$arm[1])

  # get simulation data
  data.simul <- lapply(seq_len(n.simulation), function(sim.id){
    if (verbose)
      cat(sprintf("Simulate %ith Dataset\n", sim.id))
    st <- proc.time()
    sim.data.df <- mapply(function(arm, mean.vec, cov.mat) {
      con.sim <- TRUE
      while (con.sim) {
        # Step 3. Using multivariate mean and dependence structure to simulate multivariate normal distribution
        sim.data <- rmvnorm(n.patient,
                            mean = mean.vec,
                            sigma = cov.mat,
                            method = rmvnorm.matrix.decomp.method) %>%
        # Step 4. recover the simulated data to the range of empirical via quantile transformation
          empirical.ecdf.inv(arm, data.transform, mean.vec, cov.mat)
        if (validation.type == "none") {
          con.sim <- FALSE
        } else {
          # step 4-1. make a hypothesis test to compare the difference between emperical data and the simulated data
          p.value <- data.diff.test(data.input[arm.vec == arm, ], sim.data, validation.type)$p.value
          if(verbose)
            cat(sprintf("p.value for %s test: %.4f\n", validation.type, p.value))
          con.sim <- p.value < validation.sig.lvl
        }
      }
      # output the simulated data if no need for validation or pass the hypothesis test
      output.df <- cbind(
        1:n.patient,
        arm,
        rep(seq.int(1L, length(mean.vec)), each = n.patient)
      ) %>%
        set_colnames(c("id", "arm", "col.num")) %>%
        as_tibble %>%
        left_join(colname.df, by = "col.num")
      output.df$data.sim <- as.vector(sim.data)
      return(output.df)
    },
    arm.unique.vec,
    data.norm.mean.list,
    data.norm.cov.mat.list,
    SIMPLIFY = FALSE
    ) %>% bind_rows
    sim.data.df$sim.id <- sim.id
    if (verbose)
      cat(sprintf("Compelete simulating %ith Dataset in %.3f seconds\n", sim.id, (proc.time() - st)["elapsed"]))
    return(sim.data.df)
  }) %>% bind_rows

  # make a object
  res <- list(
    data.input = data.input,
    data.input.long = data.df,
    data.transform = data.transform,
    data.simul = data.simul
  )
  class(res) <- "copula.sim"
  return(res)
}

#' Converting data.simul in a copula.sim object into a list of wide-form matrices
#'
#' @param object A copula object.
#' @return A list of matrices for simulated data.
#' @export
#' @importFrom dplyr distinct
#' @importFrom magrittr set_colnames set_names
extract.data.sim <- function(object) {
  if (!("copula.sim" %in% class(object))) {
    stop("Only accept copula.sim object.")
  }
  matrix.colnames <- object$data.simul %>%
    distinct(.data$col.num, .data$col.name) %>%
    select(.data$col.name) %>% unlist %>% set_names(NULL)
  sim.data.split <- object$data.simul %>%
    group_by(.data$sim.id) %>%
    group_split
  output.list <- sim.data.split %>%
    lapply(function(df){
      data.split.arm <- df %>%
        group_by(.data$arm) %>%
        group_split
      outlist <- data.split.arm %>%
        lapply(function(subdf){
          subdf %>%
            arrange(.data$col.num, .data$id) %>%
            (function(temp.df) matrix(temp.df$data.sim, ncol = max(temp.df$col.num))) %>%
            set_colnames(matrix.colnames)
        })
      names(outlist) <- sapply(data.split.arm, function(df) sprintf("arm=%i", df$arm[1]))
      return(outlist)
    })
  names(output.list) <- sprintf("sim.id=%i", sapply(sim.data.split, function(df) df$sim.id[1]))
  return(output.list)
}
