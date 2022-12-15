#' Simulate a track from a Langevin diffusion model
#'
#' @param xlim, ylim The approximate bounding box for movement
#' @param boundary_sharpness How penalized should movement outside the bounding box be?
#' @param boundary_limit A scaling factor applied to xlim and ylim to determine where the penalty be applied.
#' @param cv_code 0 = exponential (don't use), 1 = gaussian, 2 = matern, 3 = matern32
#' @param cv_pars c(marginal std. dev., range, smoothness)
#' @param nt How many true movement locations should there be?
#' @param nping How many location pings should there be?
#' @param gamma Speed parameter.
#' @param ping_tau Std. dev.s for ping observation error.
#' @param ping_cor Correlation between ping observation error coordinates.
#' @param loc_class_probs Named probability factor giving location quality class frequencies.
#' @param seed Optional simulation seed.
#'
#' @return A list:
#'   - true An sf data.frame with the true location track
#'   - pings An sf data.frame with the observation pings
#'
#' @export
simulate_track<- function(
    xlim = c(-2, 2),
    ylim = c(-2, 2),
    boundary_sharpness = 2.5,
    boundary_limit = 0.8,
    cv_code = 1,
    cv_pars = c(2, 0.5, 6.5),
    nt = 100,
    nping = 40,
    gamma = 0.08,
    ping_tau = 0.1 * c(1, 1),
    ping_cor = 0.3,
    loc_class_probs = c(
      "G" = 0.026,
      "3" = 0.040,
      "2" = 0.035,
      "1" = 0.020,
      "0" = 0.065,
      "A" = 0.126,
      "B" = 0.688
    ),
    seed
  ) {
  if( !missing(seed) ) {
    set.seed(seed)
  } else {}

  g<- make_nn_graph(
    x = xlim,
    y = ylim,
    cv_pars = cv_pars,
    cv_code = cv_code
  )

  true_time<- sort(cumsum(runif(nt, 0, 2)))
  track_idx<- sort(sample(seq(nt), nping, replace = TRUE))
  loc_class<- sample(loc_class_K$q, nping, replace = TRUE, prob = loc_class_probs)

  simobj<- MakeADFun(
    data = list(
      model = "langevin_diffusion",
      cv_code = cv_code,
      g = list(
        st_get_dimension_values(g$stars, "x"),
        st_get_dimension_values(g$stars, "y"),
        lapply(lapply(g$graph, `[[`, 1), `+`, -1),
        lapply(lapply(g$graph, `[[`, 2), `+`, -1)
      ),
      pwg = list(
        var = numeric(0),
        coord = matrix(0, nrow = 0, ncol = 2),
        parents = list()
      ),
      true_time = true_time,
      pings = list(
        coords = matrix(0, nrow = nping, ncol = 2),
        loc_class = as.numeric(loc_class) - 1,
        track_idx = track_idx - 1,
        K = as.matrix(loc_class_K[, c("x", "y")])
      )
    ),
    para = list(
      boundary_x = boundary_limit * xlim,
      boundary_y = boundary_limit * ylim,
      working_boundary_sharpness = log(boundary_sharpness),
      working_cv_pars = log(cv_pars),
      w = g$stars$w,
      pw = numeric(0),
      true_coord = matrix(0, nrow = nt, ncol = 2),
      log_gamma = log(gamma),
      working_ping_cov_pars = c(
        log(ping_tau[[1]]),
        qlogis(0.5 + 0.5 * ping_cor),
        log(ping_tau[[2]])
      )
    ),
    map = list(
      boundary_x = as.factor(c(NA, NA)),
      boundary_y = as.factor(c(NA, NA)),
      working_boundary_sharpness = as.factor(NA)
    ),
    random = c("w", "pw", "true_coord"),
    DLL = "npmlangevin_TMB"
  )
  sim<- simobj$simulate()

  return(
    list(
      true = st_as_sf(
        data.frame(
          x = sim$true_coord,
          t = true_time
        ),
        coords = c(1, 2)
      ),
      pings = st_as_sf(
        data.frame(
          x = sim$sim_pings,
          t = true_time[track_idx],
          q = loc_class
        ),
        coords = c(1, 2)
      )
    )
  )
}