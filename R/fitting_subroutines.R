#' Pre-filter a track using a random walk model
#'
#' @param locations An n x 2 sf data.frame (t, q, geom) giving the observed projected coordinates (point geometries), observation time, and location quality class of a movement path. The time column should be either a POSIXt column or a numeric column.
#' @param delta_t Time between locations in estimated true path, see ?seq.POSIXt
#'
#' @return A list
#'   - pings A time-sorted copy of the passed in locations
#'   - track An sf data.frame with columns
#'     - x_se, y_se Standard errors for the x and y coordinate, respectively
#'     - t Time for the location estimate
#'     - geom Point geometries giving the estimated locations
#'   - parameters parameter estimates
#'
#' @export
fit_rw<- function(locations, delta_t = NA) {
  locations<- locations[order(locations$t), , drop = FALSE]
  if( is.na(delta_t) ) {
    regular_t<- NULL
  } else {
    regular_t<- seq(min(locations$t), max(locations$t), by = delta_t)
  }
  true_time<- unique(c(locations$t, regular_t))
  true_time<- sort(true_time)

  data<- list(
    model = "random_walk",
    true_time = unname(true_time),
    pings = list(
      coords = unname(sf::st_coordinates(locations)),
      loc_class = as.numeric(locations$q) - 1,
      track_idx = match(locations$t, true_time) - 1,
      K = as.matrix(loc_class_K[, c("x", "y")])
    )
  )
  para<- list(
    true_loc = matrix(0, nrow = length(true_time), ncol = 2),
    log_gamma = 0,
    working_obs_cov_pars = numeric(3)
  )

  obj<- TMB::MakeADFun(
    data = data,
    para = para,
    random = "true_loc",
    DLL = "npmlangevin_TMB"
  )
  opt<- nlminb(obj$par, obj$fn, obj$gr)
  sdr<- sdreport(obj, opt$par)

  return(
    list(
      pings = locations,
      track = sf::st_set_crs(
        sf::st_as_sf(
          data.frame(
            x = as.list(sdr, "Est")$true_loc[, 1],
            x_se = as.list(sdr, "Std")$true_loc[, 1],
            y = as.list(sdr, "Est")$true_loc[, 2],
            y_se = as.list(sdr, "Std")$true_loc[, 2],
            t = true_time
          ),
          coords = c("x", "y")
        ),
        sf::st_crs(locations)
      ),
      parameters = opt$par
    )
  )
}

#' Fit a utilization distribution using a filtered movement track estimate
#'
#' @param track_estimate The output of fit_rw or fit_langevin, e.g. a list
#'   - pings An sf data.frame with columns
#'     - t Time for the location ping
#'     - q Location quality class for the location ping
#'     - geom Point geometries giving the observed location pings
#'   - track An sf data.frame with columns
#'     - x_se, y_se Standard errors for the x and y coordinate, respectively
#'     - t Time for the location estimate
#'     - geom Point geometries giving the estimated locations
#'   - parameters parameter estimates
#' @param cv_code 0 = exponential (don't use), 1 = gaussian, 2 = matern, 3 = matern32
#' @param The initial covariance function parameters used to build the nearest neighbour graph. Can either be a length 3 vector (std. dev., range, smoothness) or a k x 3 matrix. If a matrix is supplied, each row will be tested to find which gives the highest likelihood value.
#' @param bbox A named vector with elements "xmin", "xmax", "ymin", and "ymax"
fit_utilization_distribution<- function(
    track_estimate,
    cv_code = 1,
    initial_cv_pars = c(1, 1, 1.5),
    bbox = 1.1 * sf::st_bbox(track_estimate$pings)
  ) {
  # Create nearest neighbour graph
  g<- make_nn_graph(
    x = bbox[c("xmin", "xmax")],
    y = bbox[c("ymin", "ymax")],
    cv_pars = initial_cv_pars,
    cv_code = cv_code
  )
  track_nn<- find_nearest_four(
    track_estimate$track,
    g
  )

  data<- list(
    model = "langevin_diffusion",
    cv_code = cv_code,
    g = list(
      stars::st_get_dimension_values(g$stars, "x"),
      stars::st_get_dimension_values(g$stars, "y"),
      lapply(lapply(g$graph, `[[`, 1), `+`, -1),
      lapply(lapply(g$graph, `[[`, 2), `+`, -1)
    ),
    pwg = list(
      var = integer(0),
      coord = matrix(0, nrow = 0, ncol = 2),
      parents = list()
    ),
    field_neighbours = lapply(track_nn, `+`, -1),
    true_time = track_estimate$track$t,
    pings = list(
      coords = unname(sf::st_coordinates(track_estimate$pings)),
      loc_class = as.numeric(track_estimate$pings$q) - 1,
      track_idx = match(track_estimate$pings$t, track_estimate$track$t) - 1,
      K = as.matrix(loc_class_K[, c("x", "y")])
    )
  )
  para<- list(
    boundary_x = bbox[c("xmin", "xmax")],
    boundary_y = bbox[c("ymin", "ymax")],
    working_boundary_sharpness = log(0),
    working_cv_pars = log(initial_cv_pars),
    w = g$stars$w,
    true_coord = sf::st_coordinates(track_estimate$track),
    log_gamma = track_estimate$parameters["log_gamma"],
    working_ping_cov_pars = track_estimate$parameters[names(track_estimate$parameters) == "working_obs_cov_pars"]
  )
  map<- list(
    boundary_x = as.factor(
      c(NA, NA)
    ),
    boundary_y = as.factor(
      c(NA, NA)
    ),
    working_boundary_sharpness = as.factor(
      NA
    ),
    working_cv_pars = as.factor(
      c(1, 2, NA)
    ),
    true_coord = as.factor(
      matrix(NA, nrow = nrow(para$true_coord), ncol = ncol(para$true_coord))
    ),
    log_gamma = as.factor(NA),
    working_ping_cov_pars  = as.factor(
      rep(NA, length(para$working_ping_cov_pars))
    )
  )
  obj<- TMB::MakeADFun(
    data = data,
    para = para,
    map = map,
    random = c("w", "true_coord"),
    DLL = "npmlangevin_TMB"
  )
  opt<- nlminb(obj$par, obj$fn, obj$gr)
  sdr<- TMB::sdreport(obj, opt$par)


  ###
  ### Refit with updated graph
  ###
  updated_cv_pars<- c(exp(opt$par), 1.5)
  g<- make_nn_graph(
    x = bbox[c("xmin", "xmax")],
    y = bbox[c("ymin", "ymax")],
    cv_pars = exp(opt$par),
    cv_code = cv_code
  )
  track_nn<- find_nearest_four(
    track_estimate$track,
    g
  )

  data<- list(
    model = "langevin_diffusion",
    cv_code = cv_code,
    g = list(
      stars::st_get_dimension_values(g$stars, "x"),
      stars::st_get_dimension_values(g$stars, "y"),
      lapply(lapply(g$graph, `[[`, 1), `+`, -1),
      lapply(lapply(g$graph, `[[`, 2), `+`, -1)
    ),
    pwg = list(
      var = integer(0),
      coord = matrix(0, nrow = 0, ncol = 2),
      parents = list()
    ),
    field_neighbours = lapply(track_nn, `+`, -1),
    true_time = track_estimate$track$t,
    pings = list(
      coords = unname(sf::st_coordinates(track_estimate$pings)),
      loc_class = as.numeric(track_estimate$pings$q) - 1,
      track_idx = match(track_estimate$pings$t, track_estimate$track$t) - 1,
      K = as.matrix(loc_class_K[, c("x", "y")])
    )
  )
  para<- list(
    boundary_x = bbox[c("xmin", "xmax")],
    boundary_y = bbox[c("ymin", "ymax")],
    working_boundary_sharpness = log(1),
    working_cv_pars = log(initial_cv_pars),
    w = g$stars$w,
    pw = numeric(0),
    true_coord = unname(sf::st_coordinates(track_estimate$track)),
    log_gamma = track_estimate$parameters["log_gamma"],
    working_ping_cov_pars = track_estimate$parameters[names(track_estimate$parameters) == "working_obs_cov_pars"]
  )
  map<- list(
    boundary_x = as.factor(
      c(NA, NA)
    ),
    boundary_y = as.factor(
      c(NA, NA)
    ),
    working_boundary_sharpness = as.factor(
      NA
    ),
    working_cv_pars = as.factor(
      c(1, 2, NA)
    ),
    true_coord = as.factor(
      matrix(NA, nrow = nrow(para$true_coord), ncol = ncol(para$true_coord))
    ),
    log_gamma = as.factor(NA),
    working_ping_cov_pars  = as.factor(
      rep(NA, length(para$working_ping_cov_pars))
    )
  )
  obj<- TMB::MakeADFun(
    data = data,
    para = para,
    map = map,
    random = c("w", "pw", "true_coord"),
    DLL = "npmlangevin_TMB"
  )
  opt<- nlminb(obj$par, obj$fn, obj$gr)

  #
  sdr<- TMB::sdreport(obj, opt$par)
  sdr_est<- as.list(sdr, "Est")
  sdr_se<- as.list(sdr, "Std")
}
