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
#' @param filtered_locations The output of fit_rw, e.g. a list
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
#' @param max.edge The maximum edge length for the mesh
#' @param ... Additional arguments to pass to make_starve_graph
#'
#' @return A list with the following elements:
#'   - filtered_locations: A copy of the filtered_locations argument
#'   - graph: The mesh used for the random field
#'   - track_graph: The parents used to connect the track and field
#'   - opt: The output of nlminb
#'   - sdr: The output of sdreport
#'   - cv_code: The covariance function code
#'   - mesh_predictions: A data.frame containing the field predictions for the mesh.
#' 
#' @export
fit_utilization_distribution<- function(
    filtered_locations,
    cv_code = 1,
    max.edge = 1,
    ...
  ) {
  pings<- filtered_locations$pings
  pings$dx<- c(
    tail(sf::st_coordinates(pings)[, 1], -1) - head(sf::st_coordinates(pings)[, 1], -1),
    NA
  )
  pings$dy<- c(
    tail(sf::st_coordinates(pings)[, 2], -1) - head(sf::st_coordinates(pings)[, 2], -1),
    NA
  )
  filtered_locations$pings<- pings

  graph<- make_starve_graph(
    filtered_locations$track,
    max.edge = max.edge,
    ...
  )
  track_graph<- make_starve_pred_graph(
    pred_coordinates = filtered_locations$track,
    field_coordinates = graph$coordinates
  )
  
  data<- list(
    model = "starve_npmlangevin",
    cv_code = cv_code,
    g = list(
      sf::st_coordinates(graph$coordinates),
      lapply(lapply(graph$edge_list, `[[`, 1), `+`, -1),
      lapply(lapply(graph$edge_list, `[[`, 2), `+`, -1)
    ),
    pwg = list(
      coord = matrix(0, nrow = 0, ncol = 2),
      parents = list()
    ),
    coordinates = sf::st_coordinates(filtered_locations$track),
    field_neighbours = lapply(track_graph$parents, `+`, -1),
    time = filtered_locations$track$t,
    location_differences = as.matrix(
      head(
        filtered_locations$pings[, c("dx", "dy"), drop = TRUE],
        -1
      )
    ),
    location_quality_class = as.numeric(filtered_locations$pings$q) - 1,
    K = as.matrix(loc_class_K[, c("x", "y")])
  )
  para<- list(
    working_cv_pars = log(c(1, 1)),
    w = matrix(0, nrow = nrow(graph$coordinates), ncol = 2),
    random_walk = matrix(0, nrow = nrow(filtered_locations$track) - 1, ncol = 2),
    log_gamma = 0.1 * filtered_locations$parameters[["log_gamma"]],
    working_ping_cov_pars = filtered_locations$parameters[
      names(filtered_locations$parameters) %in% c("working_obs_cov_pars")
    ]
  )
  obj<- TMB::MakeADFun(
    data = data,
    para = para,
    random = c("w", "random_walk"),
    DLL = "npmlangevin_TMB"
  )
  opt<- nlminb(
    obj$par,
    obj$fn,
    obj$gr
  )
  sdr<- TMB::sdreport(
    obj,
    opt$par
  )

  pwg<- make_starve_gg_pred_graph(
    pred_coordinates = graph$coordinates,
    field_coordinates = graph$coordinates,
    cv_pars = as.list(sdr, "Est", report = TRUE)$cv_pars,
    cv_code = cv_code,
    k = 1
  )
  data$pwg[[1]]<- sf::st_coordinates(pwg$coordinates)
  data$pwg[[2]]<- lapply(pwg$parents, `+`, -1)
  obj<- TMB::MakeADFun(
    data = data,
    para = para,
    random = c("w", "random_walk"),
    DLL = "npmlangevin_TMB"
  )
  obj$fn(opt$par)
  sdr<- sdreport(
    obj,
    opt$par
  )
  mesh_predictions<- sf::st_as_sf(
    data.frame(
      as.list(sdr, "Est", report = TRUE)$pw,
      as.list(sdr, "Est")$w,
      as.list(sdr, "Std", report = TRUE)$pw,
      as.list(sdr, "Std")$w,
      graph$coordinates
    )
  )
  colnames(mesh_predictions)[1:6]<- c("g", "dx", "dy", "g_se", "dx_se", "dy_se")

  return(
    list(
      filtered_locations = filtered_locations,
      graph = graph,
      track_graph = track_graph,
      opt = opt,
      sdr = sdr,
      cv_code = cv_code,
      mesh_predictions = mesh_predictions
    )
  )
}


#' Use a fitted langevin diffusion model to predict the utilization distribution
#' 
#' @param prediction_locations An sf object with point geometries
#' @param fitted_model The output of fit_utilization_distribution
#' @param k The number of parents in each direction
#' 
#' @return An sf object with predictions for the log-utilization distribution
#' 
#' @export
predict_utilization_distribution<- function(
    prediction_locations,
    fitted_model,
    k = 1,
    ...
  ) {
  fm<- fitted_model
  pwg<- make_starve_gg_pred_graph(
    pred_coordinates = prediction_locations,
    field_coordinates = fm$graph$coordinates,
    cv_pars = as.list(fm$sdr, "Est", report = TRUE)$cv_pars,
    cv_code = fm$cv_code,
    k = k
  )
  data<- list(
    model = "starve_npmlangevin",
    cv_code = fm$cv_code,
    g = list(
      sf::st_coordinates(fm$graph$coordinates),
      lapply(lapply(fm$graph$edge_list, `[[`, 1), `+`, -1),
      lapply(lapply(fm$graph$edge_list, `[[`, 2), `+`, -1)
    ),
    pwg = list(
      sf::st_coordinates(pwg$coordinates),
      lapply(pwg$parents, `+`, -1)
    ),
    coordinates = sf::st_coordinates(fm$filtered_locations$track),
    field_neighbours = lapply(fm$track_graph$parents, `+`, -1),
    time = fm$filtered_locations$track$t,
    location_differences = as.matrix(
      head(
        fm$filtered_locations$pings[, c("dx", "dy"), drop = TRUE],
        -1
      )
    ),
    location_quality_class = as.numeric(fm$filtered_locations$pings$q) - 1,
    K = as.matrix(loc_class_K[, c("x", "y")])
  )
  para<- list(
    working_cv_pars = c(0, 0),
    w = matrix(0, nrow = nrow(fm$graph$coordinates), ncol = 2),
    random_walk = matrix(0, nrow = nrow(fm$filtered_locations$track) - 1, ncol = 2),
    log_gamma = 0.1 * fm$filtered_locations$parameters[["log_gamma"]],
    working_ping_cov_pars = filtered_locations$parameters[
      names(fm$filtered_locations$parameters) %in% c("working_obs_cov_pars")
    ]
  )
  obj<- TMB::MakeADFun(
    data = data,
    para = para,
    random = c("w", "random_walk"),
    DLL = "npmlangevin_TMB"
  )
  obj$fn(fm$opt$par)
  sdr<- sdreport(
    obj,
    fm$opt$par
  )
  predictions<- sf::st_as_sf(
    data.frame(
      g = as.list(sdr, "Est", report = TRUE)$pw,
      g_se = as.list(sdr, "Std", report = TRUE)$pw,
      pwg$coordinates
    )
  )
  return(predictions)
}

