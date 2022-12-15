#' Pre-filter a track using a random walk model
#'
#' @param locations An n x 2 sf data.frame (t, q, geom) giving the observed projected coordinates (point geometries), observation time, and location quality class of a movement path. The time column should be either a POSIXt column or a numeric column.
#' @param delta_t Time between locations in estimated true path, see ?seq.POSIXt
#'
#' @return A list
#'   - Locations A time-sorted copy of the passed in locations
#'   - Estimated An sf data.frame with columns
#'     - x_se, y_se Standard errors for the x and y coordinate, respectively
#'     - t Time for the location estimate
#'     - geom Point geometries giving the estimated locations
#'   - Parameters parameters estimates
#'
#' @export
fit_rw<- function(locations, delta_t) {
  locations<- locations[order(locations$t), , drop = FALSE]
  regular_t<- seq(min(locations$t), max(locations$t), by = delta_t)
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
      Locations = locations,
      Estimated = sf::st_set_crs(
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
      Parameters = opt$par
    )
  )
}

fit_utilization_distribution<- function() {

}
