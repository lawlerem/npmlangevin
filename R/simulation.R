#' Simulate a utilization distribution and track from a Langevin diffusion model
#'
#' @param xlim, ylim The approximate bounding box for movement
#' @param boundary_sharpness How penalized should movement outside the bounding box be?
#' @param boundary_limit A scaling factor applied to xlim and ylim to determine where the penalty be applied.
#' @param cv_code 0 = exponential (don't use), 1 = gaussian, 2 = matern, 3 = matern32
#' @param cv_pars c(marginal std. dev., range, smoothness)
#' @param pred_loc_delta Raster cell size for interpolated utilization distribution raster. If equal to NA, predictions won't be made.
#' @param nt How many true movement locations should there be?
#' @param nping How many location pings should there be?
#' @param gamma Speed parameter.
#' @param ping_tau Std. dev.s for ping observation error.
#' @param ping_cor Correlation between ping observation error coordinates.
#' @param loc_class_probs Named probability factor giving location quality class frequencies.
#' @param seed Optional simulation seed.
#'
#' @return A list:
#'   - nn_graph The output of make_nn_graph
#'   - field A stars object with the simulated utilization distribution
#'   - pred_field An upscaled / interpolated copy of field
#'   - track An sf data.frame with the true location track and field gradient values
#'   - pings An sf data.frame with the observation pings
#'   - tmap A list of tmap objects for the above simulated values
#'   - data A copy of data passed to TMB
#'   - para A copy of parameters passed to TMB
#'
#' @export
simulate<- function(
    xlim = c(-2, 2),
    ylim = c(-2, 2),
    boundary_sharpness = 2.5,
    boundary_limit = 0.8,
    cv_code = 1,
    cv_pars = c(2, 0.5, 6.5),
    pred_loc_delta = 0.1,
    nt = 100,
    nping = 40,
    gamma = 0.08,
    ping_tau = 0.01 * c(1, 1),
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
  if( is.na(pred_loc_delta) ) {
    pred_locs<- st_as_sf(
      expand.grid(
        x = numeric(1),
        y = numeric(1),
        v = 1:3
      ),
      coords = c("x", "y")
    )
    pwg<- list(
      v = numeric(0),
      coord = matrix(0, nrow = 0, ncol = 2),
      parents = list()
    )
  } else {
    pred_locs<- st_as_sf(
      expand.grid(
        x = seq(min(xlim), max(xlim), by = pred_loc_delta),
        y = seq(min(ylim), max(ylim), by = pred_loc_delta),
        v = 1:3
      ),
      coords = c("x", "y")
    )
    pwg<- make_pred_graph(
      pred_locs,
      g
    )
  }


  true_time<- sort(cumsum(runif(nt, 0, 2)))
  track_idx<- sort(sample(seq(nt), nping, replace = TRUE))
  loc_class<- sample(loc_class_K$q, nping, replace = TRUE, prob = loc_class_probs)

  data<- list(
    model = "langevin_diffusion",
    cv_code = cv_code,
    g = list(
      st_get_dimension_values(g$stars, "x"),
      st_get_dimension_values(g$stars, "y"),
      lapply(lapply(g$graph, `[[`, 1), `+`, -1),
      lapply(lapply(g$graph, `[[`, 2), `+`, -1)
    ),
    pwg = pred_graph_to_cpp(pwg),
    field_neighbours = lapply(true_time, function(x) {
      return(matrix(0, nrow = 1, ncol = 2))
    }),
    true_time = true_time,
    pings = list(
      coords = matrix(0, nrow = nping, ncol = 2),
      loc_class = as.numeric(loc_class) - 1,
      track_idx = track_idx - 1,
      K = as.matrix(loc_class_K[, c("x", "y")])
    )
  )
  para <- list(
    boundary_x = boundary_limit * xlim,
    boundary_y = boundary_limit * ylim,
    working_boundary_sharpness = log(boundary_sharpness),
    working_cv_pars = log(cv_pars),
    w = g$stars$w,
    true_coord = matrix(0, nrow = nt, ncol = 2),
    log_gamma = log(gamma),
    working_ping_cov_pars = c(
      log(ping_tau[[1]]),
      qlogis(0.5 + 0.5 * ping_cor),
      log(ping_tau[[2]])
    )
  )

  simobj<- MakeADFun(
    data = data,
    para = para,
    DLL = "npmlangevin_TMB"
  )
  sim<- simobj$simulate()
  g$stars$w<- sim$w

  if( !is.na(pred_loc_delta) ) {
    pred<- st_sf(
      data.frame(
        w = sim$pw,
        pred_locs
      )
    )
    pred<- split(pred, pred$v)
    pred<- lapply(
      pred,
      function(x) {
        x<- st_as_stars(x["w"])
        x<- st_sfc2xy(x)
        return( x )
      }
    )
    pred<- do.call(c, c(pred, list(along = "v")))
    st_dimensions(pred)$v$values<- c("gg", "dxdx", "dydy")
    names(st_dimensions(pred))[1:2]<- c("x", "y")
    attr(st_dimensions(pred), "raster")$dimensions<- c("x", "y")
  } else {
    pred<- g$stars
  }

  track<- st_as_sf(
    data.frame(
      x = sim$true_coord,
      dx = sim$track_gradient[, 1],
      dy = sim$track_gradient[, 2],
      t = true_time
    ),
    coords = c(1, 2)
  )
  pings<- st_as_sf(
    data.frame(
      x = sim$sim_pings,
      t = true_time[track_idx],
      q = loc_class
    ),
    coords = c(1, 2)
  )


  field_tm<- tm_shape(g$stars["w"]) +
    tm_raster(
      style = "cont",
      midpoint = 0,
      interpolate = FALSE,
      # palette = "viridis"
      palette = "PRGn"
    ) +
    tm_facets(nrow = 1, ncol = 3, free.scales = FALSE) +
    tm_layout(legend.outside = FALSE)

  pred_field_tm<- tm_shape(pred["w"]) +
    tm_raster(
      style = "cont",
      midpoint = 0,
      interpolate = TRUE,
      # palette = "viridis"
      palette = "PRGn"
    ) +
    tm_facets(nrow = 1, ncol = 3, free.scales = FALSE) +
    tm_layout(legend.outside = FALSE)

  sim_util<- exp(g$stars["w", , , "gg"])
  sim_util<- sim_util / sum(sim_util[["w"]])
  util_tm<- tm_shape(sim_util) +
    tm_raster(
      style = "cont",
      interpolate = FALSE,
      palette = "Greens"
    )

  pred_util<- exp(pred["w", , , "gg"])
  pred_util<- pred_util / sum(pred_util[["w"]])
  pred_util_tm<- tm_shape(pred_util) +
    tm_raster(
      style = "cont",
      interpolate = TRUE,
      palette = "Greens"
    )

  track_tm<- tm_shape(st_cast(st_combine(track), "LINESTRING")) + tm_lines(col = "black")
  pings_tm<- tm_shape(pings) + tm_dots(col = "q", palette = "YlOrBr", size = 0.2)


  return(
    list(
      nn_graph = g,
      field = g$stars,
      pred_field = pred,
      pred_locs = pred_locs,
      track = track,
      pings = pings,
      tmap = list(
        field = field_tm,
        pred_field = pred_field_tm,
        util = util_tm,
        pred_util = pred_util_tm,
        track = track_tm,
        pings = pings_tm
      ),
      data = data,
      para = para
    )
  )
}