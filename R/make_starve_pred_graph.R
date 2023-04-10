#' Create a prediction graph
#'
#' @param pred_coordinates An sf data.frame with point geometries. The
#'   locations for which predictions will be made.
#' @param field_coordinates An sf data.frame with point geometries. The
#'   locations used in the mesh.
#' @param k The number of parents to use.
#'
#' @return A list with the prediction coordinates and list of parents.
#'
#' @export
make_starve_pred_graph<- function(pred_coordinates, field_coordinates, k = 3) {
    nn<- nngeo::st_nn(
        pred_coordinates,
        field_coordinates,
        sparse = TRUE,
        k = k,
        returnDist = FALSE
    )
    return(
        list(
            coordinates = pred_coordinates,
            parents = nn
        )
    )
}


#' Create a prediction graph used to predict the (log) utilization distribution.
#'
#' @param pred_coordinates An sf data.frame with point geometries. The
#'   locations for which predictions will be made.
#' @param field_coordinates An sf data.frame with point geometries. The
#'   locations used in the mesh.
#' @param cv_pars Vector of length 2 or 3 - [std. dev., range, smoothness]
#' @param cv_code 0 = exponential, 1 = gaussian, 2 = matern, 3 = matern (nu = 1.5)
#' @param k Number of parents in each direction
#'
#' @return A list with the prediction coordinates and list of parents.
#'
#' @export
make_starve_gg_pred_graph<- function(
      pred_coordinates,
      field_coordinates,
      cv_pars,
      cv_code = 1,
      k = 1
    ) {
  # 1.) For each pred_coordinate, find the field_coordinate to the (left/top/right/bottom)
  #   that are closest to being distance r away, where r is inflection point of covariance function
  #   For boundaries, take the closest point to being distance r away that doesn't cross a boundary.
  shift_neighbour<- function(shift, pred_coordinates, field_coordinates) {
    crs<- sf::st_crs(pred_coordinates)
    pred_coordinates<- sf::st_coordinates(pred_coordinates)
    pred_coordinates<- data.frame(t(t(pred_coordinates) + shift))
    pred_coordinates<- sf::st_as_sf(pred_coordinates, coords = c(1, 2))
    sf::st_crs(pred_coordinates)<- crs
    nn<- nngeo::st_nn(
      pred_coordinates,
      field_coordinates,
      sparse = TRUE,
      k = k,
      returnDist = FALSE
    )

    return(do.call(c, nn))
  }
  gr_obj<- TMB::MakeADFun(
    data = list(
      model = "covariance_1d_deriv",
      cv_code = cv_code,
      cv_pars = cv_pars
    ),
    para = list(
      x = 0
    ),
    DLL = "npmlangevin_TMB"
  )
  gr_opt<- nlminb(gr_obj$par, gr_obj$fn, gr_obj$gr)
  inflection_point<- abs(c(gr_opt$par))
  shifts<- as.list(
    as.data.frame(
      rbind(
        c(-1, 1, 0, 0) * inflection_point,
        c(0, 0, -1, 1) * inflection_point
      )
    )
  )
  parents<- lapply(
    shifts,
    shift_neighbour,
    pred_coordinates = pred_coordinates,
    field_coordinates = field_coordinates
  )
  parents<- as.list(
    as.data.frame(
      do.call(
        rbind,
        parents
      )
    )
  )
  parents<- lapply(
    parents,
    function(x) {
      return(
        cbind(
          i = x,
          v = rep(c(2, 3), each = 2 * k)
        )
      )
    }
  )
  names(parents)<- NULL
  return(
    list(
      coordinates = pred_coordinates,
      parents = parents
    )
  )
}