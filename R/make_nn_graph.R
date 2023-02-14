#' Make a nearest neighbour grid optimized for a covariance function
#'
#' @param xlim, ylim Vectors of length 2 giving the desired bounding box.
#' @param cv_pars Vector of length 3 - [marginal std. dev., range, smoothness]
#' @param cv_code 0 = exponential, 1 = gaussian, 2 = matern, 3 = matern (nu = 1.5)
#'
#' @return A named list
#'   - stars: A stars object with the raster locations, variables, and values (w and se)
#'   - graph: A directed acyclic graph
#'
#' @export
make_nn_graph<- function(xlim, ylim, cv_pars = c(1, 0.3, 2.5), cv_code = 1) {
  gr_obj<- MakeADFun(
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

  x<- seq(xlim[[1]], xlim[[2]], by = inflection_point)
  xl<- length(x)
  y<- seq(ylim[[1]], ylim[[2]], by = inflection_point)
  yl<- length(y)
  v<- factor(c("gg", "dxdx", "dydy"))
  vl<- length(v)

  a<- st_as_stars(
    list(
      w = array(0, dim = c(xl, yl, vl)),
      se = array(0, dim = c(xl, yl, vl))
    ),
    dimensions = st_dimensions(
      x = x,
      y = y,
      v = v,
      .raster = c("x", "y")
    )
  )

  parent_finder<- function(index) {
    i<- index[[1]]
    j<- index[[2]]
    k<- index[[3]]
    if( k == 1 ) {
      # gg
      parents<- rbind(
        c(i - 2, j, 2),
        c(i - 1, j, 2),
        c(i + 1, j, 2),
        c(i + 2, j, 2),
        c(i, j - 2, 3),
        c(i, j - 1, 3),
        c(i, j + 1, 3),
        c(i, j + 2, 3)
      )
    } else if( k == 2 ) {
      # dx_dx
      parents<- rbind(
        c(i - 2, j, 2),
        c(i, j - 1, 2),
        c(i, j - 2, 2),
        c(i - 1, j - 1, 3),
        c(i - 1, j + 1, 3)
      )
    } else {
      # dy_dy
      parents<- rbind(
        c(i, j - 2, 3),
        c(i - 1, j, 3),
        c(i - 2, j, 3),
        c(i - 1, j - 1, 2),
        c(i - 1, j + 1, 2)
      )
    }
    colnames(parents)<- c("i", "j", "k")
    parents<- parents[
      1 <= parents[, 1] & parents[, 1] <= xl &
      1 <= parents[, 2] & parents[, 2] <= yl,
      ,
      drop = FALSE
    ]
    return( parents )
  }

  # Need to put graph in the right order so that simulations are done correctly.
  # If [i, j, k] is a parent of [i', j', k'] then [i, j, k] needs to come before
  # [i', j', k'] in the graph
  idx<- expand.grid(
    i = seq(xl),
    j = seq(yl),
    k = 2:3
  )
  idx<- idx[with(idx, order(i, j, k)), ]
  idx<- rbind(
    idx,
    expand.grid(
      i = seq(xl),
      j = seq(yl),
      k = 1
    )
  )
  idx<- as.matrix(idx)
  rownames(idx)<- NULL

  idx_order<- array(0, dim = c(xl, yl, vl))
  for( row in seq(nrow(idx)) ) {
    idx_order[idx[row, "i"], idx[row, "j"], idx[row, "k"]]<- row
  }

  graph_ordered<- TRUE
  nn_list<- lapply(seq(nrow(idx)), function(i) {
    node<- list(
      to = idx[i, , drop = FALSE],
      from = parent_finder(idx[i, ])
    )
    if( !all(idx_order[node$from] < idx_order[node$to]) ) {
      graph_ordered<<- FALSE
    }

    return( node )
  })

  return( list(stars = a, graph = nn_list, graph_ordered = graph_ordered) )
}