#' Make a nearest neighbour graph for a set of prediction locations
#'
#' @param pred_coordinates The prediction locations as an sf data.frame with point geometries, with an additional column named "v" for the variable (1 = original, 2 = dx, 3 = dy).
#' @param nn_graph The output of \code{make_nn_graph}.
#'
#' @return A named list
#'   - v An integer vector specifying the field to be predicted.
#'   - coord A matrix of size n x 2 giving the prediction locations.
#'   - parents A list of integer matrices giving the index of locations used for the prediction.
#'
#' @export
make_pred_graph<- function(pred_coordinates, nn_graph) {
  nn_sf<- st_geometry(st_as_sf(nn_graph[[1]], as_points = TRUE))
  nn_cells<- st_nn(pred_coordinates, nn_sf, sparse = TRUE, k = 4, returnDist = FALSE)
  nn_raster<- as(nn_graph[[1]]["w", , , 1, drop = TRUE], "Raster")
  parents<- lapply(seq_along(nn_cells), function(i) {
    ans<- rowColFromCell(nn_raster, nn_cells[[i]])[, 2:1]
    ans<- cbind(ans, pred_coordinates$v[[i]])
    colnames(ans)<- c("i", "j", "k")
    return( ans )
  })

  return(list(
    v = pred_coordinates$v,
    coord = sf::st_coordinates(pred_coordinates),
    parents = parents
  ))
}

#' Find the four nearest neighbours of the spatial field graph
#'
#' @param pred_coordinates An sf object with point geometries
#' @param nn_graph The output of make_nn_graph
#'
#' @return A list of matrices giving the neighbour indices for each prediction location
#'
#' @export
find_nearest_four<- function(pred_coordinates, nn_graph) {
  pred_coordinates<- sf::st_coordinates(pred_coordinates)
  nearest_four<- lapply(seq(nrow(pred_coordinates)), function(i) {
    xcoord_d<- abs(stars::st_get_dimension_values(nn_graph$stars, "x") - pred_coordinates[i, 1])
    ycoord_d<- abs(stars::st_get_dimension_values(nn_graph$stars, "x") - pred_coordinates[i, 1])

    nn<- as.matrix(
      expand.grid(
        order(xcoord_d)[1:2],
        order(ycoord_d)[1:2]
      )
    )
    return(nn)
  })
  return(nearest_four)
}

#' @export
pred_graph_to_cpp<- function(x) {
  x$v<- x$v - 1
  x$parents<- lapply(x$parents, `+`, -1)
  return(x)
}