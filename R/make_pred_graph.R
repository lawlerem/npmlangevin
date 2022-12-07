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

#' @export
pred_graph_to_cpp<- function(x) {
  x$v<- x$v - 1
  x$parents<- lapply(x$parents, `+`, -1)
  return(x)
}