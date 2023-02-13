#' Create a prediction graph
#'
#' @param pred_coordinates An sf data.frame with point geometries. The
#'   locations for which predictions will be made.
#' @param field_coordinates An sf data.frame with point geometries. The
#'   locations used in the mesh.
#'
#' @return A list with the prediction coordinates and list of parents.
#'
#' @export
make_starve_pred_graph<- function(pred_coordinates, field_coordinates) {
    nn<- st_nn(
        pred_coordinates,
        field_coordinates,
        sparse = TRUE,
        k = 3,
        returnDist = FALSE
    )
    return(
        list(
            coordinates = pred_coordinates,
            parents = nn
        )
    )
}