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