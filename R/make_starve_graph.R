#' Make a starve-style graph using an INLA mesh.
#'
#' @param x An sf object with point locations
#' @param max.edge The largest allowed triangle edge length. See INLA::inla.mesh.2d.
#' @param ... Additional options to pass to INLA::inla.mesh.2d.
#'
#' @return A list with the mesh, node locations, and graph.
#'
#' @export
make_starve_graph<- function(
    x,
    max.edge = 1,
    ...
    ) {
    mesh<- INLA::inla.mesh.2d(
        loc = sf::st_coordinates(x),
        max.edge = max.edge,
        ...
    )
    mesh_nodes<- sf::st_as_sf(
        as.data.frame(mesh$loc[, c(1, 2)]),
        coords = c(1, 2)
    )
    mesh_graph<- as.matrix(mesh$graph$vv)
    o<- order_adjacency_matrix(mesh_graph) + 1
    mesh_nodes<- mesh_nodes[o, ]
    mesh_graph<- mesh_graph[o, o]
    mesh_graph[lower.tri(mesh_graph)]<- 0

    n_init<- 1
    edge_list<- vector(
        mode = "list",
        length = nrow(mesh_nodes) - n_init + 1
    )
    edge_list[[1]]<- list(
        to = seq(n_init),
        from = numeric(0)
    )
    edge_list[2:length(edge_list)]<- lapply(
        (n_init + 1):nrow(mesh_nodes),
        function(i) {
            return(
                list(
                    to = i,
                    from = which(mesh_graph[, i] != 0)
                )
            )
        }
    )

    return(
        list(
            mesh = mesh,
            coordinates = mesh_nodes,
            edge_list = edge_list
        )
    )
}