//  Read in graph
//  list(
//      to = c(1, 2, 3),
//      from = c(1, 2, 3)
//  )
template<class Type>
struct directed_graph {
    vector<vector<vector<int> > > dag;
    directed_graph(SEXP edge_list) {
        dag.resize(LENGTH(edge_list));
        for(int i = 0; i < LENGTH(edge_list); i++) {
            SEXP v = VECTOR_ELT(edge_list, i);
            dag(i).resize(2);

            vector<int> to = asVector<int>(VECTOR_ELT(v, 0));
            dag(i)(0).resizeLike(to);
            dag(i)(0) = to;

            vector<int> from = asVector<int>(VECTOR_ELT(v, 1));
            dag(i)(1).resizeLike(from);
            dag(i)(1) = from;
        }
    }
};

//  Read in animal track
//  list(
//      lat = 1,
//      lon = 1,
//      time = 1,
//      class = a,
//      parents = c(1, 2, 3),
//  )