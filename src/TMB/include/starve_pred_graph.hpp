template<class Type>
struct starve_pred_graph {
    matrix<Type> coord;
    vector<matrix<int> > parents;

    starve_pred_graph(SEXP r_list) :
            coord(asMatrix<Type>(VECTOR_ELT(r_list, 0))) {
        SEXP r_parents = VECTOR_ELT(r_list, 1);
        parents.resize(LENGTH(r_parents));
        for(int i = 0; i < LENGTH(r_parents); i++) {
            parents(i) = asMatrix<int>(VECTOR_ELT(r_parents, i));
        }
    };
};