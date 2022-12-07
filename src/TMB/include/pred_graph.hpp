template<class Type>
struct pred_graph {
  vector<int> var;
  matrix<Type> coord;
  vector<matrix<int> > parents;

  pred_graph(SEXP r_list) :
      var(asVector<int>(VECTOR_ELT(r_list, 0))),
      coord(asMatrix<Type>(VECTOR_ELT(r_list, 1))) {
    SEXP r_parents = VECTOR_ELT(r_list, 2);
    parents.resize(LENGTH(r_parents));
    for(int i = 0; i < LENGTH(r_parents); i++) {
      parents(i) = asMatrix<int>(VECTOR_ELT(r_parents, i));
    }
  };
};