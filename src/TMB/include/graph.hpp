template<class Type>
class nngp_graph {
  private:
    vector<Type> x_coordinates;
    vector<Type> y_coordinates;
    vector<matrix<int> > to_list;
    vector<matrix<int> > from_list;
  public:
    nngp_graph(
      const vector<Type>& x_coordinates,
      const vector<Type>& y_coordinates,
      const vector<matrix<int> >& to_list,
      const vector<matrix<int> >& from_list
    ) : x_coordinates(x_coordinates), to_list(to_list), from_list(from_list) {};
    nngp_graph(SEXP r_list) :
        x_coordinates(asVector<Type>(VECTOR_ELT(r_list, 0))),
        y_coordinates(asVector<Type>(VECTOR_ELT(r_list, 1))) {
      SEXP r_to = VECTOR_ELT(r_list, 2);
      to_list.resize(LENGTH(r_to));

      SEXP r_from = VECTOR_ELT(r_list, 3);
      from_list.resize(LENGTH(r_from));

      for(int i = 0; i < LENGTH(r_to); i++) {
        to_list(i) = asMatrix<int>(VECTOR_ELT(r_to, i));
        from_list(i) = asMatrix<int>(VECTOR_ELT(r_from, i));
      }
    };
    nngp_graph() = default;

    int size() { return to_list.size(); }

    // Get x & y coordinates
    vector<Type> get_x_coordinates() { return x_coordinates; };
    vector<Type> get_y_coordinates() { return y_coordinates; };

    // Get coordinates for specific vertex
    vector<Type> coordinates(const vector<int>& idx) {
      vector<Type> cc(2);
      cc(0) = x_coordinates(idx(0));
      cc(1) = y_coordinates(idx(1));
      return cc;
    };

    // Get to / from matrices
    matrix<int> to(int i) { return to_list(i); }
    matrix<int> from(int i) { return from_list(i); }
    matrix<int> operator() (int i) {
      matrix<int> ans(to(i).rows() + from(i).rows(), to(i).cols());
      ans << to(i), from(i);
      return ans;
    }
};