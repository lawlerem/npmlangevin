template<class Type>
struct vmint {
  vector<matrix<int> > x;
  vmint(SEXP r_list) {
    x.resize(LENGTH(r_list));
    for(int i = 0; i < x.size(); i++) {
      x(i) = asMatrix<int>(VECTOR_ELT(r_list, i));
    }
  }
};

template<class Type>
struct vvint {
  vector<vector<int> > x;
  vvint(SEXP r_list) {
    x.resize(LENGTH(r_list));
    for(int i = 0; i < x.size(); i++) {
      x(i) = asVector<int>(VECTOR_ELT(r_list, i));
    }
  }
};

template<class Type>
struct refSorter {
  refSorter(const Eigen::Matrix<Type, Dynamic, 1> &d) : d_(d) {}
  bool operator () (const int a, const int b) {
    return d_(a) < d_(b);
  }
  const Eigen::Matrix<Type, Dynamic, 1> d_;
};
