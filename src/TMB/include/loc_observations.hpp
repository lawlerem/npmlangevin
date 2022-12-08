template<class Type>
struct loc_observations {
  matrix<Type> coords;
  vector<int> loc_class;
  vector<int> track_idx;
  matrix<Type> K;

  loc_observations(SEXP r_list) :
    coords(asMatrix<Type>(VECTOR_ELT(r_list, 0))),
    loc_class(asVector<int>(VECTOR_ELT(r_list, 1))),
    track_idx(asVector<int>(VECTOR_ELT(r_list, 2))),
    K(asMatrix<Type>(VECTOR_ELT(r_list, 3)))
    {};

  matrix<Type> conjugate(matrix<Type> Sigma, int q) {
    matrix<Type> diagK(2, 2);
    diagK << K(q, 0), 0.0, 0.0, K(q, 1);
    return diagK * Sigma * diagK;
  };

  Type loglikelihood(const matrix<Type>& Sigma, loc_track<Type>& true_loc);
  matrix<Type> simulate(const matrix<Type>& Sigma, loc_track<Type>& true_loc);
  int size() { return coords.rows(); };
};

template<class Type>
Type loc_observations<Type>::loglikelihood(
    const matrix<Type>& Sigma,
    loc_track<Type>& true_loc) {
  Type ans = 0.0;

  vector<MVNORM_t<Type> > obs_mvns(K.rows());
  for( int q = 0; q < obs_mvns.size(); q++) {
    obs_mvns(q) = MVNORM_t<Type>(conjugate(Sigma, q));
  }

  for(int i = 0; i < coords.rows(); i++) {
    ans -= obs_mvns(loc_class(i))(vector<Type>(coords.row(i)) - true_loc(track_idx(i)));
  }

  return ans;
}

template<class Type>
matrix<Type> loc_observations<Type>::simulate(
    const matrix<Type>& Sigma,
    loc_track<Type>& true_loc) {
  vector<MVNORM_t<Type> > obs_mvns(K.rows());
  for( int q = 0; q < obs_mvns.size(); q++) {
    obs_mvns(q) = MVNORM_t<Type>(conjugate(Sigma, q));
  }

  for(int i = 0; i < coords.rows(); i++) {
    coords.row(i) = obs_mvns(loc_class(i)).simulate() + true_loc(track_idx(i));
  }

  return coords;
}