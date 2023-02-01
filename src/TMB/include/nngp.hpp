template<class Type>
class nngp {
  private:
    nngp_graph<Type> g;
    array<Type> w;
    boundary_mean<Type> boundary;
    covariance<Type> cv;

    vector<Type> w_node(int idx);
    vector<Type> meanvec(int idx);
    matrix<Type> covmat(int idx);
  public:
    nngp(
      const nngp_graph<Type>& g,
      const array<Type>& w,
      const boundary_mean<Type>& boundary,
      const covariance<Type>& cv
    ) : g(g), w(w), boundary(boundary), cv(cv) {};
    nngp() = default;

    Type loglikelihood();
    array<Type> simulate();
    Type predict(int var, const vector<Type> coords, const matrix<int> parents);
    matrix<int> find_nearest_four(vector<Type> coord);
};

template<class Type>
vector<Type> nngp<Type>::w_node(int idx) {
  matrix<int> vertices = g(idx);
  vector<Type> ans(vertices.rows());
  for(int i = 0; i < ans.size(); i++) {
    ans(i) = w(vertices(i, 0), vertices(i, 1), vertices(i, 2));
  }
  return ans;
}

template<class Type>
matrix<Type> nngp<Type>::covmat(int idx) {
  matrix<int> vertices = g(idx);
  matrix<Type> ss(vertices.rows(), vertices.rows());
  for(int i = 0; i < ss.rows(); i++) {
    for(int j = 0; j < ss.cols(); j++) {
      ss(i, j) = cv(
        g.coordinates(vector<int>(vertices.row(i))),
        g.coordinates(vector<int>(vertices.row(j))),
        vertices(i, 2),
        vertices(j, 2)
      );
    }
    ss(i, i) *= 1.001; // Add small number to main diagonal for numerical stability
  }
  return ss;
}

template<class Type>
vector<Type> nngp<Type>::meanvec(int idx) {
  matrix<int> vertices = g(idx);
  vector<Type> mm(vertices.rows());
  for(int i = 0; i < mm.size(); i++) {
    mm(i) = boundary(
      g.coordinates(vector<int>(vertices.row(i))),
      vertices(i, 2)
    );
  }
  return mm;
}

template<class Type>
Type nngp<Type>::loglikelihood() {
  Type ll = 0.0;
  for(int i = 0; i < g.size(); i++) {
    vector<Type> this_w = w_node(i);
    vector<Type> mu = meanvec(i);
    matrix<Type> Sigma = covmat(i);
    conditional_normal<Type> cmvn(Sigma, g.from(i).rows());
    ll += cmvn.loglikelihood(this_w, mu);
  }
  return ll;
}

template<class Type>
array<Type> nngp<Type>::simulate() {
  for(int i = 0; i < g.size(); i++) {
    vector<Type> this_w = w_node(i);
    vector<Type> mu = meanvec(i);
    matrix<Type> Sigma = covmat(i);
    conditional_normal<Type> cmvn(Sigma, g.from(i).rows());
    this_w = cmvn.simulate(this_w, mu);
    for(int j = 0; j < g.to(i).rows(); j++) {
      w(g.to(i)(j, 0), g.to(i)(j, 1), g.to(i)(j, 2)) = this_w(j);
    }
  }
  return w;
}

template<class Type>
Type nngp<Type>::predict(
      int var,
      const vector<Type> coords,
      const matrix<int> parents
    ) {
  vector<Type> full_w(1 + parents.rows());
  for(int i = 0; i < parents.rows(); i++) {
    full_w(i + 1) = w(parents(i, 0), parents(i, 1), parents(i, 2));
  }
  vector<Type> mu(full_w.size());
  for(int i = 0; i < mu.size(); i++) {
    if( i == 0 ) {
      mu(i) = boundary(coords, var);
    } else {
      mu(i) = boundary(
        g.coordinates(vector<int>(parents.row(i - 1))),
        parents(i - 1, 2)
      );
    }
  }

  matrix<Type> Sigma(full_w.size(), full_w.size());
  for( int i = 0; i < Sigma.rows(); i++ ) {
    for( int j = 0; j < Sigma.cols(); j++ ) {
      vector<Type> c1(2);
      vector<Type> c2(2);
      int v1;
      int v2;
      if( i == 0 ) {
        c1 = coords;
        v1 = var;
      } else {
        c1 = g.coordinates(vector<int>(parents.row(i - 1)));
        v1 = parents(i - 1, 2);
      }
      if( j == 0 ) {
        c2 = coords;
        v2 = var;
      } else {
        c2 = g.coordinates(vector<int>(parents.row(j - 1)));
        v2 = parents(j - 1, 2);
      }

      Sigma(i, j) = cv(c1, c2, v1, v2);
    }
  }
  conditional_normal<Type> cmvn(Sigma, parents.rows());
  full_w(0) = cmvn.conditional_mean(full_w, mu)(0);

  return full_w(0);
}

template<class Type>
matrix<int> nngp<Type>::find_nearest_four(vector<Type> coord) {
  // Get nearest x coordinates
  Eigen::Matrix<Type, Dynamic, 1> xcoord_d = abs(g.get_x_coordinates() - coord(0));
  Eigen::VectorXi xind = Eigen::VectorXi::LinSpaced(
    g.get_x_coordinates().size(),
    0,
    g.get_x_coordinates().size() - 1
  );
  std::partial_sort(
    xind.data(),
    xind.data() + 4,
    xind.data() + xind.size(),
    refSorter<Type>(xcoord_d)
  );

  // Get nearest y coordinates
  Eigen::Matrix<Type, Dynamic, 1> ycoord_d = abs(g.get_y_coordinates() - coord(1));
  Eigen::VectorXi yind = Eigen::VectorXi::LinSpaced(
    g.get_y_coordinates().size(),
    0,
    g.get_y_coordinates().size() - 1
  );
  std::partial_sort(
    yind.data(),
    yind.data() + 4,
    yind.data() + yind.size(),
    refSorter<Type>(ycoord_d)
  );

  matrix<int> nn(4, 2);
  nn(0, 0) = xind(0);
  nn(1, 0) = xind(0);
  nn(2, 0) = xind(1);
  nn(3, 0) = xind(1);

  nn(0, 1) = yind(0);
  nn(1, 1) = yind(1);
  nn(2, 1) = yind(0);
  nn(3, 1) = yind(1);

  return nn;
}