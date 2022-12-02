template<class Type>
class nngp {
  private:
    nngp_graph<Type> g;
    array<Type> w;
    covariance<Type> cv;

    vector<Type> w_node(int idx);
    matrix<Type> covmat(int idx);
  public:
    nngp(
      const nngp_graph<Type>& g,
      const array<Type>& w,
      const covariance<Type>& cv
    ) : g(g), w(w), cv(cv) {};
    nngp() = default;

    Type loglikelihood();
    array<Type> simulate();
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
      // () (x, y), (x, y)
      // (x, y), (x, y) [g_g]
      //
      // .gradient (x, y), (x, y)
      // (0) = (dx, y), (x, y) [dx_g]
      // (1) = (x, dy), (x, y) [dy_g]
      // (2) = (x, y), (dx, y) [g_dx]
      // (3) = (x, y), (x, dy) [g_dy]
      //
      // .hessian (x, y), (x, y)
      // (0, 2) = (dx, y), (dx, y) [dx_dx]
      // (0, 3) = (dx, y), (x, dy) [dx_dy]
      // (1, 2) = (x, dy), (dx, y) [dy_dx]
      // (1, 3) = (x, dy), (x, dy) [dy_dy]
      if( vertices(i, 2) == 0 & vertices(j, 2) == 0 ) {
        // g_g
        ss(i, j) = cv(
          g.coordinates(vector<int>(vertices.row(i))),
          g.coordinates(vector<int>(vertices.row(j)))
        );
      } else if( vertices(i, 2) == 0 & vertices(j, 2) == 1 ) {
        // g_dx
        ss(i, j) = cv.gradient(
          g.coordinates(vector<int>(vertices.row(i))),
          g.coordinates(vector<int>(vertices.row(j)))
        )(2);
      } else if( vertices(i, 2) == 0 & vertices(j, 2) == 2 ) {
        // g_dy
        ss(i, j) = cv.gradient(
          g.coordinates(vector<int>(vertices.row(i))),
          g.coordinates(vector<int>(vertices.row(j)))
        )(3);
      } else if( vertices(i, 2) == 1 & vertices(j, 2) == 0 ) {
        // dx_g
        ss(i, j) = cv.gradient(
          g.coordinates(vector<int>(vertices.row(i))),
          g.coordinates(vector<int>(vertices.row(j)))
        )(0);
      } else if( vertices(i, 2) == 1 & vertices(j, 2) == 1 ) {
        // dx_dx
        ss(i, j) = cv.hessian(
          g.coordinates(vector<int>(vertices.row(i))),
          g.coordinates(vector<int>(vertices.row(j)))
        )(0, 2);
      } else if( vertices(i, 2) == 1 & vertices(j, 2) == 2 ) {
        // dx_dy
        ss(i, j) = cv.hessian(
          g.coordinates(vector<int>(vertices.row(i))),
          g.coordinates(vector<int>(vertices.row(j)))
        )(0, 3);
      } else if( vertices(i, 2) == 2 & vertices(j, 2) == 0 ) {
        // dy_g
        ss(i, j) = cv.gradient(
          g.coordinates(vector<int>(vertices.row(i))),
          g.coordinates(vector<int>(vertices.row(j)))
        )(1);
      } else if( vertices(i, 2) == 2 & vertices(j, 2) == 1 ) {
        // dy_dx
        ss(i, j) = cv.hessian(
          g.coordinates(vector<int>(vertices.row(i))),
          g.coordinates(vector<int>(vertices.row(j)))
        )(1, 2);
      } else if( vertices(i, 2) == 2 & vertices(j, 2) == 2 ) {
        //dy_dy
        ss(i, j) = cv.hessian(
          g.coordinates(vector<int>(vertices.row(i))),
          g.coordinates(vector<int>(vertices.row(j)))
        )(1, 3);
      }
    }
  }
  return ss;
}

template<class Type>
Type nngp<Type>::loglikelihood() {
  Type ll = 0.0;
  for(int i = 0; i < g.size(); i++) {
    vector<Type> this_w = w_node(i);
    vector<Type> mu = 0.0 * this_w;
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
    vector<Type> mu = 0.0 * this_w;
    matrix<Type> Sigma = covmat(i);
    conditional_normal<Type> cmvn(Sigma, g.from(i).rows());
    this_w = cmvn.simulate(this_w, mu);
    for(int j = 0; j < g.to(i).rows(); j++) {
      w(g.to(i)(j, 0), g.to(i)(j, 1), g.to(i)(j, 2)) = this_w(j);
    }
  }
  return w;
}