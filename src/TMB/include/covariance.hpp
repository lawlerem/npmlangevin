template<class Type>
class covariance {
  private:
    vector<Type> pars;
    int covar_code; // Which covariance function to use?

  public:
    // Constructor
    covariance(const vector<Type>& pars, const int& covar_code) :
      pars{pars}, covar_code{covar_code} {};
    covariance() : pars{vector<Type>()}, covar_code(0) {};

    // Squared distance function
    template<typename T> T d(const vector<T>& x1, const vector<T>& x2);

    // Covariance function, gradient, and hessian
    template<typename T> T operator() (const vector<T>& x); // x = c(x1, x2)
    template<typename T> T operator() (const vector<T>& x1, const vector<T>& x2);
    template <typename T> vector<T> gradient(const vector<T>& x1, const vector<T>& x2);
    template<typename T> matrix<T> hessian(const vector<T>& x1, const vector<T>& x2);

    template<typename T> T operator() (const vector<T>& x1, const vector<T>& x2, int v1, int v2);
};

template<class Type>
template <typename T>
T covariance<Type>::operator() (const vector<T>& x) {
  int n = x.size() / 2;
  vector<T> x1 = x.segment(0, n);
  vector<T> x2 = x.segment(n, n);

  return operator()(x1, x2);
}

template<class Type>
template<typename T>
T covariance<Type>::d(const vector<T>& x1, const vector<T>& x2) {
  vector<T> diff = x1 - x2;
  return sqrt((diff * diff).sum() + pow(10, -6));
}

template<class Type>
template<typename T>
T covariance<Type>::operator() (const vector<T>& x1, const vector<T>& x2) {
  T d = this->d(x1, x2);
  vector<T> pars = this->pars.template cast<T>();
  switch(covar_code) {
    // Exponential [sd, range]
    // case 0 : return pow(pars(0), 2) * pars(1) * exp( -d / pars(1) );
    case 0 : return pow(pars(0), 2) * exp( -d / pars(1) );
    // Gaussian [marg_sd, range]
    case 1 : return pow(pars(0), 2) * exp( -0.5 * (1 / pow(pars(1), 2)) * pow(d, 2) );
    // Matern [sd, range, nu]
    // case 2 : return pow(pars(0), 2) * pow(pars(1), 2 * pars(2)) * matern(d, pars(1), pars(2));
    case 2 : return pow(pars(0), 2) * matern(d, pars(1), pars(2));
    // Matern32 [sd, range]
    // case 3 : return pow(pars(0), 2) * pow(pars(1), 3) * (1 + d / pars(1)) * exp( -d / pars(1) );
    case 3 : return pow(pars(0), 2) * (1 + d / pars(1)) * exp( -d / pars(1) );
    // Gaussian [marg_sd, range]
    default : return pow(pars(0), 2) * exp( -0.5 * (1 / pow(pars(1), 2)) * pow(d, 2) );
  }
}

template<class Type>
template<typename T>
vector<T> covariance<Type>::gradient(const vector<T>& x1, const vector<T>& x2) {
  vector<T> x1x2(x1.size() + x2.size());
  x1x2 << x1, x2;

  return autodiff::gradient(*this, x1x2);
}

template<class Type>
template <typename T>
matrix<T> covariance<Type>::hessian(const vector<T>& x1, const vector<T>& x2) {
  vector<T> x1x2(x1.size() + x2.size());
  x1x2 << x1, x2;

  return autodiff::hessian(*this, x1x2);
}

template<class Type>
template<typename T>
T covariance<Type>::operator() (const vector<T>& x1, const vector<T>& x2, int v1, int v2) {
  T ans;
  // operator() (x, y), (x, y)
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
  if( v1 == 0 & v2 == 0 ) {
    // g_g
    ans = operator()(x1, x2);
  } else if( v1 == 0 & v2 == 1 ) {
    // g_dx
    ans = gradient(x1, x2)(2);
  } else if( v1 == 0 & v2 == 2 ) {
    // g_dy
    ans = gradient(x1, x2)(3);
  } else if( v1 == 1 & v2 == 0 ) {
    // dx_g
    ans = gradient(x1, x2)(0);
  } else if( v1 == 1 & v2 == 1 ) {
    // dx_dx
    ans = hessian(x1, x2)(0, 2);
  } else if( v1 == 1 & v2 == 2 ) {
    // dx_dy
    ans = hessian(x1, x2)(0, 3);
  } else if( v1 == 2 & v2 == 0 ) {
    // dy_g
    ans = gradient(x1, x2)(1);
  } else if( v1 == 2 & v2 == 1 ) {
    // dy_dx
    ans = hessian(x1, x2)(1, 2);
  } else if( v1 == 2 & v2 == 2 ) {
    //dy_dy
    ans = hessian(x1, x2)(1, 3);
  }
  return ans;
}
