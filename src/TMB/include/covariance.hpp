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
