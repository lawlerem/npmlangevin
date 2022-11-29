#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type covariance_exploration(objective_function<Type>* obj) {
  DATA_MATRIX(x); // Matrix of coordinates [(lon1, lat1), (lon2, lat2), etc.]
  DATA_INTEGER(cv_code);
  PARAMETER_VECTOR(working_cv_pars);
  PARAMETER(dummy);

  vector<Type> cv_pars = exp(working_cv_pars);
  covariance<Type> f {cv_pars, cv_code};

  vector<Type> x1(2);
  x1 << 0.0, 0.0;


  vector<Type> gg(x.rows());
  vector<Type> g_dx(x.rows());
  vector<Type> g_dy(x.rows());
  vector<Type> dx_g(x.rows());
  vector<Type> dx_dx(x.rows());
  vector<Type> dx_dy(x.rows());
  vector<Type> dy_g(x.rows());
  vector<Type> dy_dx(x.rows());
  vector<Type> dy_dy(x.rows());
  for(int i = 0; i < x.rows(); i++) {
    vector<Type> x2 = x.row(i);
    gg(i) = f(x1, x2);
    g_dx(i) = f.gradient(x1, x2)(x.cols() + 0);
    g_dy(i) = f.gradient(x1, x2)(x.cols() + 1);
    dx_g(i) = f.gradient(x1, x2)(0);
    dx_dx(i) = f.hessian(x1, x2)(0, x.cols() + 0);
    dx_dy(i) = f.hessian(x1, x2)(0, x.cols() + 1);
    dy_g(i) = f.gradient(x1, x2)(1);
    dy_dx(i) = f.hessian(x1, x2)(1, x.cols() + 0);
    dy_dy(i) = f.hessian(x1, x2)(1, x.cols() + 1);
  }
  REPORT(gg);
  REPORT(g_dx);
  REPORT(g_dy);
  REPORT(dx_g);
  REPORT(dx_dx);
  REPORT(dx_dy);
  REPORT(dy_g);
  REPORT(dy_dx);
  REPORT(dy_dy);

  return pow(dummy, 2);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
