template<class Type>
class loc_track {
  private:
    matrix<Type> coords;
    matrix<Type> field_values;
    vector<matrix<int> > field_neighbours;
    vector<Type> time;
    Type gamma;
  public:
    loc_track(
      const matrix<Type>& coords,
      const vector<Type>& time,
      Type gamma
    ) : coords(coords), time(time), gamma(gamma) {};
    loc_track(
      const matrix<Type>& coords,
      const matrix<Type>& field_values,
      const vector<matrix<int> >& field_neighbours,
      const vector<Type>& time,
      Type gamma
    ) : coords(coords), field_values(field_values), field_neighbours(field_neighbours), time(time), gamma(gamma) {};
    loc_track() = default;

    // Return location i
    vector<Type> operator() (int i) {
      return vector<Type>(coords.row(i));
    };

    // If no field, just simulate a random walk
    Type loglikelihood();
    Type loglikelihood(const nngp<Type>& field);

    // If field, then use Langevin diffusion
    matrix<Type> simulate();
    matrix<Type> simulate(const nngp<Type>& field);
};

template<class Type>
Type loc_track<Type>::loglikelihood() {
  Type ans = 0.0;
  for(int t = 1; t < coords.rows(); t++) {
    for(int v = 0; v < coords.cols(); v++) {
      ans += dnorm(
        coords(t, v),
        coords(t - 1, v),
        gamma * pow(time(t) - time(t - 1), 0.5),
        true
      );
    }
  }
  return ans;
}

template<class Type>
matrix<Type> loc_track<Type>::simulate() {
  for(int t = 1; t < coords.rows(); t++) {
    for(int v = 0; v < coords.cols(); v++) {
      coords(t, v) = rnorm(
        coords(t - 1, v),
        gamma * pow(time(t) - time(t - 1), 0.5)
      );
    }
  }
  return coords;
}

template<class Type>
Type loc_track<Type>::loglikelihood(const nngp<Type>& field) {
  Type ans = 0.0;
  for(int t = 1; t < coords.rows(); t++) {
    for(int v = 0; v < coords.cols(); v++) {
      Type grad = field.predict(
        field_values(t - 1, v),
        v + 1, // 0 = gg, 1 = dxdx, 2 = dydy
        vector<Type>(coords.row(t - 1)),
        field_neighbours(t - 1),
        ans
      );
      ans += dnorm(
        coords(t, v),
        coords(t - 1, v) + 0.5 * (time(t) - time(t - 1)) * pow(gamma, 2) * grad,
        gamma * pow(time(t) - time(t - 1), 0.5),
        true
      );
    }
  }
  return ans;
}

template<class Type>
matrix<Type> loc_track<Type>::simulate(const nngp<Type>& field) {
  Type foo = 0.0;
  for(int t = 1; t < coords.rows(); t++) {
    for(int v = 0; v < coords.cols(); v++) {
      Type grad = field.predict(
        field_values(t - 1, v),
        v + 1, // 0 = gg, 1 = dxdx, 2 = dydy
        vector<Type>(coords.row(t - 1)),
        field_neighbours(t - 1),
        foo,
        true // Only return mean prediction
      );
      coords(t, v) = rnorm(
        coords(t - 1, v) + 0.5 * (time(t) - time(t - 1)) * pow(gamma, 2) * grad,
        gamma * pow(time(t) - time(t - 1), 0.5)
      );
    }
  }
  return coords;
}
