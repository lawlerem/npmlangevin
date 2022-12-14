template<class Type>
class loc_track {
  private:
    matrix<Type> coords;
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
      const vector<matrix<int> >& field_neighbours,
      const vector<Type>& time,
      Type gamma
    ) : coords(coords), field_neighbours(field_neighbours), time(time), gamma(gamma) {};
    loc_track(
      const matrix<Type>& coords,
      const vmint<Type>& field_neighbours,
      const vector<Type>& time,
      Type gamma
    ) : coords(coords), field_neighbours(field_neighbours.x), time(time), gamma(gamma) {};
    loc_track() = default;

    // Return location i
    vector<Type> operator() (int i) {
      return vector<Type>(coords.row(i));
    };

    // If no field, then use random walk
    Type loglikelihood();
    matrix<Type> simulate();

    // If field, then use Langevin diffusion
    Type loglikelihood(nngp<Type>& field);
    matrix<Type> simulate(nngp<Type>& field);
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
Type loc_track<Type>::loglikelihood(nngp<Type>& field) {
  Type ans = 0.0;
  Type foo = 0.0;
  for(int t = 1; t < coords.rows(); t++) {
    matrix<int> neighbours = field.find_nearest_four(coords.row(t - 1));
    for(int v = 0; v < coords.cols(); v++) {
      matrix<int> this_neighbours(2 * neighbours.rows(), 3);
      for(int i = 0; i < neighbours.rows(); i++) {
        // dxdx neighbours
        this_neighbours(i, 0) = neighbours(i, 0);
        this_neighbours(i, 1) = neighbours(i, 1);
        this_neighbours(i, 2) = 2;

        // dydy neighbours
        this_neighbours(i + neighbours.rows(), 0) = neighbours(i, 0);
        this_neighbours(i + neighbours.rows(), 1) = neighbours(i, 1);
        this_neighbours(i + neighbours.rows(), 2) = v;
      }
      Type grad = field.predict(
        foo,
        v + 1, // 0 = gg, 1 = dxdx, 2 = dydy
        vector<Type>(coords.row(t - 1)),
        this_neighbours,
        ans,
        true // Only return mean prediction
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
matrix<Type> loc_track<Type>::simulate(nngp<Type>& field) {
  Type foo = 0.0;
  Type bar = 0.0;
  for(int t = 1; t < coords.rows(); t++) {
    matrix<int> neighbours = field.find_nearest_four(coords.row(t - 1));
    for(int v = 0; v < coords.cols(); v++) {
      matrix<int> this_neighbours(2 * neighbours.rows(), 3);
      for(int i = 0; i < neighbours.rows(); i++) {
        // dxdx neighbours
        this_neighbours(i, 0) = neighbours(i, 0);
        this_neighbours(i, 1) = neighbours(i, 1);
        this_neighbours(i, 2) = 2;

        // dydy neighbours
        this_neighbours(i + neighbours.rows(), 0) = neighbours(i, 0);
        this_neighbours(i + neighbours.rows(), 1) = neighbours(i, 1);
        this_neighbours(i + neighbours.rows(), 2) = v;
      }
      Type grad = field.predict(
        foo,
        v + 1, // 0 = gg, 1 = dxdx, 2 = dydy
        vector<Type>(coords.row(t - 1)),
        this_neighbours,
        bar,
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
