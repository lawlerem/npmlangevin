template<class Type>
class loc_track {
  private:
    matrix<Type> coords;
    vector<matrix<int> > field_neighbours;
    vector<Type> time;
    Type gamma;
  public:
    matrix<Type> track_gradient;
    loc_track(
      const matrix<Type>& coords,
      const vector<Type>& time,
      Type gamma
    ) : coords(coords), time(time), gamma(gamma) {};
    loc_track(
      const matrix<Type>& coords,
      const vmint<Type>& field_neighbours,
      const vector<Type>& time,
      Type gamma
    ) : coords(coords), field_neighbours(field_neighbours.x), time(time), gamma(gamma) {
      track_gradient = 0.0 * coords;
    };
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
  for(int t = 0; t < coords.rows(); t++) {
    for(int v = 0; v < coords.cols(); v++) {
      if( t > 0 ) {
        ans += dnorm(
          coords(t, v),
          coords(t - 1, v) + 0.5 * (time(t) - time(t - 1)) * track_gradient(t - 1, v),
          gamma * pow(time(t) - time(t - 1), 0.5),
          true
        );
      } else {}
    }
    matrix<int> nn = field_neighbours(t);
    for(int v = 0; v < coords.cols(); v++) {
      matrix<int> this_nn(2 * nn.rows(), 3);
      for(int i = 0; i < nn.rows(); i++) {
        // dxdx neighbours
        this_nn(i, 0) = nn(i, 0);
        this_nn(i, 1) = nn(i, 1);
        this_nn(i, 2) = 1;

        // dydy neighbours
        this_nn(i + nn.rows(), 0) = nn(i, 0);
        this_nn(i + nn.rows(), 1) = nn(i, 1);
        this_nn(i + nn.rows(), 2) = 2;
      }

      track_gradient(t, v) = field.predict(
        v + 1, // 0 = gg, 1 = dxdx, 2 = dydy
        vector<Type>(coords.row(t)),
        this_nn
      );
    }
  }
  return ans;
}

template<class Type>
matrix<Type> loc_track<Type>::simulate(nngp<Type>& field) {
  for(int t = 0; t < coords.rows(); t++) {
    for(int v = 0; v < coords.cols(); v++) {
      if( t > 0 ) {
        coords(t, v) = rnorm(
          coords(t - 1, v) + 0.5 * (time(t) - time(t - 1)) * track_gradient(t - 1, v),
          gamma * pow(time(t) - time(t - 1), 0.5)
        );
      } else {}
    }
    matrix<int> nn = field.find_nearest_four(coords.row(t));
    for(int v = 0; v < coords.cols(); v++) {
      matrix<int> this_nn(2 * nn.rows(), 3);
      for(int i = 0; i < nn.rows(); i++) {
        // dxdx neighbours
        this_nn(i, 0) = nn(i, 0);
        this_nn(i, 1) = nn(i, 1);
        this_nn(i, 2) = 1;

        // dydy neighbours
        this_nn(i + nn.rows(), 0) = nn(i, 0);
        this_nn(i + nn.rows(), 1) = nn(i, 1);
        this_nn(i + nn.rows(), 2) = 2;
      }
      track_gradient(t, v) = field.predict(
        v + 1, // 0 = gg, 1 = dxdx, 2 = dydy
        vector<Type>(coords.row(t)),
        this_nn
      );
    }
  }
  return coords;
}
