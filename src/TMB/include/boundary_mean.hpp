template<class Type>
class boundary_mean {
  private:
    vector<Type> xlim;
    vector<Type> ylim;
    Type sharpness;
  public:
    boundary_mean(
      const vector<Type>& xlim,
      const vector<Type>& ylim,
      Type sharpness
    ) : xlim(xlim), ylim(ylim), sharpness(sharpness) {};
    boundary_mean() = default;

    template<typename T> T operator() (const vector<T>& x);
    template<typename T> vector<T> gradient(const vector<T>& x);
    template<typename T> T operator() (const vector<T>& x, int v);
};

template<class Type>
template<typename T>
T boundary_mean<Type>::operator() (const vector<T>& x) {
  vector<T> txlim = xlim.template cast<T>();
  vector<T> tylim = ylim.template cast<T>();

  // T xcomp = pow((2.0 * x(0) - txlim(1) - txlim(0)) / (txlim(1) - txlim(0)), (T)sharpness);
  // T ycomp = pow((2.0 * x(1) - tylim(1) - tylim(0)) / (tylim(1) - tylim(0)), (T)sharpness);

  // return -1.0 * (xcomp + ycomp);

  T xlower = -1.0 * exp( -1.0 * (T)sharpness * (x(0) - txlim(0)) );
  T xupper = -1.0 * exp( (T)sharpness * (x(0) - txlim(1)) );
  T ylower = -1.0 * exp( -1.0 * (T)sharpness * (x(1) - tylim(0)) );
  T yupper = -1.0 * exp( (T)sharpness * (x(1) - tylim(1)) );

  return xlower + xupper + ylower + yupper;
}

template<class Type>
template<typename T>
vector<T> boundary_mean<Type>::gradient(const vector<T>& x) {
  return autodiff::gradient(*this, x);
}

template<class Type>
template<typename T>
T boundary_mean<Type>::operator() (const vector<T>& x, int v) {
  T ans;
  if( v == 0 ) {
    ans = operator()(x);
  } else {
    ans = gradient(x)(v - 1);
  }
  return ans;
}