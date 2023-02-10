template<class Type>
class starve_nngp {
    private:
        starve_graph<Type> g;
        array<Type> w;
        covariance<Type> cv;

        vector<Type> w_node(int idx, int v);
        vector<Type> meanvec(int idx, int v);
        matrix<Type> covmat(int idx, int v);
    public:
        starve_nngp(
            const starve_graph<Type>& g,
            const array<Type>& w,
            const covariance<Type>& cv
        ) : g(g), w(w), cv(cv) {};
        starve_nngp() = default;

        Type loglikelihood();
        Type predict(
            int var,
            const vector<Type> coords,
            const vector<int> parents
        );
}

template<class Type>
vector<Type> starve_nngp<Type>::w_node(int idx, int v) {
    vector<int> vertices = g(idx);
    vector<Type> ans(vertices.size());
    for(int i = 0; i < ans.size(); i++) {
        ans(i) = w(vertices(i, 0), v);
    }
    return ans;
}

template<class Type>
matrix<Type> starve_nngp<Type>::covmat(int idx, int v) {
    vector<int> vertices = g(idx);
    matrix<Type> ss(vertices.size(), vertices.size());
    for(int i = 0; i < ss.rows(); i++) {
        for(int j = 0; j < ss.cols(); j++) {
            ss(i, j) = cv(
                g.coordinates(vertices(i)),
                g.coordinates(vertices(j)),
                v + 1, // 0 = gg, 1 = dxdx, 2 = dydy
                v + 1 // 0 = gg, 1 = dxdx, 2 = dydy
            )
        }
        ss(i, i) *= 1.001; // Add small number to main diagonal for numerical stability
    }
    return ss;
}

template<class Type>
vector<Type> starve_nngp<Type>::meanvec(int idx, int v) {
    vector<int> vertices = g(idx);
    vector<Type> mm(vertices.size());
    mm.setZero();
    return mm;
}

template<class Type>
Type starve_nngp<Type>::loglikelihood() {
    Type ll = 0.0;
    for(int i = 0; i < g.size(); i++) {
        for(int v = 0; v < w.cols(); v++) {
            vector<Type> this_w = w_node(i, v);
            vector<Type> mu = meanvec(i, v);
            matrix<Type> Sigma = covmat(i, v);
            conditional_normal<Type> cmvn(Sigma, g.from(i).size());
            ll += cmvn.loglikelihood(this_w, mu);
        }
    }
}

template<class Type>
Type starve_nngp<Type>::predict(
        int var,
        const vector<Type> coords,
        const vector<int> parents
    ) {
    vector<Type> full_w(1 + parents.size());
    for(int i = 0; i < parents.size(); i++) {
        full_w(i + 1) = w(parents(i), v);
    }
    vector<Type> mu(full_w.size());
    mu.setZero();

    matrix<Type> Sigma(full_w.size(), full_w.size());
    for(int i = 0; i < Sigma.rows(); i++) {
        for(int j = 0; j < Sigma.cols(); j++) {
            vector<Type> c1(2);
            vector<Type> c2(2);
            if( i = 0 ) {
                c1 = coords;
            } else {
                c1 = g.coordinates(parents(i - 1));
            }
            if( j = 0 ) {
                c2 = coords;
            } else {
                c2 = g.coordinates(parents(j - 1));
            }

            Sigma(i, j) = cv(c1, c2, v1 + 1, v2 + 1);
        }
    }
    conditional_normal<Type> cmvn(Sigma, parents.size());
    full_w(0) = cmvn.conditional_mean(full_w, mu)(0);

    return full_w(0);
}