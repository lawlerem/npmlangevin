template<class Type>
class starve_graph {
    private:
        matrix<Type> coordinates;
        vector<matrix<int> > to_list;
        vector<matrix<int> > from_list;
    public:
        starve_graph(
            const matrix<Type>& coordinates,
            const vector<vector<int> >& to_list,
            const vector<vector<int> >& from_list
        ) : coordinates(coordinates), to_list(to_list), from_list(from_list) {};
        starve_graph(SEXP r_list) :
                coordinates(asMatrix<Type>(VECTOR_ELT(r_list, 0))) {
            SEXP r_to = VECTOR_ELT(r_list, 1);
            to_list.resize(LENGTH(r_to));

            SEXP r_from = VECTOR_ELT(r_list, 2);
            from_list.resize(LENGTH(r_from));

            for(int i = 0; i < LENGTH(r_to); i++) {
                to_list(i) = asVector<int>(VECTOR_ELT(r_to, i));
                from_list(i) = asVector<int>(VECTOR_ELT(r_from, i));
            }
        };
        starve_graph() = default;

        int size() { return to_list.size(); }

        matrix<Type> get_coordinates() { return coordinates; };
        vector<Type> get_coordinates(int i) { return vector<Type>(coordinates.row(i)); };

        vector<int> to(int i) { return to_list(i); }
        vector<int> from(int i) { return from_list(i); }
        vector<int> operator() (int i) {
            vector<int> ans(to(i).rows() + from(i).rows());
            ans << to(i), from(i);
            return ans;
        }
};