#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export("order_adjacency_matrix")]]
Eigen::VectorXi order_adjacency_matrix(Eigen::MatrixXi &m) {
  Eigen::VectorXi order = Eigen::VectorXi::LinSpaced(
    m.rows(),
    0,
    m.rows() - 1
  );
  for(int i = 0; i < m.rows(); i++) {
    Eigen::VectorXi n_parents = m
      .topRightCorner(i, m.cols() - i)
      .colwise()
      .sum();
    int next_vertex;
    n_parents.maxCoeff(&next_vertex);
    next_vertex += i;
    m.row(i).swap(m.row(next_vertex));
    m.col(i).swap(m.col(next_vertex));
    std::iter_swap(order.data() + i, order.data() + next_vertex);
  }
  return order;
}