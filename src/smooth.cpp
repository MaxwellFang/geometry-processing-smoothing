#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include "igl/edge_lengths.h"
#include <Eigen/SparseCholesky>

void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  Eigen::MatrixXd l(F.rows(), 3);
  igl::edge_lengths(V, F, l);
  Eigen::SparseMatrix<double> L(F.maxCoeff() + 1, F.maxCoeff() + 1);
  cotmatrix(l, F, L);
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> M;
  massmatrix(l, F, M);
  Eigen::MatrixXd b = M * G;
  Eigen::SparseMatrix<double> A = -lambda * L;
  Eigen::VectorXd diag = M.diagonal();
  for (int i = 0; i < A.rows(); i++) {
    A.coeffRef(i, i) += diag[i];
  }
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  U = solver.solve(b);
}
