#include "massmatrix.h"

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
  Eigen::VectorXd diag = Eigen::VectorXd::Zero(F.maxCoeff() + 1);
  for (int i = 0; i < F.rows(); i++) {
    int v0 = F(i, 0);
    int v1 = F(i, 1);
    int v2 = F(i, 2);
    double a = l(i, 0);
    double b = l(i, 1);
    double c = l(i, 2);
    double s = (a + b + c) / 2;
    double area = sqrt(s * (s - a) * (s - b) * (s - c));
    diag[v0] += area / 3;
    diag[v1] += area / 3;
    diag[v2] += area / 3;
  }
  M = diag.asDiagonal();
  
    

}

