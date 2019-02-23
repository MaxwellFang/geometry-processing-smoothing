#include "cotmatrix.h"

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  // l: F by 3 lists of edges [1, 2], [2, 0], [0, 1]
  // F: faces
  // Assume that L has already been constructed
  L.resize(F.maxCoeff() + 1, F.maxCoeff() + 1);
  for (int i = 0; i < F.rows(); i++) {
    int v0 = F(i, 0);
    int v1 = F(i, 1);
    int v2 = F(i, 2);
    double a = l(i, 0);
    double b = l(i, 1);
    double c = l(i, 2);
    double s = (a + b + c) / 2;
    double area = sqrt(s * (s - a) * (s - b) * (s - c));
    double sinA = 2 * area / (b * c);
    double sinB = 2 * area / (a * c);
    double sinC = 2 * area / (a * b);
    double cosA = (pow(b, 2) + pow(c, 2) - pow(a, 2)) / (2 * b * c);
    double cosB = (pow(a, 2) + pow(c, 2) - pow(b, 2)) / (2 * a * c);
    double cosC = (pow(a, 2) + pow(b, 2) - pow(c, 2)) / (2 * a * b);
    double cotA = cosA / sinA;
    double cotB = cosB / sinB;
    double cotC = cosC / sinC;
    L.coeffRef(v1, v2) += 0.5 * cotA;
    L.coeffRef(v0, v2) += 0.5 * cotB;
    L.coeffRef(v0, v1) += 0.5 * cotC;
    L.coeffRef(v2, v1) += 0.5 * cotA;
    L.coeffRef(v2, v0) += 0.5 * cotB;
    L.coeffRef(v1, v0) += 0.5 * cotC;
    L.coeffRef(v0, v0) -= (0.5 * cotB + 0.5 * cotC);
    L.coeffRef(v1, v1) -= (0.5 * cotA + 0.5 * cotC);
    L.coeffRef(v2, v2) -= (0.5 * cotA + 0.5 * cotB);

  }
}

