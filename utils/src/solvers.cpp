#include "solvers.h"
#include "geometrycentral/numerical/linear_solvers.h"
/*
 * Compute the inverse of a sparse diagonal matrix.
 *
 * Input: A sparse diagonal matrix <M>.
 * Returns: The inverse of M, which is also a sparse diagonal matrix.
 */
SparseMatrix<double> sparseInverseDiagonal(SparseMatrix<double>& M) {

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    SparseMatrix<double> inv(M.rows(), M.cols());
    for (int i = 0; i < M.rows(); i++) {
        tripletList.push_back(T(i, i, 1.0 / M.coeffRef(i, i)));
    }
    inv.setFromTriplets(tripletList.begin(), tripletList.end());
    return inv;
}

/*
 * Computes the residual of Ax - λx, where x has unit norm and λ = x.Ax.
 *
 * Input: <A>, the complex sparse matrix whose eigendecomposition is being computed; and <x>, the current guess for the
 * smallest eigenvector
 * Returns: The residual
 */
double residual(const SparseMatrix<std::complex<double>>& A, const Vector<std::complex<double>>& x) {

    // TODO
    
    return (A*x - (x.transpose()*A*x)*x).norm(); // placeholder
}

/*
 * Solves Ax = λx, where λ is the smallest nonzero eigenvalue of A, and x is the corresponding eigenvector.
 *
 * Input: <A>, the complex positive definite sparse matrix whose eigendecomposition is being computed.
 * Returns: The smallest eigenvector of A.
 */
Vector<std::complex<double>> solveInversePowerMethod(SparseMatrix<std::complex<double>>& A) {

    // TODO
    Vector<std::complex<double>> y = Vector<std::complex<double>>::Random(A.cols()).normalized();
    //std::cout << y << "y" << std::endl;
    PositiveDefiniteSolver<std::complex<double>> solver(A);
    double r = residual(A, y);

    while (residual(A, y) > 1e-4)
    {
        r = residual(A, y);
        std::cout << "r = " <<  r << std::endl;
        y = solver.solve(y);
        std::complex<double> center{0., 0.};
        for (int i = 0; i < y.size(); ++i) center += y[i];
        center /= y.rows();
        for (int i = 0; i < y.size(); ++i) y[i] -= center;
        
        y = y.normalized();
    }
    std::cout << y << "y" << std::endl;

    return y;
}