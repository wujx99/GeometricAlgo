#include "PoissonReconstruct.h"
#include "Eigen/src/SparseLU/SparseLU.h"
 

void PoissonRec::poisson_rec(const Eigen::MatrixXd& P, const Eigen::MatrixXd& N, Eigen::Vector3i& dims,
                             std::pair<Eigen::Vector3d, Eigen::Vector3d>& bbox, double& iso, Eigen::MatrixXd& x,
                             Eigen::VectorXd& g) {

    // number of input points
    const int n = P.rows();
    // Grid dimensions
    int nx, ny, nz;
    // Maximum extent (side length of bounding box) of points
    double max_extent = (P.colwise().maxCoeff() - P.colwise().minCoeff()).maxCoeff();
    // padding: number of cells beyond bounding box of input points
    const double pad = 8;
    // choose grid spacing (h) so that shortest side gets 30+2*pad samples
    double h = max_extent / double(30 + 2 * pad);
    // Place bottom-left-front corner of grid at minimum of points minus padding
    Eigen::RowVector3d corner = P.colwise().minCoeff().array() - pad * h;
    // Grid dimensions should be at least 3
    nx = std::max((P.col(0).maxCoeff() - P.col(0).minCoeff() + (2. * pad) * h) / h, 3.);
    ny = std::max((P.col(1).maxCoeff() - P.col(1).minCoeff() + (2. * pad) * h) / h, 3.);
    nz = std::max((P.col(2).maxCoeff() - P.col(2).minCoeff() + (2. * pad) * h) / h, 3.);
    // for render
    dims.x() = nx;
    dims.y() = ny;
    dims.z() = nz;
    bbox.first.x() = corner.x();
    bbox.first.y() = corner.y();
    bbox.first.z() = corner.z();
    bbox.second = bbox.first + Eigen::Vector3d{(nx - 1) * h, (ny - 1) * h, (nz - 1) * h};

    // Compute positions of grid nodes
    x.resize(nx * ny * nz, 3);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                // Convert subscript to index
                const auto ind = i + nx * (j + k * ny);
                x.row(ind) = corner + h * Eigen::RowVector3d(i, j, k);
            }
        }
    }

    int rowSize = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1);
    int colSize = nx * ny * nz;
    Eigen::SparseMatrix<double, Eigen::RowMajor> G(rowSize, colSize);
    fd_grad(nx, ny, nz, h, G);

    Eigen::RowVector3d cornerX = corner + Eigen::RowVector3d{h / 2., 0, 0};
    Eigen::RowVector3d cornerY = corner + Eigen::RowVector3d{0, h / 2., 0};
    Eigen::RowVector3d cornerZ = corner + Eigen::RowVector3d{0., 0, h / 2.};
    Eigen::SparseMatrix<double> Wx, Wy, Wz, W;
    fd_interpolate(nx - 1, ny, nz, h, cornerX, P, Wx);
    fd_interpolate(nx, ny - 1, nz, h, cornerY, P, Wy);
    fd_interpolate(nx, ny, nz - 1, h, cornerZ, P, Wz);
    fd_interpolate(nx, ny, nz, h, corner, P, W);
    Eigen::VectorXd vx = Wx.transpose() * N.col(0);
    Eigen::VectorXd vy = Wy.transpose() * N.col(1);
    Eigen::VectorXd vz = Wz.transpose() * N.col(2);
    Eigen::VectorXd v(rowSize);
    v << vx, vy, vz;
    Eigen::BiCGSTAB<SparseMatrix<double>> solver;
    solver.compute(G.transpose() * G);
    g = solver.solve(G.transpose() * v);
    iso = (W * g).sum() / g.rows();
}

void PoissonRec::fd_grad(const int nx, const int ny, const int nz, const double h,
                         Eigen::SparseMatrix<double, Eigen::RowMajor>& G) 
{
    int sizeX = (nx - 1) * ny * nz;
    int sizeY = nx * (ny - 1) * nz;
    int sizeZ = nx * ny * (nz - 1);
    int sizeP = nx * ny * nz;
    if (G.rows() != (sizeX+sizeY+sizeZ))
    {
        std::cout << "G row size don't equal to the result size" << std::endl;
    }
    Eigen::SparseMatrix<double> gradX(sizeX, sizeP);
    Eigen::SparseMatrix<double> gradY(sizeY, sizeP);
    Eigen::SparseMatrix<double> gradZ(sizeZ, sizeP);
    fd_partial_derivative(nx, ny, nz, h, 0, gradX);
    fd_partial_derivative(nx, ny, nz, h, 1, gradY);
    fd_partial_derivative(nx, ny, nz, h, 2, gradZ);
    G.topRows(sizeX) = gradX;
    G.middleRows(sizeX, sizeY) = gradY;
    G.bottomRows(sizeZ) = gradZ;
}

void PoissonRec::fd_interpolate(const int nx, const int ny, const int nz, const double h,
                                const Eigen::RowVector3d& corner, const Eigen::MatrixXd& P,
                                Eigen::SparseMatrix<double>& W) 
{
    assert(P.cols() == 3); 
    W.resize(P.rows(), nx * ny * nz);
    std::vector<Eigen::Triplet<double>> tripleLists;
    for (int i = 0; i < P.rows(); ++i)
    {
        auto relP = P.row(i) - corner;
        int idxX = relP.x() / h;
        int idxY = relP.y() / h;
        int idxZ = relP.z() / h;
        bool inInterval =
            (relP.x() >= 0 )&& (relP.y() >= 0 )&& (relP.z() >= 0 )&& (idxX < nx - 1 )&& (idxY < ny - 1 )&& (idxZ < nz - 1);
        if (inInterval)
        {
            int blfIdx = idxX +     idxY * nx +     idxZ * nx * ny;
            int brfIdx = idxX + 1 + idxY * nx +     idxZ * nx * ny;
            int blbIdx = idxX +     (idxY + 1)* nx +idxZ * nx * ny;
            int brbIdx = idxX + 1 + (idxY + 1)* nx +idxZ * nx * ny;
            int tlfIdx = idxX +     idxY * nx +     (idxZ + 1) * nx * ny;
            int trfIdx = idxX + 1 + idxY * nx +     (idxZ + 1) * nx * ny;
            int tlbIdx = idxX +     (idxY +1)* nx + (idxZ + 1) * nx * ny;
            int trbIdx = idxX + 1 + (idxY +1)* nx + (idxZ + 1) * nx * ny;
            
            double residualX = relP.x() - idxX * h;
            double residualY = relP.y() - idxY * h;
            double residualZ = relP.z() - idxZ * h;
            // tri-linear Interpolation
            double wx = residualX / h, wy = residualY / h, wz = residualZ / h;

            double wblf = (1. - wx) * (1. - wy) * (1. - wz); // weight of bottom-left-front
            double wbrf =  wx       * (1. - wy) * (1. - wz); 
            double wblb = (1. - wx) *    wy     * (1. - wz);
            double wbrb =  wx       *    wy     * (1. - wz); 
            double wtlf = (1. - wx) * (1. - wy) *  wz; 
            double wtrf =  wx       * (1. - wy) *  wz;
            double wtlb = (1. - wx) *    wy     *  wz;
            double wtrb =  wx       *    wy     *  wz; 
            
            tripleLists.push_back({i, blfIdx, wblf});
            tripleLists.push_back({i, brfIdx, wbrf});
            tripleLists.push_back({i, blbIdx, wblb});
            tripleLists.push_back({i, brbIdx, wbrb});
            tripleLists.push_back({i, tlfIdx, wtlf});
            tripleLists.push_back({i, trfIdx, wtrf});
            tripleLists.push_back({i, tlbIdx, wtlb});
            tripleLists.push_back({i, trbIdx, wtrb});
            
        }
    }
    W.setFromTriplets(tripleLists.begin(), tripleLists.end());
}

void PoissonRec::fd_partial_derivative(const int nx, const int ny, const int nz, const double h, const int dir,
                                       Eigen::SparseMatrix<double>& D) 
{
    std::vector<Eigen::Triplet<double>> tripletList;
    // determind which dim subtract 1 form partial derivative direction (dir) 
    // all the logic about direction are solve here
    int sub_nx = 0, sub_ny = 0, sub_nz = 0;
    int nextIndexDiff = 0; // the index difference of origin grid

    
    switch (dir) {
    case 0:
        sub_nx = 1;
        nextIndexDiff = 1;
        break;
    case 1:
        sub_ny = 1;
        nextIndexDiff = nx;
        break;
    case 2:
        sub_nz = 1;
        nextIndexDiff = nx * ny;
        break;
    default:
        std::cout << "unknown partial derivate direction!" << std::endl;
    }
    for (int i = 0; i < nx-sub_nx; ++i)
    {
        for (int j = 0; j < ny-sub_ny; ++j)
        {
            for (int k = 0; k < nz-sub_nz; ++k)
            {
                int idx = i + j * (nx - sub_nx) + k * (nx - sub_nx) * (ny - sub_ny);
                int firstIdx = i + j * nx + k * nx * ny;
                int nextIdx = firstIdx + nextIndexDiff;
                tripletList.push_back({idx, firstIdx, -1./h});
                tripletList.push_back({idx, nextIdx, 1./h});
            }
        }
    }
    D.setFromTriplets(tripletList.begin(), tripletList.end());
}
