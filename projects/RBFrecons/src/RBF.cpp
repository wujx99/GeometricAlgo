#include "RBF.h"
#include "geometrycentral/numerical/linear_solvers.h"



RBFrecons::RBFrecons(PointCloud* c, PointPositionGeometry* g) {
    addOffSurfacePoints(c, g);
    computeCoeff();

}

double RBFrecons::slove(Vector3 point) {
    double ret = 0.;
    for (size_t i = 0; i < coeff.size(); ++i)
    {
        ret += coeff[i].second * norm(point - coeff[i].first);
    }
    ret += constCoeff;
    ret += dot(polyCoeff, point);
    return ret;
}

RBFrecons::~RBFrecons() {
    
}

void RBFrecons::addOffSurfacePoints(PointCloud* c, PointPositionGeometry* g,
                                    double dis, size_t ratio){
    
    Vector<PointValPair> allPairs(3 * c->nPoints());
    g->requireNormals();
    assert(g->normals.size() == g->positions.size());
    for (Point p : c->points()) {
        size_t idx = p.getIndex();

        allPairs[idx].first = g->positions[idx];
        allPairs[idx].second = 0;
        allPairs[c->nPoints() + idx].first = g->positions[idx] + dis * g->normals[idx];
        allPairs[c->nPoints() +idx].second = dis;
        allPairs[2 * c->nPoints() + idx].first = g->positions[idx] - dis * g->normals[idx];
        allPairs[2 * c->nPoints() + idx].second = -dis;
        
    }
    // because of the dim of matrix is so big, so we select some of them 
    size_t ratioSize = 3 * c->nPoints() / ratio;
    
    allPointValPairs.resize(ratioSize);
    for (size_t i = 0; i < ratioSize; ++i)
    {
        allPointValPairs[i] = allPairs[i * ratio];
    }
}

void RBFrecons::computeCoeff() 
{
    
	size_t size = allPointValPairs.size();
    // construct linear equation
    DenseMatrix<double> A = DenseMatrix<double>::Zero(size, size);
    DenseMatrix<double> P = DenseMatrix<double>::Zero(size, 4);
    Vector<double> b = Vector<double>::Zero(size+4);

    for (size_t i = 0; i <  size; ++i)
    {
        Vector3 pos = allPointValPairs[i].first;
        for (size_t j = 0; j <  size; ++j)
        {
            A(i, j) = norm(pos - allPointValPairs[j].first);
        }
       
        b[i] = allPointValPairs[i].second;


        P(i, 0) = 1;
        P(i, 1) = pos.x;
        P(i, 2) = pos.y;
        P(i, 3) = pos.z;
    }


    DenseMatrix<double> B(size + 4, size + 4);
    B << A, P, P.transpose(), DenseMatrix<double>::Zero(4, 4);
    
    Vector<double> slo = B.householderQr().solve(b);
    
    coeff.resize(size);
    for (size_t i = 0; i < size; ++i)
    {
        coeff[i].first = allPointValPairs[i].first;
        coeff[i].second = slo[i];
    }
    constCoeff = slo[size];
    polyCoeff.x = slo[size+1];
    polyCoeff.y = slo[size+2];
    polyCoeff.z = slo[size+3];

}
