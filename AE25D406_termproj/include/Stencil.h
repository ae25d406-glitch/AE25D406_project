#ifndef STENCIL_H
#define STENCIL_H

// 5-point stencil equation on the periodic O-grid    
//
// Same convention as PressureEquation in PressureCorrection.h:
//      aP*x_P = aW*x_W + aE*x_E + aS*x_S + aN*x_N + b
// aS = 0 at J = 0 and aN = 0 at J = Nj-1 (boundaries go into aP and b).

#include <vector>

typedef std::vector<double> Vec;
struct Stencil
{
    int Ni = 0, Nj = 0;
    Vec aW, aE, aS, aN, aP, b;

    void resize(int Ni_,int Nj_);
    void zero();//initialise value of the coeff
    int n() const { return Ni*Nj; }

    // y = A x   with  (A x)_P = aP x_P - sum a_nb x_nb
    void matvec(const Vec& x, Vec& y) const;
    // r = b - A x   (uses the vector rhs, not the member b)
    void residual(const Vec& x, const Vec& rhs, Vec& r) const;

    // one Gauss-Seidel / SOR sweep (forward or backward)
    void sweep(Vec& x, const Vec& rhs, double omega, bool forward) const;
};

double dot(const Vec& a, const Vec& b);//CG
double norm2(const Vec& a);

#endif
