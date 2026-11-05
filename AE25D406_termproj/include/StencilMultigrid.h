#ifndef STENCILMULTIGRID_H
#define STENCILMULTIGRID_H

//solveMultigrid()/coarsen()/vCycle() in PressureCorrection.cpp:
// 2x2 cells are merged into one coarse cell 
//symmetric V-cycle (forward GS before, backward GS after) and a fixed
//coarse solve, so one V-cycle is a valid preconditioner for CG

#include "Stencil.h"

class StencilMultigrid
{
public:
//nPre-Gauss–Seidel sweeps before to course.
//nPost-sweeps after fine (that is after course to fine)
    StencilMultigrid(const Stencil& A, int nPre=2, int nPost=2, bool scaleCorrection=false);

    void vCycle(const Vec& rhs, Vec& x);                 // one cycle, x is updated
    int  solve(const Vec& rhs, Vec& x, double relTol, int maxCycles, double& relres);
    int  levels() const { return (int)level.size(); }

private:
    std::vector<Stencil> level;//xs[L] the unknown on level L. On the fine level this is the pressure; on the coarse levels it is the correction e.
    std::vector<Vec> xs, bs, rs;//solution, right-hand side, residual
    int nPre, nPost;
    bool scale;
    std::vector<double> lu; std::vector<int> piv; int nc = 0;   // dense LU of coarsest

    static Stencil coarsen(const Stencil& f);// build the coarse matrix from a fine one
    void cycle(int L);//the recursive V-cycle starting at level L
    void coarseSolve();//exact solve on the coarsest level
};

#endif
