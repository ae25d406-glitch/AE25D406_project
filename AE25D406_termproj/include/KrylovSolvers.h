#ifndef KRYLOVSOLVERS_H
#define KRYLOVSOLVERS_H

// Krylov subspace solvers for the 5-point stencil  
//  CG and BiCGSTAB 
// Convergence: ||b - A x|| / ||b|| < tol

#include <string>
#include "Stencil.h"
#include "Preconditioners.h"

struct SolveInfo { int iters = 0; double relres = 0; };

SolveInfo cgSolve      (const Stencil& A, const Vec& b, Vec& x, Preconditioner& M, double tol, int maxit);
SolveInfo bicgstabSolve(const Stencil& A, const Vec& b, Vec& x, Preconditioner& M, double tol, int maxit);

// solver = cg | bicgstab | mg (multigrid alone);  prec = none|jacobi|sgs|ilu0|mg
SolveInfo solveSystem(const std::string& solver, const std::string& prec,
                      const Stencil& A, const Vec& b, Vec& x,
                      double tol, int maxit);

#endif
