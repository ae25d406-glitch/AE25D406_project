#ifndef PRECONDITIONERS_H
#define PRECONDITIONERS_H

// Preconditioners for the 5-point stencil  
// Same list as get_M() in assignment_7.cpp (none, jacobi, sgs, ilu0),
// rewritten for the sparse periodic stencil, plus one multigrid V-cycle.
// apply(r, z) computes z = M^{-1} r.

#include <memory>
#include <string>
#include "Stencil.h"
#include "StencilMultigrid.h"

struct Preconditioner
{
    virtual ~Preconditioner() {}
    virtual void apply(const Vec& r, Vec& z) = 0;
};

std::unique_ptr<Preconditioner> makePreconditioner(const std::string& name, const Stencil& A);

#endif
