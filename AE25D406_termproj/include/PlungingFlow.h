#ifndef PLUNGINGFLOW_H
#define PLUNGINGFLOW_H

// Incompressible flow past a plunging body on the conformal O-grid
//   - collocated variables (u, v, p at cell centres) with Rhie-Chow face fluxes
//   - central differencing in space
//   - Adams-Bashforth 2 (convection, explicit) + Crank-Nicolson (diffusion, implicit)
//   - foil-fixed frame: the grid never moves, the inflow is (U, -hdot(t)) and the
//     frame acceleration enters the v-equation as the body force -hddot
//   - linear systems solved with the Krylov / multigrid solvers of KrylovSolvers.h
// Non-dimensional: chord (or diameter) = 1, U = 1, rho = 1, nu = 1/Re.

#include <string>
#include <vector>
#include <fstream>
#include "Mesh.h"
#include "OGrid.h"
#include "Stencil.h"

struct PlungeSettings
{
    double Re=100.0;
    double h0=0.0;            // plunge amplitude
    double k=0.0;             // reduced frequency, h = h0 sin(k t)
    double aoaDeg=0.0;        // mean angle of attack
    double perturbTime=0.0;   // tilt the inflow for t < perturbTime (starts shedding)
    double perturbDeg=0.0;
    double rampTime=0.0;      // plunge amplitude ramped smoothly from 0 over this time (0 = no ramp)
    double dt=0.01;
    double tEnd=10.0;
    int    nCorr=2;           // PISO correctors
    double blend=0.0;         // upwind fraction in the convective face value (0 = pure central)
    std::string pSolver="cg", pPrec="mg";         // pressure system
    std::string mSolver="cg", mPrec="jacobi";     // momentum system
    double pTol=1e-8, mTol=1e-10;
    int    maxIt=2000;
    int    printEvery=100;
    int    vtkEvery=0;        // 0 = no field output during the run
    std::string outDir="results";
};

class PlungingFlow
{
public:
    PlungingFlow(const Mesh& mesh, const OGrid& grid, const PlungeSettings& set);// set up
    void solve();
    void writeVTK(const std::string& file) const;// write u, v, p, vorticity for ParaView

    double time=0.0;//current time t
    long   steps=0;//number of steps taken

private:
    const Mesh& mesh;
    const OGrid& G;
    PlungeSettings set;
    int N;
    double nu;

    std::vector<double> u, v, p, un, vn;//un, vn-u, v from the previous step
    std::vector<double> Cu, Cv, CuOld, CvOld;// convection C = sum F u_f at steps n and n-1 (AB2)
    std::vector<double> F_e, F_n;// face volume fluxes F = u.S (east, north face of each cell)
    std::vector<double> b_u, b_v;// source b of the momentum equation (without pressure)
    std::vector<double> uHat, vHat, d;// pseudo-velocity and d = dV/a_P; u = û − d ∇p
    std::vector<char>   inflow;
    Stencil momEq, pEqn;// coefficients a_P, a_E, a_W, a_N, a_S, b of momentum / pressure eq.
    double uInf=1.0, vInf=0.0;// far-field velocity seen by the foil (1, -hdot)
    bool haveOld=false; // C^{n−1} (false on the very first step)

    std::ofstream fForces, fStats;
    int itU=0, itV=0, itP[2]={0,0};

    double h(double t)  const;
    double hd(double t) const;
    double ramp(double t, double& dR) const;

    void setInflow(double t);
    void gradP(const std::vector<double>& pf, std::vector<double>& gx, std::vector<double>& gy) const;
    void convection();
    void buildMomentum(double dt, double bodyForceY);
    void computeUhat();          // uHat and d
    int  pressureCorrector();
    void initialProjection();
    void forces(double& Fx, double& Fy) const;
    void vorticity(std::vector<double>& w) const;
    void step();
};

#endif
