#include "KrylovSolvers.h"
#include <cmath>
#include <stdexcept>

using namespace std;

//conjugate gradient (preconditioned)
//alpha =(r'z)/(p'A p)
//x=x+alpha p,  r_new =r-alpha A p,  z = M^-1 r_new
//beta =(r_new'z_new)/(r'z),  p=z+beta p
SolveInfo cgSolve(const Stencil& A,const Vec& b,Vec& x,Preconditioner& M,double tol,int maxit)
{
    int n=A.n();// number of cells (unknowns)
    SolveInfo info;
    double nb=norm2(b); if(nb==0) nb=1;// |b|,used to make the residual relative; avoid ÷0
    Vec r(n),z(n),p(n),Ap(n);// residual, preconditioned residual, direction, A·p
    A.residual(x,b,r);// r₀ = b − A x₀
    info.relres=norm2(r)/nb;// |r₀|/|b|
    if(info.relres<tol) return info;
    M.apply(r,z); p=z;// z₀ = M⁻¹ r₀ ;  first direction p₀ = z₀
    double rz=dot(r,z);
    for(int k=1;k<=maxit;k++)
    {
        A.matvec(p,Ap);// Ap(the only multiplication by A per iteration)
        double den=dot(p,Ap);// p·A p
        if(fabs(den)<1e-300) throw runtime_error("CG fails");// p·Ap = 0 would mean ÷0
        double alpha=rz/den;// α = (r·z)/(p·A p)
        for(int i=0;i<n;i++){ x[i]+=alpha*p[i]; r[i]-=alpha*Ap[i]; }// move x; update r
        info.iters=k; //save iteration count
        info.relres=norm2(r)/nb;// |r|/|b|
        if(info.relres<tol) break;
        M.apply(r,z); // z_new = M⁻¹ r_new preconditioner
        double rzNew=dot(r,z);// r_new·z_new
        double beta=rzNew/rz; rz=rzNew;// β = (r_new·z_new)/(r·z) how far old direction in new dorection 
        for(int i=0;i<n;i++) p[i]=z[i]+beta*p[i];// p_new = z_new + β p
    }
    return info;
}

//BiCGSTAB (right preconditioned)
//rho= r_hat' r,  beta=(rho/rho_old)(alpha/omega),  p=r+beta(p-omega v)
//v=A M^-1 p,  alpha=rho/(r_hat' v),  s=r-alpha v
//t=A M^-1 s,  omega=(t's)/(t't),  x=x+alpha M^-1 p+omega M^-1 s,  r=s-omega t
SolveInfo bicgstabSolve(const Stencil& A,const Vec& b,Vec& x,Preconditioner& M,double tol,int maxit)
{
    int n=A.n();
    SolveInfo info;
    double nb=norm2(b); if(nb==0) nb=1;
    Vec r(n),rhat(n),p(n,0.0),v(n,0.0),s(n),t(n),ph(n),sh(n);
    A.residual(x,b,r); rhat=r;// r₀ = b − A x₀ ;  r̂ = r₀ 
    info.relres=norm2(r)/nb;
    if(info.relres<tol) return info;
    double rhoOld=1,alpha=1,omega=1;// start values so the first β = ρ 
    for(int k=1;k<=maxit;k++)
    {
        double rho=dot(rhat,r);// ρ = r̂·r
        if(fabs(rho)<1e-300) break;
        double beta=(rho/rhoOld)*(alpha/omega);// β = (ρ/ρ_old)(α/ω)
        for(int i=0;i<n;i++) p[i]=r[i]+beta*(p[i]-omega*v[i]);// p = r + β(p − ω v)new direction
        M.apply(p,ph); A.matvec(ph,v);// p̂ = M⁻¹p ;  v = A p̂
        alpha=rho/dot(rhat,v);// α = ρ/(r̂·v)
        for(int i=0;i<n;i++) s[i]=r[i]-alpha*v[i];// s = r − α v residual after half-step
        info.iters=k;
        if(norm2(s)/nb<tol)//if s is already small enough:
        {
            for(int i=0;i<n;i++) x[i]+=alpha*ph[i];//x = x + α p̂   and stop
            info.relres=norm2(s)/nb; break;
        }
        //half step2 
        M.apply(s,sh); A.matvec(sh,t);// ŝ = M⁻¹s ;  t = A ŝ
        omega=dot(t,s)/dot(t,t);//ω = (t·s)/(t·t): the step that makes |s − ωt| smallest
        //ph-preconditioned direction sh-preconditioned leftover residual
        for(int i=0;i<n;i++){ x[i]+=alpha*ph[i]+omega*sh[i]; r[i]=s[i]-omega*t[i]; }// x = x + α p̂ + ω ŝ ,// r = s − ω t
        rhoOld=rho;
        info.relres=norm2(r)/nb;
        if(info.relres<tol || omega==0) break;
    }
    return info;
}

SolveInfo solveSystem(const string& solver,const string& prec,const Stencil& A,const Vec& b,Vec& x,
                      double tol,int maxit)
{
    if(solver=="mg")
    {
        StencilMultigrid mg(A,2,2,false);
        SolveInfo info;
        info.iters=mg.solve(b,x,tol,maxit,info.relres);
        return info;
    }
    auto M=makePreconditioner(prec,A);
    if(solver=="cg")       return cgSolve(A,b,x,*M,tol,maxit);
    if(solver=="bicgstab") return bicgstabSolve(A,b,x,*M,tol,maxit);
    throw runtime_error("unknown solver "+solver);
}
