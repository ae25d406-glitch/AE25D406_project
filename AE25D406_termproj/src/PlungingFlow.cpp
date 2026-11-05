#include "PlungingFlow.h"
#include "KrylovSolvers.h"
#include <cmath>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <algorithm>

using namespace std;
using clk=chrono::steady_clock;

// One PISO time step (collocated grid, AB2/CN, central differencing)
// dV = cell volume, F = face volume flux u.S, D = nu*A/dx (diffusion conductance), a_P, a_nb, b = coefficients.
// 1. inflow u_inf = (1,-hdot) at t^{n+1}; mark outer faces as inflow/outflow
// 2. convection C^n = sum_f F_f u_f  with u_f = (u_P+u_nb)/2                       (eq. 3.3, 3.4)
// 3. predictor: a_P u*_P = sum a_nb u*_nb + b - dV (dp/dx)^n                        (eq. 3.9-3.11, 8.1)
//    a_nb = D/2, a_P = dV/dt + sum a_nb,  b = dV/dt u^n + D/2 (old diffusion) - (3/2 C^n - 1/2 C^{n-1})
// 4. corrector (x2): uHat = (sum a_nb u_nb + b)/a_P, d = dV/a_P                    (eq. 6.1)
//    pressure eq. a_P p_P = sum a_nb p_nb + b,  a_nb = dbar*A/dx,  b = Fhat_in - Fhat_out   (eq. 7.3-7.5)
//    then F = Fhat - a(p_nb - p_P) and u = uHat - d grad p                          (eq. 6.3, 6.1)

PlungingFlow::PlungingFlow(const Mesh& mesh_,const OGrid& G_,const PlungeSettings& s)
    : mesh(mesh_), G(G_), set(s) // store the mesh, grid and settings
{
    N=G.n();// number of cells
    nu=1.0/set.Re;// ν = 1/Re
    for(auto* f : {&u,&v,&p,&un,&vn,&Cu,&Cv,&CuOld,&CvOld,&F_e,&F_n,&b_u,&b_v,&uHat,&vHat,&d})
        f->assign(N,0.0);//size N filled with 0
    inflow.assign(G.Ni,0);
    momEq.resize(G.Ni,G.Nj);// momentum matrix
    pEqn.resize(G.Ni,G.Nj);// pressure matrix

    setInflow(0.0);//far-field velocity at t = 0
    for(int k=0;k<N;k++){ u[k]=uInf; v[k]=vInf; }//free stream velocity
    initialProjection();
}

// plunge h = h0 R(t) sin(kt).  R = amplitude ramp: R = (1 - cos(pi t/Tr))/2 for t < Tr, then 1 (Tr = rampTime).
// rampTime = 0 (default): R = 1, the plunge starts at full amplitude.
double PlungingFlow::ramp(double t,double& dR) const
{
    dR=0.0;
    if(set.rampTime<=0.0 || t>=set.rampTime) return 1.0;
    double a=M_PI/set.rampTime;// a = π/T_r
    dR=0.5*a*sin(a*t);// dR/dt = ½ (π/T_r) sin(πt/T_r)
    return 0.5*(1.0-cos(a*t));// R = ½ (1 − cos(πt/T_r))
}
double PlungingFlow::h(double t)  const { double dR; return set.h0*ramp(t,dR)*sin(set.k*t); }//h = h₀ R sin kt
double PlungingFlow::hd(double t) const { double dR; double R=ramp(t,dR); return set.h0*(dR*sin(set.k*t)+R*set.k*cos(set.k*t)); }// derivative of above h//ḣ = h₀ (Ṙ sin kt + R k cos kt)

//inflow seen by the foil: (U cos a, U sin a - hdot)
void PlungingFlow::setInflow(double t)
{
    double a=set.aoaDeg*M_PI/180.0;// mean angle of attack,degrees -radians
    if(t<set.perturbTime) a+=set.perturbDeg*M_PI/180.0;//small tilt at the start only (cylinder)
    uInf=cos(a);// u∞ = cos a
    vInf=sin(a)-hd(t);// v∞ = sin a − ḣ 
    int J=G.Nj-1;
    for(int I=0;I<G.Ni;I++)
    {
        int k=G.id(I,J);
        inflow[I]=(uInf*G.Snx[k]+vInf*G.Sny[k]<0.0);// u∞·S_n< 0 air inflow 
    }
}


// Green-Gauss cell gradient of pressure
// face value: interior average; wall and inflow: p_P (zero gradient); outflow: 0
void PlungingFlow::gradP(const vector<double>& pf,vector<double>& gx,vector<double>& gy) const
{
    fill(gx.begin(),gx.end(),0.0); fill(gy.begin(),gy.end(),0.0);
    for(int I=0;I<G.Ni;I++)
        for(int J=0;J<G.Nj;J++)
        {
            int k=G.id(I,J), kE=G.id(G.ip(I),J);
            double f=0.5*(pf[k]+pf[kE]);// p_e = (p_P + p_E)/2
            gx[k]+=f*G.Sex[k]; gy[k]+=f*G.Sey[k];// gradient of pressure at x and y 
            gx[kE]-=f*G.Sex[k]; gy[kE]-=f*G.Sey[k];// for east neghbour the face vector points in to P cell so neg
            if(J<G.Nj-1)// not the top row 
            {
                int kN=k+1; // north neighbour
                f=0.5*(pf[k]+pf[kN]); // p_n = (p_P + p_N)/2
                gx[k]+=f*G.Snx[k]; gy[k]+=f*G.Sny[k];
                gx[kN]-=f*G.Snx[k]; gy[kN]-=f*G.Sny[k];
            }
            else
            {// top row
                f=inflow[I]? pf[k] : 0.0; // BC:(pressure gradient =0) inflow: p_f = p_P ,(pressure is zero)outflow: p_f = 0
                gx[k]+=f*G.Snx[k]; gy[k]+=f*G.Sny[k];
            }
            if(J==0){ gx[k]+=pf[k]*G.Swx[I]; gy[k]+=pf[k]*G.Swy[I]; }//pressure grad is 0 and p_wall=p_P
        }
    for(int k=0;k<N;k++){ gx[k]/=G.dV[k]; gy[k]/=G.dV[k]; }
}

// C^n = sum over faces of F_f * u_f  (central differencing)
//F conv is computed in pressure corrector 
void PlungingFlow::convection()
{
    fill(Cu.begin(),Cu.end(),0.0); fill(Cv.begin(),Cv.end(),0.0);
    const double beta=set.blend;
    auto face=[beta](double F,double aP,double aN){
        double up=(F>=0.0)? aP : aN;
        return (1.0-beta)*0.5*(aP+aN)+beta*up;
    };
    for(int I=0;I<G.Ni;I++)
        for(int J=0;J<G.Nj;J++)
        {// east face of P
            int k=G.id(I,J), kE=G.id(G.ip(I),J);// this cell P ,east neighbour E 
            double F=F_e[k];// F_e = u.S_e
            double uf=face(F,u[k],u[kE]), vf=face(F,v[k],v[kE]);// u_e,v_e
            Cu[k]+=F*uf; Cv[k]+=F*vf; Cu[kE]-=F*uf; Cv[kE]-=F*vf;// P= + F_e u_e ,E= - F_e u_e
            
            //north face of P
            if(J<G.Nj-1)
            {
                int kN=k+1;
                F=F_n[k]; // F_n = u.S_n
                uf=face(F,u[k],u[kN]); vf=face(F,v[k],v[kN]);// u_n,v_n
                Cu[k]+=F*uf; Cv[k]+=F*vf; Cu[kN]-=F*uf; Cv[kN]-=F*vf;//(leaves) P= + F_n u_n,(enters)  N= - F_n u_n
            }
            else// top row
            {
                F=F_n[k];
                uf=inflow[I]? uInf : u[k];// inflow u_f = u_inf
                vf=inflow[I]? vInf : v[k];// outflow u_f = u_P (du/dn = 0)
                Cu[k]+=F*uf; Cv[k]+=F*vf;
            }
        }
}

// momentum matrix (time + 1/2 diffusion, symmetric) and RHS without pressure
void PlungingFlow::buildMomentum(double dt,double bodyForceY)
{
    momEq.zero();
    bool ab2=haveOld;

    for(int I=0;I<G.Ni;I++)
    {
        int kWi=G.im(I), kEi=G.ip(I);//wraps-I=0 -> Ni-1
        for(int J=0;J<G.Nj;J++)
        {
            int k=G.id(I,J), kW=G.id(kWi,J), kE=G.id(kEi,J);//cell P,W,E
            double aP0=G.dV[k]/dt;                                   // a_P^0 = dV/dt

            // diffusion conductances D = nu*A/dx = nu*g; wall: distance centre->wall
            double De=nu*G.gE[k], Dw=nu*G.gE[kW];//east and west face

            // south face (J>0)= north face of the cell below   and south face (J=0): the WALL
            double Dn=nu*G.gN[k], Ds=(J>0)? nu*G.gN[k-1] : nu*G.gW[I];

            // Crank-Nicolson: half implicit (a_nb = D/2), half explicit with old u (goes into b)
            double aE=0.5*De, aW=0.5*Dw, aN=0.5*Dn, aS=0.5*Ds;
            momEq.aE[k]=aE; momEq.aW[k]=aW;
            double oldDiffU=aE*(u[kE]-u[k])+aW*(u[kW]-u[k]);     // explicit half of diffusion
            double oldDiffV=aE*(v[kE]-v[k])+aW*(v[kW]-v[k]);
            double bcU=0, bcV=0;

            //south nb
            if(J>0){ momEq.aS[k]=aS; oldDiffU+=aS*(u[k-1]-u[k]); oldDiffV+=aS*(v[k-1]-v[k]); }
            else   { momEq.aP[k]+=aS; oldDiffU+=aS*(0.0-u[k]); oldDiffV+=aS*(0.0-v[k]); }        // J=0

            //north
            if(J<G.Nj-1){ momEq.aN[k]=aN; oldDiffU+=aN*(u[k+1]-u[k]); oldDiffV+=aN*(v[k+1]-v[k]); }
            
            else if(inflow[I])//top                                            // u = u_inf
            {
                momEq.aP[k]+=aN;//implici thalf
                bcU=aN*uInf; bcV=aN*vInf;
                oldDiffU+=aN*(uInf-u[k]); oldDiffV+=aN*(vInf-v[k]);
            }                                                                      // outflow: zero gradient

            // a_P = dV/dt + a_E + a_W + a_S + a_N
            momEq.aP[k]+=aP0+momEq.aE[k]+momEq.aW[k]+momEq.aS[k]+momEq.aN[k];

//conv AB2(explicit)
            double convU= ab2? 1.5*Cu[k]-0.5*CuOld[k] : Cu[k];
            double convV= ab2? 1.5*Cv[k]-0.5*CvOld[k] : Cv[k];

//RHS b (known terms)
            b_u[k]=aP0*u[k]+oldDiffU+bcU-convU;
            b_v[k]=aP0*v[k]+oldDiffV+bcV-convV-G.dV[k]*bodyForceY;// foil accelaration
        }
    }
}


// pseudo-velocity 
//uhat velocity without pressure grad
//momentum:  a_P u_P = sum a_nb u_nb + b - dV dp/dx
//  u_P = (sum a_nb u_nb + b)/a_P  -  (dV/a_P) dp/dx  =  uHat_P - d_P dp/dx
void PlungingFlow::computeUhat()// reads momEq, b_u, b_v, u, v
{
    for(int I=0;I<G.Ni;I++)
    {
        int Iw=G.im(I), Ie=G.ip(I);// west / east column (periodic)
        for(int J=0;J<G.Nj;J++)
        {
            int k=G.id(I,J), kW=G.id(Iw,J), kE=G.id(Ie,J);
            // neighbour sum + source, WITHOUT the pressure gradient 
            //su = b + a_W u_W + a_E u_E (+ a_S u_S + a_N u_N)
            double su=b_u[k]+momEq.aW[k]*u[kW]+momEq.aE[k]*u[kE];
            double sv=b_v[k]+momEq.aW[k]*v[kW]+momEq.aE[k]*v[kE];

            if(J>0){ su+=momEq.aS[k]*u[k-1]; sv+=momEq.aS[k]*v[k-1]; }
            if(J<G.Nj-1){ su+=momEq.aN[k]*u[k+1]; sv+=momEq.aN[k]*v[k+1]; }
            uHat[k]=su/momEq.aP[k]; // uHat_P = (sum a_nb u_nb + b)/a_P   velocity without pressure
            vHat[k]=sv/momEq.aP[k]; 
            d[k]=G.dV[k]/momEq.aP[k];// d_P = dV/a_P
        }
    }
}

//pressure corrector
// pressure equation from continuity with Rhie-Chow fluxes, then correct fluxes and u
int PlungingFlow::pressureCorrector()
{
    pEqn.zero();//pressure matrix
    vector<double> Fhat_e(N,0.0), Fhat_n(N,0.0);//pressure free flux
    for(int I=0;I<G.Ni;I++)
        for(int J=0;J<G.Nj;J++)
        {
            int k=G.id(I,J), kE=G.id(G.ip(I),J);//cell P, E

            // east face: a_face = dbar*A/dx = dbar*g,  Fhat = uHat_face . S 
            double aFace=0.5*(d[k]+d[kE])*G.gE[k]; 
            // Fhat_e = uHat_e . S_e
            double FhatFace=0.5*(uHat[k]+uHat[kE])*G.Sex[k]+0.5*(vHat[k]+vHat[kE])*G.Sey[k];
            pEqn.aE[k]=aFace; 
            pEqn.aW[kE]=aFace; 
            pEqn.aP[k]+=aFace; 
            pEqn.aP[kE]+=aFace;
            pEqn.b[k]-=FhatFace; 
            pEqn.b[kE]+=FhatFace; 
            Fhat_e[k]=FhatFace;   // b = Fhat_in - Fhat_out
            
            //north
            if(J<G.Nj-1)
            {
                int kN=k+1;
                aFace=0.5*(d[k]+d[kN])*G.gN[k];
                FhatFace=0.5*(uHat[k]+uHat[kN])*G.Snx[k]+0.5*(vHat[k]+vHat[kN])*G.Sny[k];
                pEqn.aN[k]=aFace; pEqn.aS[kN]=aFace; pEqn.aP[k]+=aFace; pEqn.aP[kN]+=aFace;
                pEqn.b[k]-=FhatFace; pEqn.b[kN]+=FhatFace; Fhat_n[k]=FhatFace;
            }
            //top row inflow face
            else if(inflow[I])
            {
                Fhat_n[k]=uInf*G.Snx[k]+vInf*G.Sny[k];
                pEqn.b[k]-=Fhat_n[k];
            }
            else  //top row outflow face                                      // outflow: p = 0 on the face
            {
                aFace=d[k]*G.gN[k];
                Fhat_n[k]=uHat[k]*G.Snx[k]+vHat[k]*G.Sny[k];
                pEqn.aP[k]+=aFace; pEqn.b[k]-=Fhat_n[k];
            }
        }


    //solve apPp
    //CG+ one multigrid V cycle as preconditioner
    SolveInfo si=solveSystem(set.pSolver,set.pPrec,pEqn,pEqn.b,p,set.pTol,set.maxIt);

    // CORRECT THE FACE FLUXES:  F = Fhat - a_f (p_nb - p_P)
    // these F satisfy  sum_f F = 0  in every cell 
    for(int I=0;I<G.Ni;I++)
        for(int J=0;J<G.Nj;J++)
        {
            int k=G.id(I,J), kE=G.id(G.ip(I),J);
            F_e[k]=Fhat_e[k]-0.5*(d[k]+d[kE])*G.gE[k]*(p[kE]-p[k]);
            if(J<G.Nj-1)      F_n[k]=Fhat_n[k]-0.5*(d[k]+d[k+1])*G.gN[k]*(p[k+1]-p[k]);
            else if(inflow[I])F_n[k]=Fhat_n[k];
            else              F_n[k]=Fhat_n[k]+d[k]*G.gN[k]*p[k];
        }
    vector<double> gx(N),gy(N);
    gradP(p,gx,gy);// Green–Gauss gradient
    for(int k=0;k<N;k++){ u[k]=uHat[k]-d[k]*gx[k]; v[k]=vHat[k]-d[k]*gy[k]; }// u_P = uHat_P - d_P (dp/dx)_P
    return si.iters;
}

// make the uniform start field divergence-free (uHat = u, d = 1)
void PlungingFlow::initialProjection()
{
    uHat=u; vHat=v; fill(d.begin(),d.end(),1.0);
    fill(p.begin(),p.end(),0.0);
    pressureCorrector();
    fill(p.begin(),p.end(),0.0);
}

//FORCE
// force on the body: pressure + wall shear (S_w points into the body)
void PlungingFlow::forces(double& Fx,double& Fy) const
{
    Fx=Fy=0;
    for(int I=0;I<G.Ni;I++)
    {
        int k=G.id(I,0);
        Fx+=p[k]*G.Swx[I]+nu*G.gW[I]*u[k];
        Fy+=p[k]*G.Swy[I]+nu*G.gW[I]*v[k];
    }
}

void PlungingFlow::vorticity(vector<double>& w) const// writes w[k] = omega in every cell
{
    w.assign(N,0.0);
    for(int I=0;I<G.Ni;I++)
        for(int J=0;J<G.Nj;J++)
        {
            int k=G.id(I,J), kE=G.id(G.ip(I),J);
            double uf=0.5*(u[k]+u[kE]), vf=0.5*(v[k]+v[kE]);// u_e = (u_P + u_E)/2 
            double c=vf*G.Sex[k]-uf*G.Sey[k];// v_e S_x - u_e S_y  = velocity along the face x face length
            w[k]+=c; w[kE]-=c;// cell P= + c   (S_e points out of P);cell E=- c   (same face, opposite direction for E)
            
            if(J<G.Nj-1)//north
            {
                uf=0.5*(u[k]+u[k+1]); vf=0.5*(v[k]+v[k+1]);
                c=vf*G.Snx[k]-uf*G.Sny[k];
                w[k]+=c; w[k+1]-=c;// P: + c,N: - c
            }
            else// top row
            {
                uf=inflow[I]?uInf:u[k]; vf=inflow[I]?vInf:v[k];// inflow: u_inf ; outflow: u_P (zero gradient)
                w[k]+=vf*G.Snx[k]-uf*G.Sny[k];
            }
        }
    for(int k=0;k<N;k++) w[k]/=G.dV[k];
}

//ONE TIME STEP  t^n -> t^{n+1}   (predictor + PISO correctors)
void PlungingFlow::step()
{
    double dt=set.dt, tNew=time+dt;// time step
    double bodyForceY=(hd(tNew)-hd(time))/dt;      // hddot over the step: exact for the far field

    un=u; vn=v;
    convection();                                   // C^n with fluxes and inflow of step n
    setInflow(tNew);
    buildMomentum(dt,bodyForceY);

    // predictor: a_P u* - sum a_nb u*_nb = b - dV grad p^n   (eq. 8.1)
    vector<double> gx(N),gy(N),ru(N),rv(N);
    gradP(p,gx,gy);
    for(int k=0;k<N;k++){ ru[k]=b_u[k]-G.dV[k]*gx[k]; rv[k]=b_v[k]-G.dV[k]*gy[k]; }// RHS_u = b_u - dV dp/dx
    
     // solve momEq u* = RHS_u  and  momEq v* = RHS_v 
    itU=solveSystem(set.mSolver,set.mPrec,momEq,ru,u,set.mTol,set.maxIt).iters;
    itV=solveSystem(set.mSolver,set.mPrec,momEq,rv,v,set.mTol,set.maxIt).iters;

    // PISO correctors
    for(int c=0;c<set.nCorr;c++)
    {
        computeUhat();// uHat = (sum a_nb u_nb + b)/a_P from  u (u* then u^1), d = dV/a_P
        int it=pressureCorrector();// pressure equation gives p ; F = Fhat - a(p_nb - p_P) ; u = uHat - d grad p
        if(c<2) itP[c]=it;
    }

    CuOld=Cu; CvOld=Cv; haveOld=true;
    time=tNew; steps++;
}

void PlungingFlow::writeVTK(const string& file) const
{
    ofstream f(file);
    f<<"# vtk DataFile Version 3.0\nplunging foil\nASCII\nDATASET STRUCTURED_GRID\n";
    f<<"DIMENSIONS "<<G.Ni+1<<" "<<G.Nj+1<<" 1\nPOINTS "<<mesh.getNumberOfNodes()<<" float\n";
    for(auto& nd:mesh.nodes) f<<nd.x<<" "<<nd.y<<" 0\n";
    vector<double> w; vorticity(w);
    f<<"CELL_DATA "<<N<<"\n";
    auto sc=[&](const char* nm,const vector<double>& a){
        f<<"SCALARS "<<nm<<" float 1\nLOOKUP_TABLE default\n";
        for(int J=0;J<G.Nj;J++) for(int I=0;I<G.Ni;I++) f<<(float)a[G.id(I,J)]<<"\n";
    };
    sc("pressure",p); sc("vorticity",w);
    f<<"VECTORS velocity float\n";
    for(int J=0;J<G.Nj;J++) for(int I=0;I<G.Ni;I++){ int k=G.id(I,J); f<<(float)u[k]<<" "<<(float)v[k]<<" 0\n"; }
}

void PlungingFlow::solve()
{
    filesystem::create_directories(set.outDir);
    fForces.open(set.outDir+"/forces.csv");
    fStats.open(set.outDir+"/stats.csv");
    fForces<<"t,h,hdot,CD,CL,CP\n";
    fStats<<"step,t,itU,itV,itP1,itP2,maxDiv,maxCFL,cpuStep\n";
    if(set.vtkEvery>0) writeVTK(set.outDir+"/field_0000.vtk");

    auto start=clk::now();
    int nOut=1;
    while(time<set.tEnd-1e-9)
    {
        auto t0=clk::now();
        step();
        double cpu=chrono::duration<double>(clk::now()-t0).count();

        double Fx,Fy; forces(Fx,Fy);
        double CD=2*Fx, CL=2*Fy, CP=-2*Fy*hd(time);

        double maxDiv=0, maxCFL=0;
        for(int I=0;I<G.Ni;I++)
            for(int J=0;J<G.Nj;J++)
            {
                int k=G.id(I,J), kW=G.id(G.im(I),J);
                double net=F_e[k]-F_e[kW]+F_n[k]-(J>0?F_n[k-1]:0.0);
                double sumAbs=fabs(F_e[k])+fabs(F_e[kW])+fabs(F_n[k])+(J>0?fabs(F_n[k-1]):0.0);
                maxDiv=max(maxDiv,fabs(net)/G.dV[k]);
                maxCFL=max(maxCFL,0.5*sumAbs*set.dt/G.dV[k]);
            }

        if(!isfinite(CD) || !isfinite(CL)){ cout<<"*** diverged at step "<<steps<<", t = "<<time<<endl; break; }

        fForces<<setprecision(10)<<time<<","<<h(time)<<","<<hd(time)<<","<<CD<<","<<CL<<","<<CP<<"\n";
        fStats<<steps<<","<<time<<","<<itU<<","<<itV<<","<<itP[0]<<","<<itP[1]<<","<<maxDiv<<","<<maxCFL<<","<<cpu<<"\n";

        if(set.printEvery>0 && steps%set.printEvery==0)
        {
            cout<<fixed<<setprecision(3)<<"step "<<setw(6)<<steps<<"  t = "<<setw(8)<<time
                <<"  CD = "<<setw(7)<<CD<<"  CL = "<<setw(7)<<CL
                <<"  it(u,v,p1,p2) = "<<itU<<","<<itV<<","<<itP[0]<<","<<itP[1]
                <<scientific<<setprecision(1)<<"  div = "<<maxDiv
                <<fixed<<setprecision(2)<<"  CFL = "<<maxCFL<<endl;
            fForces.flush(); fStats.flush();
        }
        if(set.vtkEvery>0 && steps%set.vtkEvery==0)
        {
            char nm[32]; snprintf(nm,sizeof nm,"/field_%04d.vtk",nOut++);
            writeVTK(set.outDir+nm);
        }
    }
    double total=chrono::duration<double>(clk::now()-start).count();
    cout<<"done: "<<steps<<" steps, CPU "<<total<<" s"<<endl;
}
