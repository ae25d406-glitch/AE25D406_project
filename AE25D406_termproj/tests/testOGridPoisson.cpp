//g++ -O2 -std=gnu++17 -Iinclude tests/testOGridPoisson.cpp src/*.cpp -o tests/testOGridPoisson
// Laplace equation around a circle with exact solution (potential flow):
//   phi = x (1 + a^2/r^2),  a = 0.5,  d(phi)/dn = 0 on the circle
// Checks: grid areas, 2nd-order accuracy, and every solver/preconditioner.
#include <iostream>
#include <iomanip>
#include <cmath>
#include "OGrid.h"
#include "KrylovSolvers.h"

using namespace std;

double exact(double x,double y){ double r2=x*x+y*y; return x*(1.0+0.25/r2); }// φ = x(1 + a²/r²), a² = 0.25

int main()
{
    double Rfar=10.0, prevErr=0;// outer radius 10, // prevErr = error on the previous (coarser) grid
    cout<<setw(10)<<"grid"<<setw(14)<<"area err"<<setw(14)<<"L2 error"<<setw(8)<<"order";
    const char* combos[][2]={{"cg","none"},{"cg","jacobi"},{"cg","sgs"},{"cg","ilu0"},{"cg","mg"},
                             {"mg","-"},{"bicgstab","ilu0"},{"bicgstab","mg"}};// solver + preconditioner 
    for(auto& c:combos) cout<<setw(13)<<(string(c[0])+"+"+c[1]);
    cout<<"\n";
    for(int Ni : {32,64,128,256})// 32, 64, 128, 256 cells around the circle
    {
        int Nj=Ni/2;// 16, 32, 64, 128 cells outward
        Mesh mesh=makeConformalMesh("circle",0,Ni,Nj,Rfar);// the SAME grid code
        OGrid G; G.build(mesh);// ΔV, S, g, wall faces
        double area=0; for(double v:G.dV) area+=v;// sum of all cell areas ΔV
        double areaExact=acos(-1.0)*(Rfar*Rfar-0.25); // exact ring area π(R_far² − a²)
        
        // polygon edges are straight, so compare with the polygon area
        Stencil A; A.resize(Ni,Nj);
        for(int I=0;I<Ni;I++) for(int J=0;J<Nj;J++)
        {
            int k=G.id(I,J);
            A.aE[k]=G.gE[k]; A.aW[k]=G.gE[G.id(G.im(I),J)];
            if(J>0) A.aS[k]=G.gN[k-1];
            if(J<Nj-1) A.aN[k]=G.gN[k];
            A.aP[k]=A.aE[k]+A.aW[k]+A.aS[k]+A.aN[k];
            if(J==Nj-1)  // Dirichlet at the outer face centre
            {
                double fx=G.xc[k]+0.5*G.Snx[k]*0, fy=G.yc[k];
                // face centre: move from cell centre along the radial direction by half a cell
                double r=sqrt(G.xc[k]*G.xc[k]+G.yc[k]*G.yc[k]);
                double rf=Rfar*cos(acos(-1.0)/Ni);  // outer polygon edge midpoint radius
                fx=G.xc[k]*rf/r; fy=G.yc[k]*rf/r;
                A.aP[k]+=G.gN[k]; A.b[k]+=G.gN[k]*exact(fx,fy);
            }
        }
        double err=0, vol=0;
        Vec x0(A.n(),0.0), x(A.n());
        x=x0; SolveInfo s=solveSystem("cg","mg",A,A.b,x,1e-12,5000);
        for(int k=0;k<A.n();k++)
        {
            if(sqrt(G.xc[k]*G.xc[k]+G.yc[k]*G.yc[k])>2.0) continue;   // error near the body
            double e=x[k]-exact(G.xc[k],G.yc[k]); err+=e*e*G.dV[k]; vol+=G.dV[k];// error at the cell centre: computed − exact
        }
        err=sqrt(err/vol);// L2 error = √( Σ e² ΔV / Σ ΔV )
        cout<<setw(6)<<Ni<<"x"<<setw(3)<<Nj<<scientific<<setprecision(2)<<setw(14)<<fabs(area-areaExact)/areaExact
            <<setw(14)<<err<<fixed<<setprecision(2)<<setw(8)<<(prevErr>0?log(prevErr/err)/log(2.0):0.0);
        prevErr=err;
        for(auto& c:combos)
        {
            x=x0;
            SolveInfo si=solveSystem(c[0],c[1],A,A.b,x,1e-8,20000);
            cout<<setw(13)<<si.iters;
        }
        cout<<"\n";
    }
}
