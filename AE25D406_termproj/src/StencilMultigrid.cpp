#include "StencilMultigrid.h"
#include <cmath>
#include <algorithm>

using namespace std;

//coarse cell (Ic,Jc) = fine cells (2Ic..2Ic+1, 2Jc..2Jc+1)
Stencil StencilMultigrid::coarsen(const Stencil& f)
{
    Stencil c;//coarse matrix
    int Nic=f.Ni/2, Njc=f.Nj/2, Nj=f.Nj;//coarse size is half in each direction
    c.resize(Nic,Njc);
    auto F=[&](int I,int J){ return I*Nj+J; };//fine cell number (same k = I*Nj + J rule)

    for(int Ic=0;Ic<Nic;Ic++)
        for(int Jc=0;Jc<Njc;Jc++)//for each coarse cell (Ic,Jc)
        {
            int i0=2*Ic, j0=2*Jc;
            int a=F(i0,j0), b=F(i0+1,j0), cN=F(i0,j0+1), d=F(i0+1,j0+1);// coarse cell location
            int k=Ic*Njc+Jc;// coarse cell number

            // coarse neigbour coeff
            c.aW[k]=f.aW[a]+f.aW[cN];//west side of the block = west links of a and cN
            c.aE[k]=f.aE[b]+f.aE[d];// east side  = east links of b and d
            c.aS[k]=f.aS[a]+f.aS[b];// south side = south links of a and b
            c.aN[k]=f.aN[cN]+f.aN[d];// north side = north links of cN and d

            // sum of fine aP minus couplings that stay inside the coarse cell
            double internal = f.aE[a]+f.aW[b]+f.aE[cN]+f.aW[d]
                            + f.aN[a]+f.aS[cN]+f.aN[b]+f.aS[d];
            c.aP[k]=f.aP[a]+f.aP[b]+f.aP[cN]+f.aP[d]-internal;
        }
    return c;
}

StencilMultigrid::StencilMultigrid(const Stencil& A,int nPre_,int nPost_,bool scale_)
    : nPre(nPre_), nPost(nPost_), scale(scale_)// store the settings
{
    level.push_back(A);//level0 = the fine matrix
    while(level.back().Ni%2==0 && level.back().Nj%2==0 &&
          level.back().Ni>=4 && level.back().Nj>=4 && level.back().n()>256)
        level.push_back(coarsen(level.back()));

    for(auto& L:level){ xs.emplace_back(L.n(),0.0); bs.emplace_back(L.n(),0.0); rs.emplace_back(L.n(),0.0); }//creates x,band r vectors for every level and filled with zero

    // dense LU of the coarsest operator
    const Stencil& C=level.back();// the coarsest matrix
    nc=C.n();// 128;For 256*128 levels:256*128→128*64→64*32→32*16→1*×8=128 cells
    lu.assign(nc*nc,0.0); piv.resize(nc);
    Vec e(nc,0.0), col(nc);
    for(int j=0;j<nc;j++)
    {
        e[j]=1.0; C.matvec(e,col); e[j]=0.0;
        for(int i=0;i<nc;i++) lu[i*nc+j]=col[i];
    }
    for(int k=0;k<nc;k++)
    {
        int p=k; double mx=fabs(lu[k*nc+k]);
        for(int r=k+1;r<nc;r++) if(fabs(lu[r*nc+k])>mx){mx=fabs(lu[r*nc+k]);p=r;}
        piv[k]=p;
        if(p!=k) for(int q=0;q<nc;q++) swap(lu[k*nc+q],lu[p*nc+q]);
        double dkk=lu[k*nc+k]; if(fabs(dkk)<1e-300) dkk=1e-300;
        for(int r=k+1;r<nc;r++)
        {
            double m=lu[r*nc+k]/=dkk;
            for(int q=k+1;q<nc;q++) lu[r*nc+q]-=m*lu[k*nc+q];
        }
    }
}

void StencilMultigrid::coarseSolve()
{
    Vec& y=xs.back(); y=bs.back();
    for(int k=0;k<nc;k++)
    {
        swap(y[k],y[piv[k]]);
        for(int r=k+1;r<nc;r++) y[r]-=lu[r*nc+k]*y[k];
    }
    for(int k=nc-1;k>=0;k--)
    {
        double s=y[k];
        for(int q=k+1;q<nc;q++) s-=lu[k*nc+q]*y[q];
        y[k]=s/lu[k*nc+k];
    }
}

void StencilMultigrid::cycle(int L)
{
    if(L==(int)level.size()-1){ coarseSolve(); return; }// on the coarsest level solve excatly 

    const Stencil& A=level[L];// matrix on this level
    Vec& x=xs[L]; Vec& b=bs[L]; Vec& r=rs[L];

    for(int s=0;s<nPre;s++) A.sweep(x,b,1.0,true);       // pre-smoothing

    A.residual(x,b,r);//resi
    const Stencil& C=level[L+1];// next coarser matrix
    Vec& bc=bs[L+1];// its right-hand side
    for(int Ic=0;Ic<C.Ni;Ic++)                             // restriction (sum of 4)
        for(int Jc=0;Jc<C.Nj;Jc++)
        {
            int i0=2*Ic,j0=2*Jc;
            bc[Ic*C.Nj+Jc]=r[i0*A.Nj+j0]+r[(i0+1)*A.Nj+j0]+r[i0*A.Nj+j0+1]+r[(i0+1)*A.Nj+j0+1];
        }
    fill(xs[L+1].begin(),xs[L+1].end(),0.0);
    cycle(L+1);// solve A_c e = r_c  the correction e smooths,restrict and recurses again until corasest level is reached 

    const Vec& e=xs[L+1];//the coarse correction computed
    double s=1.0;//scale factor(1 = no scaling)
    if(scale)                                              
    {
        Vec c(A.n()),Ac(A.n());
        for(int I=0;I<A.Ni;I++) for(int J=0;J<A.Nj;J++) c[I*A.Nj+J]=e[(I/2)*C.Nj+J/2];
        A.matvec(c,Ac);
        double cAc=dot(c,Ac);
        if(cAc>0) s=dot(c,r)/cAc;
    }
    for(int I=0;I<A.Ni;I++)                                // prolongation
        for(int J=0;J<A.Nj;J++)
            x[I*A.Nj+J]+=s*e[(I/2)*C.Nj+J/2];

    for(int s2=0;s2<nPost;s2++) A.sweep(x,b,1.0,false);   // post-smoothing (backward)
}

void StencilMultigrid::vCycle(const Vec& rhs, Vec& x)
{
    bs[0]=rhs; xs[0]=x;
    cycle(0);
    x=xs[0];
}

int StencilMultigrid::solve(const Vec& rhs, Vec& x, double relTol, int maxCycles, double& relres)
{
    double nb=norm2(rhs); if(nb==0) nb=1;
    Vec r(rhs.size());
    level[0].residual(x,rhs,r);
    relres=norm2(r)/nb;// relative residual |r|/|b|
    int it=0;
    while(relres>relTol && it<maxCycles)//runs V-cycles one after another and stops when |r|/|b| < the tolerance
    {
        vCycle(rhs,x);
        level[0].residual(x,rhs,r);
        relres=norm2(r)/nb;
        it++;
    }
    return it;
}
