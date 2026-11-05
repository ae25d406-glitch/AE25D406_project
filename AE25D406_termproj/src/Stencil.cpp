#include "Stencil.h"
#include <cmath>
#include <algorithm>

void Stencil::resize(int Ni_,int Nj_)
{
    Ni=Ni_; Nj=Nj_;
    int N=Ni*Nj;// no of cells
    for(Vec* v : {&aW,&aE,&aS,&aN,&aP,&b}) v->assign(N,0.0); // asssign discritised eqn in each cell
}

void Stencil::zero()
{
    for(Vec* v : {&aW,&aE,&aS,&aN,&aP,&b}) std::fill(v->begin(),v->end(),0.0);//initialize to zero 
}

void Stencil::matvec(const Vec& x, Vec& y) const
{
    for(int I=0;I<Ni;I++)
    {
        int Iw=(I==0)?Ni-1:I-1, Ie=(I+1==Ni)?0:I+1;
        for(int J=0;J<Nj;J++)
        {
            int k=I*Nj+J;
            double s=aP[k]*x[k]-aW[k]*x[Iw*Nj+J]-aE[k]*x[Ie*Nj+J];
            if(J>0)    s-=aS[k]*x[k-1];// bc at bottom
            if(J<Nj-1) s-=aN[k]*x[k+1];// bc at top 
            y[k]=s;//store y_k = a_P·x_P − Σ a_nb·x_nb
        }
    }
}

void Stencil::residual(const Vec& x, const Vec& rhs, Vec& r) const
{
    matvec(x,r);//computes A·x and stores it in r and rhs=b
    for(size_t k=0;k<r.size();k++) r[k]=rhs[k]-r[k];//resi in each cell eqn 
}

// Gauss-Seidel / SOR sweep (as in PressureCorrection.cpp, periodic in I)
void Stencil::sweep(Vec& x, const Vec& rhs, double omega, bool forward) const
{
    for(int s=0;s<Ni;s++)
    {
        int I=forward? s : Ni-1-s;//f forward is true, I = s; otherwise I = Ni−1−s.
        int Iw=(I==0)?Ni-1:I-1, Ie=(I+1==Ni)?0:I+1;//Iw= the I-index of the west neighbour. Normally I−1, but if I = 0 it wraps to Ni−1 = 3.(periodic bc)
        for(int t=0;t<Nj;t++)
        {
            int J=forward? t : Nj-1-t;
            int k=I*Nj+J;
            double sum=rhs[k]+aW[k]*x[Iw*Nj+J]+aE[k]*x[Ie*Nj+J];
            if(J>0)    sum+=aS[k]*x[k-1];
            if(J<Nj-1) sum+=aN[k]*x[k+1];
            double xNew=sum/aP[k];//solve x_P = (rhs_P + Σ a_nb·x_nb) / a_P
            x[k]+=omega*(xNew-x[k]);
        }
    }
}

double dot(const Vec& a, const Vec& b)
{
    double s=0;
    for(size_t i=0;i<a.size();i++) s+=a[i]*b[i];
    return s;
}

double norm2(const Vec& a){ return std::sqrt(dot(a,a)); }//multigrid relative residual to check less the tolerance to stop Vcycle 

