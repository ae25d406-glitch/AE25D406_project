#include "Preconditioners.h"
#include <stdexcept>
#include <algorithm>

using namespace std;

// M = I
struct NonePC : Preconditioner
{
    void apply(const Vec& r, Vec& z) override { z=r; }
};

// M = D
struct JacobiPC : Preconditioner
{
    const Stencil& A;
    explicit JacobiPC(const Stencil& A_):A(A_){}
    void apply(const Vec& r, Vec& z) override
    {
        for(size_t k=0;k<r.size();k++) z[k]=r[k]/A.aP[k];
    }
};

// M = (D+L) D^{-1} (D+U)  : one forward + one backward Gauss-Seidel sweep from z = 0
struct SGSPC : Preconditioner
{
    const Stencil& A;
    explicit SGSPC(const Stencil& A_):A(A_){}
    void apply(const Vec& r, Vec& z) override
    {
        fill(z.begin(),z.end(),0.0);
        A.sweep(z,r,1.0,true);
        A.sweep(z,r,1.0,false);
    }
};

// ILU(0): same algorithm as ilu0
struct ILU0PC : Preconditioner
{
    int n;
    vector<vector<pair<int,double>>> row;   // sorted (col, value) per row
    vector<int> diag;

    explicit ILU0PC(const Stencil& A)
    {
        n=A.n();// Ni*Nj
        row.resize(n); diag.resize(n);// one empty list per row
        int Ni=A.Ni,Nj=A.Nj;
        for(int I=0;I<Ni;I++)
        {
            int Iw=(I==0)?Ni-1:I-1, Ie=(I+1==Ni)?0:I+1;//check conditio for periodic
            for(int J=0;J<Nj;J++)
            {
                int k=I*Nj+J;
                auto& R=row[k];
                R.push_back({k,A.aP[k]});
                R.push_back({Iw*Nj+J,-A.aW[k]});
                R.push_back({Ie*Nj+J,-A.aE[k]});
                if(J>0)    R.push_back({k-1,-A.aS[k]});
                if(J<Nj-1) R.push_back({k+1,-A.aN[k]});
                sort(R.begin(),R.end());
            }
        }
        auto find=[&](int r,int c)->double*{
            for(auto& e:row[r]) if(e.first==c) return &e.second;
            return nullptr;
        };
        for(int i=0;i<n;i++)
            for(size_t q=0;q<row[i].size();q++) if(row[i][q].first==i) diag[i]=(int)q;

        for(int i=0;i<n;i++)                       // IKJ variant restricted to the pattern
        {
            for(size_t kk=0;kk<row[i].size();kk++)
            {
                int k=row[i][kk].first;
                if(k>=i) break;
                double& lik=row[i][kk].second;
                lik/=row[k][diag[k]].second;
                for(size_t jj=kk+1;jj<row[i].size();jj++)
                {
                    double* ukj=find(k,row[i][jj].first);
                    if(ukj) row[i][jj].second-=lik*(*ukj);
                }
            }
        }
    }
    void apply(const Vec& r, Vec& z) override
    {
        for(int i=0;i<n;i++)                       // L y = r
        {
            double s=r[i];
            for(auto& e:row[i]){ if(e.first>=i) break; s-=e.second*z[e.first]; }
            z[i]=s;
        }
        for(int i=n-1;i>=0;i--)                    // U z = y
        {
            double s=z[i];
            for(auto& e:row[i]) if(e.first>i) s-=e.second*z[e.first];
            z[i]=s/row[i][diag[i]].second;
        }
    }
};

// one symmetric V-cycle from zero initial guess
struct MGPC : Preconditioner
{
    StencilMultigrid mg;
    explicit MGPC(const Stencil& A):mg(A,2,2,false){}//2 pre-smoothing sweeps, 2 post-smoothing sweeps;
    void apply(const Vec& r, Vec& z) override
    {
        fill(z.begin(),z.end(),0.0);
        mg.vCycle(r,z);
    }
};

unique_ptr<Preconditioner> makePreconditioner(const string& name, const Stencil& A)
{
    if(name=="none")   return make_unique<NonePC>();
    if(name=="jacobi") return make_unique<JacobiPC>(A);
    if(name=="sgs")    return make_unique<SGSPC>(A);
    if(name=="ilu0")   return make_unique<ILU0PC>(A);
    if(name=="mg")     return make_unique<MGPC>(A);
    throw runtime_error("unknown preconditioner "+name);
}
