#include "Multigrid.h"
#include <cmath>
#include <algorithm>

static constexpr double PI = 3.14159265358979323846;
//contruct
MultigridLap::MultigridLap(int N)
{
	//n+1 grid points
	N_=N;
	u_.resize(N+1,0.0);
	applyBoundaryCondi(u_);
}

//grid index
int MultigridLap::idx( int i) const
{
	return i;
}

//analytical solution u(x)=sin(pi*x), for -u''=f with f=pi^2*sin(pi*x)
double MultigridLap:: analytical(double x)
{
	return std::sin(PI*x);
}

//f(x)=pi^2*sin(pi*x)
static double rhsFunc(double x)
{
	return PI*PI*std::sin(PI*x);
}

//bc u(0)=0 and u(L)=0
void MultigridLap::applyBoundaryCondi(std::vector<double>&u)
{
	u[0]=0.0;
	u[N_]=0.0;
}

//smoothing and fd u[i] = (u[i-1] + u[i+1] + h^2)/2
void MultigridLap::smooth(std::vector<double>& u,const std::vector<double>& f, int n ,int sweeps, Smoother type)
{
	double h=1.0/n;
	for (int s=0; s<sweeps; s++)
	{
		if(type==GAUSS_SEIDEL)
		{
			for (int i=1;i<n;i++)
			{
				u[i]=0.5*(u[i-1]+u[i+1]+h*h*f[i]);
			}
		}

		else
		{
			//jacobi
			const double omega=2.0/3.0;
			std::vector<double>old=u;
			for (int i=1;i<n;i++)
			{
				double gs_update=0.5*(old[i-1]+old[i+1]+h*h*f[i]);
				u[i]=(1.0-omega)*old[i]+omega*gs_update;
			}
		}
		//bc 
		u[0]=0.0;
		u[n]=0.0;
	}
}

//residual Au=f, r=f-Au, Au=(-u[i-1] + 2u[i] - u[i+1])/h^2

void MultigridLap::computeResi(const std::vector<double>& u, const std::vector<double>& f, std::vector<double>& r, int n)
{
	double h=1.0/n;
	r.assign(n+1,0.0);
	for (int i=1;i<n;i++)
	{
		double Au= (-u[i-1]-u[i+1]+2*u[i])/(h*h);
		r[i]=f[i]-Au;
	}
	r[0]=0.0;
	r[n]=0.0;
}

//resriction 
//fine to course
//coarse[i]=fine[2i]

void MultigridLap::restrictResi(const std::vector<double>& fine, std::vector<double>& coarse, int nFine)
{
	int ncoarse=nFine/2;
	coarse.assign(ncoarse+1,0.0);
	for (int i=1;i<ncoarse;i++)
	{
		coarse[i]=fine[2*i];
	}
	coarse[0]=0.0;
	coarse[ncoarse]=0.0;

}
//prolongation  coarse to fine
//linear interpolation

void MultigridLap::prolongCorrec(const std::vector<double>& coarse, std::vector<double>& fine, int nFine)
{
	int ncoarse=nFine/2;
	//coarse pts
	for (int i=0;i<=ncoarse;i++)
	{
		fine[2*i]+=coarse[i]; //from fine to coarse which reduce by 2 F2i=Ci 
	}
	//pts btw coarse pts
	for (int i=0;i<ncoarse;i++)
	{
		fine[2*i+1]+=0.5*(coarse[i]+coarse[i+1]);//eg ;F3=C1+C2/2 this is linear interpolation
	}
}

//V cycle

void MultigridLap::vcycle(std::vector<double>& u,const std::vector<double>& f, int n, int level, int numLevels, int nu1,int nu2,Smoother smoother)
{
//coarse grid smooth iteration
	if (level==numLevels -1)
	{
		smooth(u,f,n,20,smoother);
		return;
	}

	// 1. pre smoothing
	smooth(u,f,n,nu1,smoother);
	
	// 2. calculate residual
	std::vector<double>residual;
	computeResi(u,f,residual,n);

	// 3. restrict residual to coarse grid
	int ncoarse=n/2;
	std::vector<double>coarseResidual;
	restrictResi(residual,coarseResidual,n);
	
	// 4. solve error eqn in coarse gid
	//Ae=residual

	std::vector<double> error( ncoarse+1,0.0);
	vcycle(error,coarseResidual,ncoarse,level+1,numLevels,nu1,nu2,smoother);

	// 5. Prolongate coarse correction
	std::vector<double> correction(n+1,0.0);
	prolongCorrec(error,correction,n);

	// 6. Correct fine grid solution
	for (int i=1;i<n;i++)
	{
		u[i]+=correction[i];
	}
	
	// 7. Post smoothing
	smooth(u,f,n,nu2,smoother);
}

//Residual norm

double MultigridLap::norm( const std::vector<double>&r) const
{
	double sum=0.0;
	for(double value: r)
	{
		sum+=value*value;
	}
	return std::sqrt(sum);
}

//solve multigrid
std::vector<double>MultigridLap::solve(int numLevels,Smoother smoother,int nu1,int nu2,double tolerance,int maxcycle)
{
	//start from zero
	u_.assign(N_+1,0.0);
	applyBoundaryCondi(u_);
	
	//Rhs f=1
	std::vector<double> f(N_+1,1.0);
	for (int i=0; i<=N_ ;i++)
	{
		f[i]=rhsFunc(getX(i));
	}
	f[0]=0.0;
	f[N_]=0.0;
	
	//store residual data
	std::vector<double>residualHis;
	

	//multigrid cycle
	for (int cycle=0;cycle<maxcycle;cycle++)
	{
		vcycle(u_,f,N_,0,numLevels,nu1,nu2,smoother);
	
		//calculate residual
		std::vector<double>residual;
		computeResi(u_,f,residual,N_);
		double rnorm=norm(residual);
		residualHis.push_back(rnorm);

		//convergence check
		if (rnorm<tolerance)
		{
			break;
		}
	}
	return residualHis;
}
// solve single grid


std::vector<double>MultigridLap::solveSingleGrid(Smoother smoother,double tolerance,int maxIter)
{
	//start from zero
	u_.assign(N_+1,0.0);
	applyBoundaryCondi(u_);
	
	//Rhs f=1
	std::vector<double> f(N_+1,0.0);
	for (int i=0; i<=N_ ;i++)
	{
		f[i]=rhsFunc(getX(i));
	}
	f[0]=0.0;
	f[N_]=0.0;
	
	//store residual data
	std::vector<double>residualHis;
	
	for (int it=0; it<maxIter;it++)
	{
		smooth(u_,f,N_,1,smoother);
	
		//calculate residual
		std::vector<double>residual;
		computeResi(u_,f,residual,N_);
		double rnorm=norm(residual);
		residualHis.push_back(rnorm);

		//convergence check
		if (rnorm<tolerance)
		{
			break;
		}
	}
	return residualHis;
}


//get solution

double MultigridLap::getU(int i) const
{
	return u_[idx(i)];
}
//get N
int MultigridLap::getN() const
{
	return N_;
}

//get x coordinate

double MultigridLap::getX(int i) const
{
	return static_cast<double>(i)/N_;
}

//L2 norm
double MultigridLap::computeL2Error() const
{
	double sum=0.0;
	for (int i=0; i<=N_;i++)
	{
		double x=getX(i);
		double error=getU(i)-analytical(x);
		sum+=error*error;
	}
	return std::sqrt(sum/(N_+1));

}

//L infinity error
double MultigridLap::computeLinfError() const
{
        double Maxerr=0.0;
        for (int i=0; i<=N_;i++)
        {
                double x=getX(i);
                double error= std::abs(getU(i)-analytical(x));
                Maxerr=std::max(Maxerr,error);
        }
        return Maxerr;
}


