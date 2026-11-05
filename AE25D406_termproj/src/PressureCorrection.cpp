
#include "PressureCorrection.h"
#include <cmath>
using namespace std;

PressureEquation newPressureEquation(int nx,int ny)
{
	PressureEquation eq;
	Field zero(nx,vector<double>(ny,0.0));
	eq.aW=eq.aE=eq.aS=eq.aN=eq.aP=eq.b=zero;
	return eq;
}


// when  omega =1 (gauss seidal) / SOR
static void sweep(const PressureEquation& eq,Field& x,double omega)
{
	int nx=x.size(),ny=x[0].size();

	for(int I=0;I<nx;I++)
		for(int J=0;J<ny;J++)
		{
			double sum=eq.b[I][J];
			if(I>0)sum+=eq.aW[I][J]*x[I-1][J];
			if(I<nx-1)sum+=eq.aE[I][J]*x[I+1][J];
			if(J>0)sum+=eq.aS[I][J]*x[I][J-1];
			if(J<ny-1)sum+=eq.aN[I][J]*x[I][J+1];

			double xNew=sum/eq.aP[I][J];
			x[I][J]+=omega*(xNew-x[I][J]);
		}
}

// residual r=b+sum(a_nd x_nd)-aP x
static Field residual(const PressureEquation& eq,const Field& x)
{
	int nx=x.size(),ny=x[0].size();
	Field r(nx,vector<double>(ny,0.0));

	for(int I=0;I<nx;I++)
		for(int J=0;J<ny;J++)
		{
			double sum=eq.b[I][J]-eq.aP[I][J]*x[I][J];
			if(I>0)sum+=eq.aW[I][J]*x[I-1][J];
			if(I<nx-1)sum+=eq.aE[I][J]*x[I+1][J];
			if(J>0)sum+=eq.aS[I][J]*x[I][J-1];
			if(J<ny-1)sum+=eq.aN[I][J]*x[I][J+1];
			r[I][J]=sum;
		}
	return r;
}

static double norm(const Field& r)
{
	double s=0.0;
	for(auto& col:r)for(double v:col)s+=v*v;
	return sqrt(s);
}

static void removeMean(Field& x)
{
	double mean=0.0;
	int n=0;
	for(auto& col:x)for(double v:col){mean+=v;n++;}
	mean/=n;
	for(auto& col:x)for(double& v:col)v-=mean;
}


//SOR 

int solveSOR(const PressureEquation& eq,Field& x,double omega,double relTol,int maxSweeps)
{
	double r0=norm(residual(eq,x));
	if(r0==0.0)return 0;

	int sweeps=0;
	while(sweeps<maxSweeps)
	{
		for(int k=0;k<10;k++)sweep(eq,x,omega);
		sweeps+=10;
		if(norm(residual(eq,x))<relTol*r0)break;
	}
	removeMean(x);
	return sweeps;
}


//multigrid
static PressureEquation coarsen(const PressureEquation& f)
{
	int nx=f.aP.size()/2,ny=f.aP[0].size()/2;
	PressureEquation c=newPressureEquation(nx,ny);

	for(int I=0;I<nx;I++)
		for(int J=0;J<ny;J++)
		{
			int i0=2*I,j0=2*J;

			c.aW[I][J]=f.aW[i0][j0]+f.aW[i0][j0+1];
			c.aE[I][J]=f.aE[i0+1][j0]+f.aE[i0+1][j0+1];
			c.aS[I][J]=f.aS[i0][j0]+f.aS[i0+1][j0];
			c.aN[I][J]=f.aN[i0][j0+1]+f.aN[i0+1][j0+1];
			c.aP[I][J]=c.aW[I][J]+c.aE[I][J]+c.aS[I][J]+c.aN[I][J];
		}
	return c;
}

static void vCycle(vector<PressureEquation>& level,int L,Field& x)
{
	PressureEquation& eq=level[L];
	int nx=x.size(),ny=x[0].size();

	if(L==(int)level.size()-1)
	{
		for(int k=0;k<50;k++)sweep(eq,x,1.0);
		return;
	}

	sweep(eq,x,1.0);
	sweep(eq,x,1.0);

	Field r=residual(eq,x);
	PressureEquation& ceq=level[L+1];
	for(int I=0;I<nx/2;I++)
		for(int J=0;J<ny/2;J++)
			ceq.b[I][J]=r[2*I][2*J]+r[2*I+1][2*J]+r[2*I][2*J+1]+r[2*I+1][2*J+1];

	Field e(nx/2,vector<double>(ny/2,0.0));
	vCycle(level,L+1,e);

	Field c(nx,vector<double>(ny,0.0));
	for(int I=0;I<nx;I++)
		for(int J=0;J<ny;J++)
			c[I][J]=e[I/2][J/2];

	double cr=0.0,cAc=0.0;
	for(int I=0;I<nx;I++)
		for(int J=0;J<ny;J++)
		{
			double Ac=eq.aP[I][J]*c[I][J];
			if(I>0)Ac-=eq.aW[I][J]*c[I-1][J];
			if(I<nx-1)Ac-=eq.aE[I][J]*c[I+1][J];
			if(J>0)Ac-=eq.aS[I][J]*c[I][J-1];
			if(J<ny-1)Ac-=eq.aN[I][J]*c[I][J+1];
			cr+=c[I][J]*r[I][J];
			cAc+=c[I][J]*Ac;
		}
	double s=(cAc>0.0)?cr/cAc:1.0;

	for(int I=0;I<nx;I++)
		for(int J=0;J<ny;J++)
			x[I][J]+=s*c[I][J];

	sweep(eq,x,1.0);
	sweep(eq,x,1.0);
}

int solveMultigrid(const PressureEquation& eq,Field& x,double relTol,int maxCycles)
{
	double r0=norm(residual(eq,x));
	if(r0==0.0)return 0;

	vector<PressureEquation>level;
	level.push_back(eq);
	while(level.back().aP.size()%2==0&&level.back().aP[0].size()%2==0&&
		level.back().aP.size()>=4&&level.back().aP[0].size()>=4)
		level.push_back(coarsen(level.back()));

	int cycles=0;
	while(cycles<maxCycles)
	{
		vCycle(level,0,x);
		cycles++;
		if(norm(residual(eq,x))<relTol*r0)break;
	}
	removeMean(x);
	return cycles;
}