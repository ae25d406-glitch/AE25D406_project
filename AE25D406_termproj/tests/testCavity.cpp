//g++ -std=gnu++17 -O2 -I include tests/testCavity.cpp src/*.cpp -o tests/testCavity
// ./tests/testCavity


#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <iomanip>
#include "Test.h"
#include "Point.h"
#include "Line.h"
#include "Surface.h"
#include "Mesh.h"
#include "PressureCorrection.h"
#include "CavityFlow.h"

//reference value
const double PSI_REF_32=-0.10209843;
const double PSI_GHIA=-0.103423;

Mesh makeMesh(int N)
{
	Point p0(0,0.0,0.0),p1(1,1.0,0.0),p2(2,1.0,1.0),p3(3,0.0,1.0);
	Line l0(0,&p0,&p1),l1(1,&p1,&p2),l2(2,&p2,&p3),l3(3,&p3,&p0);
	Surface s(0);
	s.addBoundary(&l0);s.addBoundary(&l1);s.addBoundary(&l2);s.addBoundary(&l3);
	Mesh m;
	m.generateCartesian(s,N,N);
	return m;
}

FlowSettings settings(const std::string& alg,double Utop,double Ubottom,double Re)
{
	FlowSettings s;
	s.Re=Re;s.Utop=Utop;s.Ubottom=Ubottom;s.algorithm=alg;
	s.tol=1e-9;s.printEvery=0;s.maxIter=20000;
	if(alg=="PISO"){s.dt=0.02;}
	if(alg=="PIMPLE"){s.dt=0.5;s.nOuter=3;}
	return s;
}

//stream function
double psiMin(const Field& u,int N)
{
	double h=1.0/N,best=0.0;
	for(int i=0;i<=N;i++)
	{
		double psi=0.0;
		for(int J=0;J<N;J++){psi+=u[i][J]*h;best=std::min(best,psi);}
	}
	return best;
}

double maxDivergence(const Field& u,const Field& v,int N)
{
	double h=1.0/N,m=0.0;
	for(int I=0;I<N;I++)
		for(int J=0;J<N;J++)
			m=std::max(m,std::fabs((u[I+1][J]-u[I][J])*h+(v[I][J+1]-v[I][J])*h));
	return m;
}

double interp(const std::vector<double>& x,const std::vector<double>& f,double xq)
{
	for(size_t k=1;k<x.size();k++)
		if(xq<=x[k]){double t=(xq-x[k-1])/(x[k]-x[k-1]);return (1-t)*f[k-1]+t*f[k];}
	return f.back();
}

int main()
{
	Test t;
	{
		int n=64;
		//pressure solver
		PressureEquation eq=newPressureEquation(n,n);
		std::mt19937 g(7);
		std::uniform_real_distribution<double> U(-1.0,1.0);
		double mean=0.0;
		for(int I=0;I<n;I++)
			for(int J=0;J<n;J++)
			{
				eq.aW[I][J]=(I>0)?1.0:0.0;eq.aE[I][J]=(I<n-1)?1.0:0.0;
				eq.aS[I][J]=(J>0)?1.0:0.0;eq.aN[I][J]=(J<n-1)?1.0:0.0;
				eq.aP[I][J]=eq.aW[I][J]+eq.aE[I][J]+eq.aS[I][J]+eq.aN[I][J];
				eq.b[I][J]=U(g);mean+=eq.b[I][J];
			}
		for(auto& col:eq.b)for(double& b:col)b-=mean/(n*n);

		Field xMG(n,std::vector<double>(n,0.0)),xSOR=xMG;
		solveMultigrid(eq,xMG,1e-10,500);
		solveSOR(eq,xSOR,1.9,1e-10,200000);
		double d=0.0;
		for(int I=0;I<n;I++)for(int J=0;J<n;J++)d=std::max(d,std::fabs(xMG[I][J]-xSOR[I][J]));
		t.expectTrue(d<1e-6,"Multigrid and SOR give the same p' solution");

		int cycles[2],sizes[2]={32,128};
		for(int k=0;k<2;k++)
		{
			int m=sizes[k];
			PressureEquation e=newPressureEquation(m,m);
			for(int I=0;I<m;I++)
				for(int J=0;J<m;J++)
				{
					e.aW[I][J]=(I>0)?1.0:0.0;e.aE[I][J]=(I<m-1)?1.0:0.0;
					e.aS[I][J]=(J>0)?1.0:0.0;e.aN[I][J]=(J<m-1)?1.0:0.0;
					e.aP[I][J]=e.aW[I][J]+e.aE[I][J]+e.aS[I][J]+e.aN[I][J];
					e.b[I][J]=std::sin(3.0*I/m)*std::cos(2.0*J/m);
				}
			double s=0.0;for(auto& col:e.b)for(double b:col)s+=b;
			for(auto& col:e.b)for(double& b:col)b-=s/(m*m);
			Field x(m,std::vector<double>(m,0.0));
			cycles[k]=solveMultigrid(e,x,1e-6,500);
		}
		std::cout<<"  multigrid cycles: 32x32 -> "<<cycles[0]<<", 128x128 -> "<<cycles[1]<<std::endl;
		t.expectTrue(std::abs(cycles[1]-cycles[0])<=3,"Multigrid cycle count does not grow with the mesh");
	}
	{
		// single lid :SIMPLE, PISO, PIMPLE
		const char* algs[3]={"SIMPLE","PISO","PIMPLE"};
		double psi[3];
		for(int a=0;a<3;a++)
		{
			Mesh mesh=makeMesh(32);
			CavityFlow flow(mesh,settings(algs[a],1.0,0.0,100.0));
			flow.solve();
			std::string name=algs[a];
			t.expectTrue(flow.converged,name+" converges (single lid, Re 100, 32x32)");
			psi[a]=psiMin(flow.getU(),32);
			std::cout<<std::setprecision(8)<<"  "<<name<<": psi_min = "<<psi[a]<<", steps = "<<flow.iterations<<std::endl;
			t.expectEqual(psi[a],PSI_REF_32,1e-5,name+" psi_min matches the stored reference");
			t.expectEqual(psi[a],PSI_GHIA,0.002,name+" psi_min within 0.002 of Ghia");
			t.expectTrue(maxDivergence(flow.getU(),flow.getV(),32)<1e-9,name+" satisfies discrete continuity");
		}
		t.expectEqual(psi[1],psi[0],1e-6,"PISO and SIMPLE reach the same steady solution");
		t.expectEqual(psi[2],psi[0],1e-6,"PIMPLE and SIMPLE reach the same steady solution");
	}
	{
		// Ghia re =100 n=64
		int N=64;
		Mesh mesh=makeMesh(N);
		CavityFlow flow(mesh,settings("SIMPLE",1.0,0.0,100.0));
		flow.solve();
		const Field& u=flow.getU();
		const Field& v=flow.getV();

		std::vector<double> y{0.0},uc{0.0},x{0.0},vc{0.0};
		for(int J=0;J<N;J++){y.push_back((J+0.5)/N);uc.push_back(u[N/2][J]);}
		y.push_back(1.0);uc.push_back(1.0);
		for(int I=0;I<N;I++){x.push_back((I+0.5)/N);vc.push_back(v[I][N/2]);}
		x.push_back(1.0);vc.push_back(0.0);

		double gy[5]={0.9531,0.7344,0.5000,0.2813,0.1016};
		double gu[5]={0.68717,0.00332,-0.20581,-0.15662,-0.06434};
		double gx[5]={0.9453,0.8047,0.5000,0.2344,0.0938};
		double gv[5]={-0.10313,-0.24533,0.05454,0.17527,0.12317};
		double eu=0.0,ev=0.0;
		for(int k=0;k<5;k++)
		{
			eu=std::max(eu,std::fabs(interp(y,uc,gy[k])-gu[k]));
			ev=std::max(ev,std::fabs(interp(x,vc,gx[k])-gv[k]));
		}
		std::cout<<"  Ghia check 64x64: max |du| = "<<eu<<", max |dv| = "<<ev<<std::endl;
		t.expectTrue(eu<0.01&&ev<0.01,"Re 100, 64x64: centre-line velocities within 0.01 of Ghia");
	}
	{
		//moving lid symmetric
		int N=32;
		//parralel lid 

		Mesh m1=makeMesh(N);
		CavityFlow par(m1,settings("SIMPLE",1.0,1.0,100.0));
		par.solve();
		const Field& u1=par.getU();
		double e1=0.0;
		for(int i=0;i<=N;i++)for(int J=0;J<N;J++)e1=std::max(e1,std::fabs(u1[i][J]-u1[i][N-1-J]));
		t.expectTrue(par.converged&&e1<1e-6,"Parallel lids: flow is a mirror image about y = 0.5");

		//antiparallel lids 
		Mesh m2=makeMesh(N);
		CavityFlow anti(m2,settings("SIMPLE",1.0,-1.0,100.0));
		anti.solve();
		const Field& u2=anti.getU();
		double e2=0.0;
		for(int i=0;i<=N;i++)for(int J=0;J<N;J++)e2=std::max(e2,std::fabs(u2[i][J]+u2[N-i][N-1-J]));
		t.expectTrue(anti.converged&&e2<1e-6,"Antiparallel lids: flow is point-symmetric about the centre");
		std::cout<<"  symmetry errors: parallel "<<e1<<", antiparallel "<<e2<<std::endl;
	}
	{
		Mesh m1=makeMesh(32);
		FlowSettings s1=settings("SIMPLE",1.0,0.0,100.0);
		s1.alphaP=1.0;s1.maxIter=2000;
		CavityFlow f1(m1,s1);
		f1.solve();
		t.expectTrue(!f1.converged,"SIMPLE without pressure under-relaxation (alphaP = 1) diverges");

		Mesh m2=makeMesh(64);
		FlowSettings s2=settings("SIMPLE",1.0,0.0,1000.0);
		s2.scheme="CENTRAL";s2.maxIter=2000;
		CavityFlow f2(m2,s2);
		f2.solve();
		t.expectTrue(!f2.converged,"CENTRAL scheme at Re 1000 (cell Peclet > 2) diverges");
	}

	t.summary();
	return 0;
}