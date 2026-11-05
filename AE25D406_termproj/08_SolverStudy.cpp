//g++ -O2 -std=gnu++17 08_SolverStudy.cpp src/*.cpp -Iinclude -o solverstudy
// ./solverstudy shape=ellipse sizes=64,128,256,512 Re=500 dt=0.002 out=results/solver_study
//
// Linear-solver study on the conformal O-grid:
//  (1) pressure system  sum g (p_P - p_N) = b  (SPD, Dirichlet on the outflow half)
//  (2) momentum system  [V/dt + 1/2 diffusion] u = b   (SPD, AB2/CN)
// For every grid size and solver/preconditioner pair: iterations and CPU time.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <cmath>
#include <filesystem>
#include "OGrid.h"
#include "KrylovSolvers.h"

using namespace std;
using clk=chrono::steady_clock;

int main(int argc,char* argv[])
{
	string shape="ellipse", sizes="64,128,256,512", outDir="results/solver_study";
	double thickness=0.12, Rfar=30, Re=500, dt=0.002, tol=1e-8, xiStretch=0.3;
	for(int a=1;a<argc;a++)
	{
		string arg=argv[a]; size_t eq=arg.find('=');
		string key=arg.substr(0,eq), val=arg.substr(eq+1);
		if(key=="shape")shape=val; else if(key=="sizes")sizes=val; else if(key=="out")outDir=val;
		else if(key=="thickness")thickness=stod(val); else if(key=="Rfar")Rfar=stod(val);
		else if(key=="Re")Re=stod(val); else if(key=="xiStretch")xiStretch=stod(val); else if(key=="dt")dt=stod(val); else if(key=="tol")tol=stod(val);
	}
	filesystem::create_directories(outDir);
	ofstream fp(outDir+"/pressure.csv"), fm(outDir+"/momentum.csv");
	fp<<"Ni,Nj,cells,solver,prec,iters,relres,cpu\n"; fm<<"Ni,Nj,cells,solver,prec,iters,relres,cpu\n";

	vector<int> list; { stringstream ss(sizes); string s; while(getline(ss,s,',')) list.push_back(stoi(s)); }
	vector<pair<string,string>> pc={{"cg","none"},{"cg","jacobi"},{"cg","sgs"},{"cg","ilu0"},{"cg","mg"},
	                                {"mg","-"},{"bicgstab","ilu0"},{"bicgstab","mg"}};
	vector<pair<string,string>> mc={{"cg","none"},{"cg","jacobi"},{"cg","ilu0"},{"bicgstab","jacobi"}};

	auto run=[&](ofstream& f,const Stencil& A,const Vec& b,const string& s,const string& p){
		Vec x(A.n(),0.0);
		auto t0=clk::now();
		SolveInfo si=solveSystem(s,p,A,b,x,tol,20000);
		double cpu=chrono::duration<double>(clk::now()-t0).count();
		f<<A.Ni<<","<<A.Nj<<","<<A.n()<<","<<s<<","<<p<<","<<si.iters<<","<<si.relres<<","<<cpu<<"\n";
		cout<<"  "<<setw(9)<<s<<" + "<<setw(6)<<p<<"  iters "<<setw(6)<<si.iters
		    <<"  cpu "<<fixed<<setprecision(3)<<cpu<<" s"<<endl;
	};

	for(int Ni:list)
	{
		int Nj=Ni/2;
		Mesh mesh=makeConformalMesh(shape,thickness,Ni,Nj,Rfar,xiStretch);
		OGrid G; G.build(mesh);
		mt19937 rng(7); uniform_real_distribution<double> U(-1,1);
		Vec b(G.n()); for(int k=0;k<G.n();k++) b[k]=U(rng)*G.dV[k];

		// (1) pressure: outflow (x>0) half of the outer circle has p = 0
		Stencil P; P.resize(Ni,Nj);
		for(int I=0;I<Ni;I++) for(int J=0;J<Nj;J++)
		{
			int k=G.id(I,J);
			P.aE[k]=G.gE[k]; P.aW[k]=G.gE[G.id(G.im(I),J)];
			if(J>0)P.aS[k]=G.gN[k-1]; if(J<Nj-1)P.aN[k]=G.gN[k];
			P.aP[k]=P.aE[k]+P.aW[k]+P.aS[k]+P.aN[k];
			if(J==Nj-1 && G.xc[k]>0) P.aP[k]+=G.gN[k];
		}
		cout<<"\n=== "<<Ni<<" x "<<Nj<<" ("<<G.n()<<" cells): pressure ===\n";
		for(auto& c:pc) run(fp,P,b,c.first,c.second);

		// (2) momentum: V/dt + 1/2 nu * diffusion
		double nu=1.0/Re;
		Stencil M; M.resize(Ni,Nj);
		for(int I=0;I<Ni;I++) for(int J=0;J<Nj;J++)
		{
			int k=G.id(I,J);
			M.aE[k]=0.5*nu*G.gE[k]; M.aW[k]=0.5*nu*G.gE[G.id(G.im(I),J)];
			if(J>0)M.aS[k]=0.5*nu*G.gN[k-1]; if(J<Nj-1)M.aN[k]=0.5*nu*G.gN[k];
			M.aP[k]=G.dV[k]/dt+M.aE[k]+M.aW[k]+M.aS[k]+M.aN[k];
			if(J==0) M.aP[k]+=0.5*nu*G.gW[I];
		}
		cout<<"=== momentum ===\n";
		for(auto& c:mc) run(fm,M,b,c.first,c.second);
	}
	return 0;
}




