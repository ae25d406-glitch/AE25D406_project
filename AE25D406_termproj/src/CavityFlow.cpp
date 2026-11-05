	#include "CavityFlow.h"
#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;


// SIMPLE
// solve mv eqn with p* to give u* and v*
//solve pressure correction eqn to  get p'
//p = p* + alphaP p',u = u* + d (p'_I-1 - p'_I)

//PISO
// predictor - solve mv  with p*
//correction 1 - p' gives  p** = p* + p',u** = u* + d (p'_I-1 - p'_I)
//correction 2- p'' gives p*** = p** + p'',u*** = u** + sum(a_nb (u**-u*)_nb)/a + d (p''_I-1 - p''_>

// PIMPLE
//outer loop (SIMPLE)
//mv predictor-first pressure correction

static Field zeros(int n1,int n2){return Field(n1,vector<double>(n2,0.0));}

// initialize set up
CavityFlow::CavityFlow(const Mesh& mesh,const FlowSettings& s):set(s)
{
	Nx=mesh.Nx;
	Ny=mesh.Ny;
	dx=mesh.dx;
	dy=mesh.dy;
	nu=1.0/set.Re;
	p=zeros(Nx,Ny);
	u=zeros(Nx+1,Ny);
	v=zeros(Nx,Ny+1);
	du=zeros(Nx+1,Ny);
	dv=zeros(Nx,Ny+1);
	uOld=u;
	vOld=v;
}

// neibhour coeff F- flux and D-diffusion
double CavityFlow::neighbourCoeff(double F,double D)const
{
	if(set.scheme=="UPWIND")return D+max(F,0.0);// upwind scheme
	if(set.scheme=="CENTRAL")return D+0.5*F;// central scheme
	return max(F,max(D+0.5*F,0.0));//hybris scheme
}

// duscritised mv eqn
void CavityFlow::buildMomentum(double alpha,double dtStep,MomentumEquation& Ue,MomentumEquation& Ve)const
{
	Ue={zeros(Nx+1,Ny),zeros(Nx+1,Ny),zeros(Nx+1,Ny),zeros(Nx+1,Ny),zeros(Nx+1,Ny),zeros(Nx+1,Ny)};
	Ve={zeros(Nx,Ny+1),zeros(Nx,Ny+1),zeros(Nx,Ny+1),zeros(Nx,Ny+1),zeros(Nx,Ny+1),zeros(Nx,Ny+1)};

	const double Dx=nu*dy/dx;
	const double Dy=nu*dx/dy;
	const double aP0=(dtStep>0.0)?dx*dy/dtStep:0.0;
//u momentum eqn
	for(int i=1;i<Nx;i++)
		for(int J=0;J<Ny;J++)
		{
			double Fw=0.5*(u[i-1][J]+u[i][J])*dy;
			double Fe=0.5*(u[i][J]+u[i+1][J])*dy;
			double Fs=0.5*(v[i-1][J]+v[i][J])*dx;
			double Fn=0.5*(v[i-1][J+1]+v[i][J+1])*dx;

			double aW=neighbourCoeff(Fw,Dx);
			double aE=neighbourCoeff(-Fe,Dx);
			double aS=0.0,aN=0.0,wall=0.0,b=0.0;

			if(J>0)aS=neighbourCoeff(Fs,Dy);
			else{wall+=2.0*Dy;b+=2.0*Dy*set.Ubottom;}

			if(J<Ny-1)aN=neighbourCoeff(-Fn,Dy);
			else{wall+=2.0*Dy;b+=2.0*Dy*set.Utop;}

			double aP=aW+aE+aS+aN+wall+(Fe-Fw+Fn-Fs)+aP0;
			b+=aP0*uOld[i][J];

			aP/=alpha;
			b+=(1.0-alpha)*aP*u[i][J];

			Ue.aW[i][J]=aW;Ue.aE[i][J]=aE;Ue.aS[i][J]=aS;Ue.aN[i][J]=aN;
			Ue.aP[i][J]=aP;Ue.b[i][J]=b;
		}

// v momentum eqn
	for(int I=0;I<Nx;I++)
		for(int j=1;j<Ny;j++)
		{
			double Fw=0.5*(u[I][j-1]+u[I][j])*dy;
			double Fe=0.5*(u[I+1][j-1]+u[I+1][j])*dy;
			double Fs=0.5*(v[I][j-1]+v[I][j])*dx;
			double Fn=0.5*(v[I][j]+v[I][j+1])*dx;

			double aS=neighbourCoeff(Fs,Dy);
			double aN=neighbourCoeff(-Fn,Dy);
			double aW=0.0,aE=0.0,wall=0.0,b=0.0;

			if(I>0)aW=neighbourCoeff(Fw,Dx);
			else wall+=2.0*Dx;

			if(I<Nx-1)aE=neighbourCoeff(-Fe,Dx);
			else wall+=2.0*Dx;

			double aP=aW+aE+aS+aN+wall+(Fe-Fw+Fn-Fs)+aP0;
			b+=aP0*vOld[I][j];

			aP/=alpha;
			b+=(1.0-alpha)*aP*v[I][j];

			Ve.aW[I][j]=aW;Ve.aE[I][j]=aE;Ve.aS[I][j]=aS;Ve.aN[I][j]=aN;
			Ve.aP[I][j]=aP;Ve.b[I][j]=b;
		}
}

// d=A/a_p zero at wall as the velocity is known 
void CavityFlow::computeD()
{
	du=zeros(Nx+1,Ny);
	dv=zeros(Nx,Ny+1);
	for(int i=1;i<Nx;i++)
		for(int J=0;J<Ny;J++)du[i][J]=dy/U.aP[i][J];
	for(int I=0;I<Nx;I++)
		for(int j=1;j<Ny;j++)dv[I][j]=dx/V.aP[I][j];
}

//gauss seidal
void CavityFlow::solveMomentum(int sweeps)
{
	for(int k=0;k<sweeps;k++)
	{
		for(int i=1;i<Nx;i++)
			for(int J=0;J<Ny;J++)
			{
				double sum=U.b[i][J]+(p[i-1][J]-p[i][J])*dy+U.aW[i][J]*u[i-1][J]+U.aE[i][J]*u[i+1][J];
				if(J>0)sum+=U.aS[i][J]*u[i][J-1];
				if(J<Ny-1)sum+=U.aN[i][J]*u[i][J+1];
				u[i][J]=sum/U.aP[i][J];
			}

		for(int I=0;I<Nx;I++)
			for(int j=1;j<Ny;j++)
			{
				double sum=V.b[I][j]+(p[I][j-1]-p[I][j])*dx+V.aS[I][j]*v[I][j-1]+V.aN[I][j]*v[I][j+1];
				if(I>0)sum+=V.aW[I][j]*v[I-1][j];
				if(I<Nx-1)sum+=V.aE[I][j]*v[I+1][j];
				v[I][j]=sum/V.aP[I][j];
			}
	}
}

PressureEquation CavityFlow::buildPressureCorrection(const Field& uStar,const Field& vStar)const
{
	PressureEquation eq=newPressureEquation(Nx,Ny);

	for(int I=0;I<Nx;I++)
		for(int J=0;J<Ny;J++)
		{
			eq.aW[I][J]=du[I][J]*dy;
			eq.aE[I][J]=du[I+1][J]*dy;
			eq.aS[I][J]=dv[I][J]*dx;
			eq.aN[I][J]=dv[I][J+1]*dx;
			eq.aP[I][J]=eq.aW[I][J]+eq.aE[I][J]+eq.aS[I][J]+eq.aN[I][J];

			eq.b[I][J]=(uStar[I][J]-uStar[I+1][J])*dy+(vStar[I][J]-vStar[I][J+1])*dx;
		}
	return eq;
}

Field CavityFlow::solvePressure(const PressureEquation& eq)const
{
	Field pc=zeros(Nx,Ny);
	if(set.pSolver=="SOR")solveSOR(eq,pc,1.8,1e-2,200000);
	else solveMultigrid(eq,pc,1e-2,200);
	return pc;
}


// u = u + d (p'_I-1,J - p'_I,J),   v = v + d (p'_I,J-1 - p'_I,J)
void CavityFlow::correctVelocity(const Field& pc)
{
	for(int i=1;i<Nx;i++)
		for(int J=0;J<Ny;J++)
			u[i][J]+=du[i][J]*(pc[i-1][J]-pc[i][J]);

	for(int I=0;I<Nx;I++)
		for(int j=1;j<Ny;j++)
			v[I][j]+=dv[I][j]*(pc[I][j-1]-pc[I][j]);
}

//sum(a_nb x_nb) / a_P for interior u face  
Field CavityFlow::neighbourSumU(const Field& x)const
{
	Field s=zeros(Nx+1,Ny);
	for(int i=1;i<Nx;i++)
		for(int J=0;J<Ny;J++)
		{
			double sum=U.aW[i][J]*x[i-1][J]+U.aE[i][J]*x[i+1][J];
			if(J>0)sum+=U.aS[i][J]*x[i][J-1];
			if(J<Ny-1)sum+=U.aN[i][J]*x[i][J+1];
			s[i][J]=sum/U.aP[i][J];
		}
	return s;
}

Field CavityFlow::neighbourSumV(const Field& x)const
{
	Field s=zeros(Nx,Ny+1);
	for(int I=0;I<Nx;I++)
		for(int j=1;j<Ny;j++)
		{
			double sum=V.aS[I][j]*x[I][j-1]+V.aN[I][j]*x[I][j+1];
			if(I>0)sum+=V.aW[I][j]*x[I-1][j];
			if(I<Nx-1)sum+=V.aE[I][j]*x[I+1][j];
			s[I][j]=sum/V.aP[I][j];
		}
	return s;
}

//SIMPLE
void CavityFlow::iterationSIMPLE()
{
	// mv eqn with p* gives u* and v*
	buildMomentum(set.alphaU,0.0,U,V);
	computeD();
	solveMomentum(2);

	// pessure correction eqn p'
	PressureEquation eq=buildPressureCorrection(u,v);
	Field pc=solvePressure(eq);

	// correction pressure
	for(int I=0;I<Nx;I++)
		for(int J=0;J<Ny;J++)
			p[I][J]+=set.alphaP*pc[I][J];
	correctVelocity(pc);
}

//PISO second correction 
void CavityFlow::secondCorrection(const Field& uStar,const Field& vStar,double alphaP)
{
	Field uDiff=zeros(Nx+1,Ny),vDiff=zeros(Nx,Ny+1);
	for(int i=0;i<=Nx;i++)for(int J=0;J<Ny;J++)uDiff[i][J]=u[i][J]-uStar[i][J];
	for(int I=0;I<Nx;I++)for(int j=0;j<=Ny;j++)vDiff[I][j]=v[I][j]-vStar[I][j];

	Field uExtra=neighbourSumU(uDiff);
	Field vExtra=neighbourSumV(vDiff);

	for(int i=0;i<=Nx;i++)for(int J=0;J<Ny;J++)u[i][J]+=uExtra[i][J];
	for(int I=0;I<Nx;I++)for(int j=0;j<=Ny;j++)v[I][j]+=vExtra[I][j];

	PressureEquation eq=buildPressureCorrection(u,v);
	Field pcc=solvePressure(eq);

	for(int I=0;I<Nx;I++)
		for(int J=0;J<Ny;J++)
			p[I][J]+=alphaP*pcc[I][J];
	correctVelocity(pcc);
}

//PISO

void CavityFlow::stepPISO()
{
	uOld=u;
	vOld=v;
// predictor mv with p*
	buildMomentum(1.0,set.dt,U,V);
	computeD();
	solveMomentum(2);
	Field uStar=u,vStar=v;

// correction 1 p' gives p**  and then corrects to get u** and v**
	PressureEquation eq=buildPressureCorrection(u,v);
	Field pc=solvePressure(eq);
	for(int I=0;I<Nx;I++)
		for(int J=0;J<Ny;J++)p[I][J]+=pc[I][J];
	correctVelocity(pc);
// corrector 2 p'' gives p*** corrects and gives u*** and v***
	secondCorrection(uStar,vStar,1.0);
}

//PIMPLE
void CavityFlow::stepPIMPLE()
{
	uOld=u;
	vOld=v;

	for(int outer=0;outer<set.nOuter;outer++)
	{
		// mv predictor 
		buildMomentum(set.alphaU,set.dt,U,V);
		computeD();
		solveMomentum(2);
		Field uStar=u,vStar=v;
        
		// first pressure correction(SIMPLE)
		PressureEquation eq=buildPressureCorrection(u,v);
		Field pc=solvePressure(eq);
		for(int I=0;I<Nx;I++)
			for(int J=0;J<Ny;J++)p[I][J]+=set.alphaP*pc[I][J];
		correctVelocity(pc);

		// inner(PISO) correction
		for(int inner=0;inner<set.nInner;inner++)
			secondCorrection(uStar,vStar,set.alphaP);
	}
}

// residual 
void CavityFlow::steadyResiduals(double& Ru,double& Rv,double& Rm)
{
	MomentumEquation Us,Vs;
	buildMomentum(1.0,0.0,Us,Vs);

	double su=0.0,sv=0.0,sm=0.0;

	for(int i=1;i<Nx;i++)
		for(int J=0;J<Ny;J++)
		{
			double r=Us.b[i][J]+(p[i-1][J]-p[i][J])*dy+Us.aW[i][J]*u[i-1][J]+Us.aE[i][J]*u[i+1][J]-Us.aP[i][J]*u[i][J];
			if(J>0)r+=Us.aS[i][J]*u[i][J-1];
			if(J<Ny-1)r+=Us.aN[i][J]*u[i][J+1];
			su+=r*r;
		}

	for(int I=0;I<Nx;I++)
		for(int j=1;j<Ny;j++)
		{
			double r=Vs.b[I][j]+(p[I][j-1]-p[I][j])*dx+Vs.aS[I][j]*v[I][j-1]+Vs.aN[I][j]*v[I][j+1]-Vs.aP[I][j]*v[I][j];
			if(I>0)r+=Vs.aW[I][j]*v[I-1][j];
			if(I<Nx-1)r+=Vs.aE[I][j]*v[I+1][j];
			sv+=r*r;
		}

	for(int I=0;I<Nx;I++)
		for(int J=0;J<Ny;J++)
		{
			double m=(u[I+1][J]-u[I][J])*dy+(v[I][J+1]-v[I][J])*dx;
			sm+=m*m;
		}

	Ru=sqrt(su/((Nx-1)*Ny))/dy;
	Rv=sqrt(sv/(Nx*(Ny-1)))/dx;
	Rm=sqrt(sm/(Nx*Ny))/dy;
}

//solve main loop 
void CavityFlow::solve()
{
	auto start=chrono::steady_clock::now();

	for(int n=1;n<=set.maxIter;n++)
	{
		if(set.algorithm=="SIMPLE")iterationSIMPLE();
		else if(set.algorithm=="PISO"){stepPISO();time+=set.dt;}
		else if(set.algorithm=="PIMPLE"){stepPIMPLE();time+=set.dt;}
		else{cerr<<"Unknown algorithm "<<set.algorithm<<endl;exit(1);}

		double Ru,Rv,Rm;
		steadyResiduals(Ru,Rv,Rm);
		double cpu=chrono::duration<double>(chrono::steady_clock::now()-start).count();

		double KE=0.0;
		for(int I=0;I<Nx;I++)
			for(int J=0;J<Ny;J++)
			{
				double uc=0.5*(u[I][J]+u[I+1][J]),vc=0.5*(v[I][J]+v[I][J+1]);
				KE+=0.5*(uc*uc+vc*vc)*dx*dy;
			}
		history.push_back({double(n),time,Ru,Rv,Rm,cpu,KE});
		iterations=n;

		if(set.printEvery>0&&(n==1||n%set.printEvery==0))
			cout<<setw(7)<<n<<"  t = "<<setw(8)<<time<<scientific<<setprecision(3)<<"   Ru = "<<Ru<<"   Rv = "<<Rv<<"   Rmass = "<<Rm<<defaultfloat<<endl;

		if(!isfinite(Ru)||Ru>1e6){cout<<"*** diverged at iteration "<<n<<endl;break;}
		if(max(Ru,max(Rv,Rm))<set.tol){converged=true;break;}
		if(time>=set.maxTime)break;
	}

	cpuSeconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();
}

static void writeMatrix(const Field& f,const string& file)
{
	ofstream out(file);
	out<<setprecision(10);
	int n1=f.size(),n2=f[0].size();
	for(int J=0;J<n2;J++)
	{
		for(int I=0;I<n1;I++)out<<f[I][J]<<(I<n1-1?",":"\n");
	}
}

void CavityFlow::writeResults(const Mesh& mesh,const string& prefix)const
{
	ofstream h(prefix+"_history.csv");
	h<<"iter,time,resU,resV,resMass,cpu,KE\n"<<setprecision(10);
	for(auto& row:history)
		h<<row[0]<<","<<row[1]<<","<<row[2]<<","<<row[3]<<","<<row[4]<<","<<row[5]<<","<<row[6]<<"\n";

	ofstream fu(prefix+"_uCentre.csv");
	fu<<"y,u\n0,"<<set.Ubottom<<"\n"<<setprecision(10);
	for(int J=0;J<Ny;J++)
		fu<<(J+0.5)*dy<<","<<0.5*(u[Nx/2][J]+u[(Nx+1)/2][J])<<"\n";
	fu<<1.0<<","<<set.Utop<<"\n";

	ofstream fv(prefix+"_vCentre.csv");
	fv<<"x,v\n0,0\n"<<setprecision(10);
	for(int I=0;I<Nx;I++)
		fv<<(I+0.5)*dx<<","<<0.5*(v[I][Ny/2]+v[I][(Ny+1)/2])<<"\n";
	fv<<"1,0\n";

	Field psi=zeros(Nx+1,Ny+1);
	for(int i=0;i<=Nx;i++)
		for(int J=0;J<Ny;J++)
			psi[i][J+1]=psi[i][J]+u[i][J]*dy;

	Field uc=zeros(Nx,Ny),vc=zeros(Nx,Ny);
	for(int I=0;I<Nx;I++)
		for(int J=0;J<Ny;J++)
		{
			uc[I][J]=0.5*(u[I][J]+u[I+1][J]);
			vc[I][J]=0.5*(v[I][J]+v[I][J+1]);
		}

	writeMatrix(psi,prefix+"_psi.csv");
	writeMatrix(p,prefix+"_p.csv");
	writeMatrix(uc,prefix+"_uc.csv");
	writeMatrix(vc,prefix+"_vc.csv");

	ofstream s(prefix+"_summary.csv");
	s<<"algorithm,scheme,Re,N,Utop,Ubottom,alphaU,alphaP,dt,nOuter,nInner,pSolver,converged,iterations,time,cpu\n";
	s<<set.algorithm<<","<<set.scheme<<","<<set.Re<<","<<Nx<<","<<set.Utop<<","<<set.Ubottom<<","<<set.alphaU<<","<<set.alphaP<<","<<set.dt<<","<<set.nOuter<<","<<set.nInner<<","<<set.pSolver<<","<<converged<<","<<iterations<<","<<time<<","<<cpuSeconds<<"\n";

	ofstream vtk(prefix+".vtk");
	vtk<<"# vtk DataFile Version 3.0\nLid-driven cavity\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	vtk<<"POINTS "<<mesh.getNumberOfNodes()<<" float\n";
	for(auto& node:mesh.nodes)vtk<<node.x<<" "<<node.y<<" 0\n";
	vtk<<"CELLS "<<mesh.getNumberOfCells()<<" "<<5*mesh.getNumberOfCells()<<"\n";
	for(auto& cell:mesh.cells)
		vtk<<"4 "<<cell.nodeIDs[0]<<" "<<cell.nodeIDs[1]<<" "<<cell.nodeIDs[2]<<" "<<cell.nodeIDs[3]<<"\n";
	vtk<<"CELL_TYPES "<<mesh.getNumberOfCells()<<"\n";
	for(int k=0;k<mesh.getNumberOfCells();k++)vtk<<"9\n";

	vtk<<"CELL_DATA "<<mesh.getNumberOfCells()<<"\nSCALARS pressure float 1\nLOOKUP_TABLE default\n";
	for(int J=0;J<Ny;J++)for(int I=0;I<Nx;I++)vtk<<p[I][J]<<"\n";
	vtk<<"VECTORS velocity float\n";
	for(int J=0;J<Ny;J++)for(int I=0;I<Nx;I++)vtk<<uc[I][J]<<" "<<vc[I][J]<<" 0\n";
	vtk<<"POINT_DATA "<<mesh.getNumberOfNodes()<<"\nSCALARS streamfunction float 1\nLOOKUP_TABLE default\n";
	for(int J=0;J<=Ny;J++)for(int i=0;i<=Nx;i++)vtk<<psi[i][J]<<"\n";
}