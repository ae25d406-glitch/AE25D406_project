
#ifndef CAVITYFLOW_H
#define CAVITYFLOW_H
#include <string>
#include <vector>
#include "Mesh.h"
#include "PressureCorrection.h"

//staggered grid:pressure p(I,J) at cell center and velocity u(i,J) v(I,j) at face 
//pressure velocity coupling: SIMPLE,PISO,PIMPLE

struct FlowSettings
{
	double Re=100.0;//reynold number
	double Utop=1.0;// top velocity
	double Ubottom=0.0;// bottom velocity
	std::string scheme="HYBRID"; // convection scheme:UPWIND,HYBRID,CENTRAL
	std::string algorithm="SIMPLE";// SIMPLE,PISO,PIMPLE
	double alphaU=0.7;// under relaxation of u and v (SIMPLE/PIMPLE)
	double alphaP=0.3;// under relax for p (SIMPLE/PIMPLE)
	double dt=0.01; // time step (PISO/PIMPLE)
	int nOuter=3;
	int nInner=1;
	std::string pSolver="MG";
	double tol=1e-6;
	int maxIter=100000;
	double maxTime=1e9;
	int printEvery=500;
};

class CavityFlow
{
public:
	CavityFlow(const Mesh& mesh,const FlowSettings& set);
	void solve();
	void writeResults(const Mesh& mesh,const std::string& prefix)const;

	bool converged=false;
	int iterations=0;
	double time=0.0;
	double cpuSeconds=0.0;
	//for regresion test
	const Field& getU()const{return u;}//u(i,J)
    const Field& getV()const{return v;}//v(I,j)
    const Field& getP()const{return p;}//p(I,J)

private:
// momentum coeff 

	struct MomentumEquation
	{
		Field aW,aE,aS,aN,aP,b;
	};

	FlowSettings set;
	int Nx,Ny;
	double dx,dy,nu;
	Field u,v,p;
	Field uOld,vOld;
	Field du,dv;
	MomentumEquation U,V;
	std::vector<std::vector<double>>history;// residual history 

	double neighbourCoeff(double F,double D)const;
	void buildMomentum(double alpha,double dtStep,MomentumEquation& Ue,MomentumEquation& Ve)const;
	void computeD();// d=A/a_p
	void solveMomentum(int sweeps);// u*,v*

//pressure correction

	PressureEquation buildPressureCorrection(const Field& uStar,const Field& vStar)const;
	Field solvePressure(const PressureEquation& eq)const;
	void correctVelocity(const Field& pc);
	Field neighbourSumU(const Field& du_)const;
	Field neighbourSumV(const Field& dv_)const;

	void iterationSIMPLE();
	void stepPISO();
	void stepPIMPLE();
	void secondCorrection(const Field& uStar,const Field& vStar,double alphaP);

//convergence 

	void steadyResiduals(double& Ru,double& Rv,double& Rm);
};

#endif

