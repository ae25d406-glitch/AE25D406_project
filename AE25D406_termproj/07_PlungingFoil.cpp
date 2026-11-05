
// ---------------------------------------------------------------------
// 1. BUILD
// ---------------------------------------------------------------------
//   g++ -O2 -std=gnu++17 07_PlungingFoil.cpp src/*.cpp -Iinclude -o plunge
//   g++ -O2 -std=gnu++17 08_SolverStudy.cpp  src/*.cpp -Iinclude -o solverstudy
//   g++ -O2 -std=gnu++17 tests/testOGridPoisson.cpp src/*.cpp -Iinclude -o tests/testOGridPoisson
//
// ---------------------------------------------------------------------
// 2. VERIFICATION (exact solution, order of accuracy)
// ---------------------------------------------------------------------
//   ./tests/testOGridPoisson > results/verify_ogrid.txt
//
// ---------------------------------------------------------------------
// 3. CASE 1: CYLINDER (validation)
// ---------------------------------------------------------------------
//   3.1 Re = 100, vortex shedding
//   ./plunge shape=circle Re=100 Ni=256 Nj=128 Rfar=30 dt=0.005 tEnd=150 perturbTime=3 perturbDeg=10 out=results/cyl_re100
//   3.2 Re = 40, steady
//   ./plunge shape=circle Re=40 Ni=128 Nj=64 Rfar=30 dt=0.02 tEnd=100 out=results/cyl_re40
//
// ---------------------------------------------------------------------
// 4. CASE 2a AND 2b: PLUNGING ELLIPSE, Re = 500, 10 cycles
// ---------------------------------------------------------------------
//   4.1 Case 2a (k = 3, kh0 = 1)
//   ./plunge shape=ellipse thickness=0.12 Re=500 k=3.0 h0=0.3333 Ni=128 Nj=64 xiStretch=0.2 dt=0.001 tEnd=21 out=results/case2a_128
//   4.2 Case 2b (k = 3.33, kh0 = 0.8)
//   ./plunge shape=ellipse thickness=0.12 Re=500 k=3.3333 h0=0.24 Ni=128 Nj=64 xiStretch=0.2 dt=0.001 tEnd=19 out=results/case2b_128
//
// ---------------------------------------------------------------------
// 5. GRID AND TIME-STEP STUDY (Case 2a)
// ---------------------------------------------------------------------
//   5.1 coarse and fine grid (128x64 is run 4.1), dt halved with the cell size
//   ./plunge shape=ellipse thickness=0.12 Re=500 k=3.0 h0=0.3333 Ni=64  Nj=32  xiStretch=0.2 dt=0.002  tEnd=21 out=results/grid_case2a/g064
//   ./plunge shape=ellipse thickness=0.12 Re=500 k=3.0 h0=0.3333 Ni=256 Nj=128 xiStretch=0.2 dt=0.0005 tEnd=21 out=results/grid_case2a/g256
//   5.2 same grid, half the time step
//   ./plunge shape=ellipse thickness=0.12 Re=500 k=3.0 h0=0.3333 Ni=128 Nj=64 xiStretch=0.2 dt=0.0005 tEnd=21 out=results/dt_case2a/dt0.0005
//
// ---------------------------------------------------------------------
// 6. CASE 2c: STROUHAL SWEEP (ellipse, k = 2, 4, 6) AND OTHER CHECKS
// ---------------------------------------------------------------------
//   6.1 all missing runs at once (sweep, grid, dt, start direction), 12 at a time
//   bash cases/run_checks.sh 12
//   6.2 one sweep run by hand: tEnd = 10 periods = 20 pi / k,  dt = 0.001 / max(1, k h0)
//   ./plunge shape=ellipse thickness=0.12 Re=500 k=4 h0=0.25 Ni=128 Nj=64 xiStretch=0.2 dt=0.001 tEnd=15.708 out=results/sweep_ellipse/k4_h0.25
//
// ---------------------------------------------------------------------
// 7. JOUKOWSKI AIRFOIL: WAKE TYPES AND START DIRECTION
// ---------------------------------------------------------------------
//   7.1 wake sweep kh = 0.2 ... 2.4 (always ramp the plunge over one period)
//   bash cases/run_airfoil_wakes.sh
//   7.2 deflected wake, start upward (h0 > 0) and downward (h0 < 0)
//   ./plunge shape=joukowski thickness=0.10 Re=500 k=8 h0=0.30  ramp=0.7854 Ni=128 Nj=64 dt=0.0001614 tEnd=6.2832 blend=0.1 out=results/airfoil_blend/k8_h0.30
//   ./plunge shape=joukowski thickness=0.10 Re=500 k=8 h0=-0.30 ramp=0.7854 Ni=128 Nj=64 dt=0.0001614 tEnd=6.2832 blend=0.1 out=results/airfoil_blend/k8_h-0.30
//
// ---------------------------------------------------------------------
// 8. LINEAR-SOLVER STUDY
// ---------------------------------------------------------------------
//   ./solverstudy shape=ellipse sizes=64,128,256,512 out=results/solver_study


#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include "Mesh.h"
#include "OGrid.h"
#include "PlungingFlow.h"

using namespace std;

int main(int argc,char* argv[])//argc-how many words were typed on the command line.
{
	PlungeSettings set;// all solver settings
	string shape="circle";// body shape
	double thickness=0.12, Rfar=30.0, xiStretch=0.0;// geometry
	int Ni=256, Nj=128;// grid 

	for(int a=1;a<argc;a++)
	{
		string arg=argv[a];
		size_t eq=arg.find('=');// position of '=' in the word
		if(eq==string::npos){cerr<<"Inputs must look like name=value: "<<arg<<endl;return 1;}// no '=' found(error so stop)
		string key=arg.substr(0,eq),val=arg.substr(eq+1); // text before '='-key = "Re",// text after '='- val = "500"

		if(key=="shape")shape=val;//circle / ellipse / joukowski
		else if(key=="thickness")thickness=stod(val);
		else if(key=="Rfar")Rfar=stod(val);
		else if(key=="xiStretch")xiStretch=stod(val);
		else if(key=="Ni")Ni=stoi(val);
		else if(key=="Nj")Nj=stoi(val);
		else if(key=="Re")set.Re=stod(val);
		else if(key=="h0")set.h0=stod(val);
		else if(key=="k")set.k=stod(val);
		else if(key=="aoa")set.aoaDeg=stod(val);
		else if(key=="perturbTime")set.perturbTime=stod(val);
		else if(key=="perturbDeg")set.perturbDeg=stod(val);
		else if(key=="ramp")set.rampTime=stod(val);
		else if(key=="dt")set.dt=stod(val);
		else if(key=="tEnd")set.tEnd=stod(val);
		else if(key=="nCorr")set.nCorr=stoi(val);
		else if(key=="blend")set.blend=stod(val);
		else if(key=="pSolver")set.pSolver=val;//pressure solver and preconditioner
		else if(key=="pPrec")set.pPrec=val;
		else if(key=="mSolver")set.mSolver=val;//mv solver and preconditioner
		else if(key=="mPrec")set.mPrec=val;
		else if(key=="pTol")set.pTol=stod(val);
		else if(key=="print")set.printEvery=stoi(val);
		else if(key=="vtk")set.vtkEvery=stoi(val);
		else if(key=="out")set.outDir=val;
		else{cerr<<"Unknown input: "<<key<<endl;return 1;}
	}

	// conformal O-grid (CoordinateMapping::conformalBody)
	Mesh mesh=makeConformalMesh(shape,thickness,Ni,Nj,Rfar,xiStretch);
	OGrid grid;
	grid.build(mesh);

	cout<<"Body: "<<shape<<"  grid "<<Ni<<" x "<<Nj<<"  Rfar = "<<Rfar
	    <<"  Re = "<<set.Re<<"  h0 = "<<set.h0<<"  k = "<<set.k
	    <<"  St = "<<set.k*set.h0/acos(-1.0)<<"  dt = "<<set.dt<<endl;

		//run solver
	PlungingFlow flow(mesh,grid,set);
	flow.solve();
	flow.writeVTK(set.outDir+"/field_final.vtk");

	ofstream info(set.outDir+"/run_info.txt");
	info<<"shape "<<shape<<"\nthickness "<<thickness<<"\nNi "<<Ni<<"\nNj "<<Nj<<"\nRfar "<<Rfar<<"\nxiStretch "<<xiStretch
	    <<"\nRe "<<set.Re<<"\nh0 "<<set.h0<<"\nk "<<set.k<<"\nramp "<<set.rampTime<<"\ndt "<<set.dt<<"\nblend "<<set.blend<<"\nbodyArea "<<grid.bodyArea<<"\n";
	return 0;
}
