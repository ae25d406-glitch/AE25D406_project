
#ifndef PRESSURECORRECTION_H
#define PRESSURECORRECTION_H
#include <vector>


//presuure correction for simple/piso.pimple
// solve - sor and multigrid(MG)

typedef std::vector<std::vector<double>> Field;//field

struct PressureEquation
{
	Field aW,aE,aS,aN,aP,b;//coeff and b
};

PressureEquation newPressureEquation(int nx,int ny);

//solve using sor /MG

int solveSOR(const PressureEquation& eq,Field& x,double omega,double relTol,int maxSweeps);
int solveMultigrid(const PressureEquation& eq,Field& x,double relTol,int maxCycles);
#endif
