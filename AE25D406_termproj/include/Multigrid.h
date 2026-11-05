#ifndef MULTIGRID_H
#define MULTIGRID_H
#include <vector>

class MultigridLap
{
public:
	enum Smoother {GAUSS_SEIDEL=0, JACOBI=1};
	MultigridLap(int N);//construct
	std::vector<double>solve //solve prob
	(
		int numLevels,
		Smoother smoother,
		int nu1,//pre smoothing
		int nu2,//post smoothing
		double tolerance,
		int maxcycle);
	
	//no multigrid
	std::vector<double>solveSingleGrid
	(
		Smoother smoother,
		double tolerance,
		int maxIter);


	//post process
	static double analytical(double x);//exact solution
	double computeL2Error() const;//error
	double computeLinfError() const;//error
	double getU(int i) const;//get solution
	int getN() const;//grid info
	double getX(int i) const;

private:
	int N_;//divi in x
	std::vector<double>u_;
	int idx(int i) const;
	void applyBoundaryCondi(std::vector<double>& u);//bc
	void smooth(std::vector<double>& u,const std::vector<double>& f,int n,int sweeps,Smoother type);//gauss smoothing /jacobi
	void computeResi(const std::vector<double>& u,const std::vector<double>&f,std ::vector<double>&r,int n);//calculate residue
	void restrictResi(const std::vector<double>& fine,std::vector<double>& coarse, int nFine);//fine grid --- coarse
	void prolongCorrec(const std::vector<double>& coarse,std::vector<double>& fine, int nFine);
	double norm(const std::vector<double>& r) const;
	void vcycle(std::vector<double>& u ,const std::vector<double>&f ,int n, int level ,int numLevels,int nu1, int nu2, Smoother smoother);
};
#endif
