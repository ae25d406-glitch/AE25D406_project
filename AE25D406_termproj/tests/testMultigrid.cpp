//g++ -std=c++17 -O2 -I include tests/testMultigrid.cpp src/Multigrid.cpp -o testMultigrid./testMultigrid

#include <iostream>
#include "Multigrid.h"

int main()
{
	//number of intervals
	int N=64;
	//create solver
	MultigridLap solver(N);
	
	//solve 4 multigrid 64-32-16-8
	
	std::vector<double>residuals=solver.solve(4,MultigridLap::GAUSS_SEIDEL,2,2,1.0e-8,100);
	//print residual
	std::cout<<"Number of cycles="<<residuals.size()<<std::endl;
	std ::cout <<std::endl;
	for (size_t i=0;i<residuals.size();i++)
	{
		std::cout<< "cycle" << i+1 << "Residual = "<< residuals[i] << std::endl;
	}

	//print error
	
	std::cout<<std::endl;
	std::cout<<"L2 error="<<solver.computeL2Error()<<std::endl;
	std::cout<<"Linf error="<<solver.computeLinfError()<<std::endl;

	//print solution
	std::cout<< std::endl;
	for (int i=0;i<=N;i+=8)
	{
		double x=solver.getX(i);
		double numerical=solver.getU(i);
		double exact=solver.analytical(x);
		std::cout<< "x=" << x << "Numerical = "<< numerical <<" exact = "<< exact << std::endl;
	}
	
	//regression test
	const double tolerance=1.0e-8;
	const double L2_limit=1.0e-3;
	const double Linf_limit=1.0e-3;

	bool converged=!residuals.empty()&& residuals.back()< tolerance;
	bool accurate= solver.computeL2Error()<L2_limit && solver.computeLinfError()<Linf_limit;
	std::cout<< std::endl;
	std::cout<< "converged (residual<" << tolerance<<")?"
			<< (converged ? "yes" : "no") << std::endl;
	std::cout << "Accurate  (errors  < " << L2_limit <<")?"
	          << (accurate ? "yes" : "no") << std::endl;
	if (converged && accurate)
	{
		std::cout<< std::endl << "REGRESSION TEST :PASS" <<std::endl;
		return 0;

	}
	else
	{
		std::cout<< std::endl << "REGRESSION TEST :FAIL" <<std::endl;
		return 1;	
	}
}