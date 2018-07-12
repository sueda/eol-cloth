#include "GeneralizedSolver.h"

#include <iostream>

#ifdef EOLC_MOSEK
#include "external\SolverWrappers\Mosek\QuadProgMosek.h"
#endif

#ifdef EOLC_GUROBI
#include "external\SolverWrappers\Gurobi\Gurobi.h"
#endif

using namespace std;
using namespace Eigen;

GeneralizedSolver::GeneralizedSolver() :
	whichSolver(GeneralizedSolver::NoSolver)
{

}

bool GeneralizedSolver::velocitySolve(const bool& fixedPoints, const bool& collisions,
	const SparseMatrix<double>& MDK, const VectorXd& b,
	const SparseMatrix<double>& Aeq, const VectorXd& beq,
	const SparseMatrix<double>& Aineq, const VectorXd& bineq,
	VectorXd& v)
{
	if (!collisions) {
		// Simplest case, a cloth with no fixed points and no collisions
		if (!fixedPoints) {
			ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
			cg.compute(MDK);
			v = cg.solve(-b);
			return true;
		}
		else {
			// If there are fixed points we can use KKT which is faster than solving a quadprog
			// the assumption is made that MDK and b have been built properly uppon making it here
			LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
			lscg.compute(MDK);
			v = lscg.solve(-b);
		}
	}
	else {
		if (whichSolver == GeneralizedSolver::NoSolver) {
			cout << "The simulation has encountered a collision, but a quadratic programming solver has not been specified." << endl;
			cout << "Please either set an external quadratic programming solver, or avoid collisions in your simulation." << endl;
			abort();
		}
		else if (whichSolver == GeneralizedSolver::Mosek) {
#ifdef EOLC_MOSEK

#else
			cout << "ERROR:" << endl;
			cout << "Attempting to use the mosek solver without mosek support enabled" << endl;
			abort();
#endif
		}
		else if (whichSolver == GeneralizedSolver::Gurobi) {
#ifdef EOLC_GUROBI

#else
			cout << "ERROR:" << endl;
			cout << "Attempting to use the gurobi solver without gurobi support enabled" << endl;
			abort();
#endif
		}
	}

	return false;
}