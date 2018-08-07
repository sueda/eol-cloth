#include "GeneralizedSolver.h"

#include "matlabOutputs.h"

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

void generateMatlab(const SparseMatrix<double>& MDK, const VectorXd& b,
	const SparseMatrix<double>& Aeq, const VectorXd& beq,
	const SparseMatrix<double>& Aineq, const VectorXd& bineq,
	VectorXd& v)
{
	mat_s2s_file(MDK, "MDK", "solver.m", true);
	vec_to_file(b, "b", "solver.m", false);
	mat_s2s_file(Aeq, "Aeq", "solver.m", false);
	vec_to_file(beq, "beq", "solver.m", false);
	mat_s2s_file(Aineq, "Aineq", "solver.m", false);
	vec_to_file(bineq, "bineq", "solver.m", false);
	vec_to_file(v, "v_input", "solver.m", false);
}

#ifdef EOLC_MOSEK
bool mosekSolve(const SparseMatrix<double>& MDK, const VectorXd& b,
	const SparseMatrix<double>& Aeq, const VectorXd& beq,
	const SparseMatrix<double>& Aineq, const VectorXd& bineq,
	VectorXd& v)
{
	QuadProgMosek *program = new QuadProgMosek();
	double inf = numeric_limits<double>::infinity();

	VectorXd xl;
	VectorXd xu;
	xl.setConstant(b.size(), -inf);
	xu.setConstant(b.size(), inf);

	program->setNumberOfVariables(b.size());
	program->setNumberOfEqualities(beq.size());
	program->setNumberOfInequalities(bineq.size());

	program->setObjectiveMatrix(MDK);
	program->setObjectiveVector(b);

	program->setInequalityMatrix(Aineq);
	program->setInequalityVector(bineq);

	program->setEqualityMatrix(Aeq);
	program->setEqualityVector(beq);

	bool success = program->solve();

	v = program->getPrimalSolution();

	return success;
}
#endif

#ifdef EOLC_GUROBI
bool gurobiSolve(SparseMatrix<double>& MDK, const VectorXd& b,
	SparseMatrix<double>& Aeq, const VectorXd& beq,
	SparseMatrix<double>& Aineq, const VectorXd& bineq,
	VectorXd& v)
{
	GurobiSparse qp(b.size(), beq.size(), bineq.size());
	qp.displayOutput(false);

	SparseMatrix<double> Sb(b.sparseView());
	SparseVector<double> Sbeq(beq.sparseView());
	SparseVector<double> Sbineq(bineq.sparseView());

	MDK.makeCompressed();
	Sb.makeCompressed();
	Aeq.makeCompressed();
	Aineq.makeCompressed();

	VectorXd XL, XU;
	double inf = numeric_limits<double>::infinity();
	XL.setConstant(b.size(), -inf);
	XU.setConstant(b.size(), inf);

	bool success = qp.solve(MDK, Sb,
		Aeq, Sbeq,
		Aineq, Sbineq,
		XL, XU);

	v = qp.result();

	return success;
}
#endif

bool GeneralizedSolver::velocitySolve(const bool& fixedPoints, const bool& collisions,
	SparseMatrix<double>& MDK, const VectorXd& b,
	SparseMatrix<double>& Aeq, const VectorXd& beq,
	SparseMatrix<double>& Aineq, const VectorXd& bineq,
	VectorXd& v)
{
	//if (true) {
	//	generateMatlab(MDK, b,
	//		Aeq, beq,
	//		Aineq, bineq,
	//		v);
	//}

	if (!collisions && false) {
		// Simplest case, a cloth with no fixed points and no collisions
		if (!fixedPoints) {
			ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
			cg.compute(MDK);
			v = cg.solve(-b);
			//cout << v << endl;
			return true;
		}
		else {
			// If there are fixed points we can use KKT which is faster than solving a quadprog
			// the assumption is made that MDK and b have been built properly uppon making it here
			LeastSquaresConjugateGradient<SparseMatrix<double> > lscg;
			lscg.compute(MDK);
			v = lscg.solve(-b);
			return true;
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
			bool success = mosekSolve(MDK, b,
				Aeq, beq,
				Aineq, bineq,
				v);
			return success;
#else
			cout << "ERROR:" << endl;
			cout << "Attempting to use the mosek solver without mosek support enabled" << endl;
			abort();
#endif
		}
		else if (whichSolver == GeneralizedSolver::Gurobi) {
#ifdef EOLC_GUROBI
			bool success = gurobiSolve(MDK, b,
				Aeq, beq,
				Aineq, bineq,
				v);
			return success;
#else
			cout << "ERROR:" << endl;
			cout << "Attempting to use the gurobi solver without gurobi support enabled" << endl;
			abort();
#endif
		}
	}

	return false;
}