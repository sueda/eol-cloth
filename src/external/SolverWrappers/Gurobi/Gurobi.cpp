// This file is part of EigenQP.
//
// EigenQP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// EigenQP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with EigenQP.  If not, see <http://www.gnu.org/licenses/>.

// associated header
#include "Gurobi.h"

#include <type_traits>

namespace Eigen
{


/**
	*												GurobiCommon
	*/


GurobiCommon::GurobiCommon():
	Q_(),
	C_(),
	Beq_(),
	Bineq_(),
	X_(),
	fail_(0),
	nrvar_(0),
	nreq_(0),
	nrineq_(0),
	iter_(0),
	env_(),
	model_(env_)
{
}


int GurobiCommon::iter() const
{
	return iter_;
}


int GurobiCommon::fail() const
{
	return fail_;
}


const VectorXd& GurobiCommon::result() const
{
	return X_;
}

GurobiCommon::WarmStatus GurobiCommon::warmStart() const
{
	int ws = model_.get(GRB_IntParam_MultiObjMethod);
	return static_cast<WarmStatus>(ws);
}

void GurobiCommon::warmStart(GurobiCommon::WarmStatus warmStatus)
{
	using ut = std::underlying_type<GurobiCommon::WarmStatus>::type;
	model_.set(GRB_IntParam_MultiObjMethod, static_cast<ut>(warmStatus));
}

void GurobiCommon::inform() const
{
	switch(fail_)
	{
		case GRB_LOADED:
			std::cout << "Model is loaded, but no solution information is available."
					  << std::endl;
			break;
		case GRB_OPTIMAL:
			std::cout << "Model was solved to optimality (subject to tolerances), and an optimal solution is available."
					  << std::endl;
			break;
		case GRB_INFEASIBLE:
			std::cout << "Model was proven to be infeasible."
					  << std::endl;
		case GRB_INF_OR_UNBD:
			std::cout << "Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion,"
						" set the DualReductions parameter to 0 and reoptimize."
					  << std::endl;
			break;
		case GRB_UNBOUNDED:
			std::cout << "Model was proven to be unbounded. "
						"Important note: an unbounded status indicates the presence of an unbounded ray that allows the objective to improve without limit. "
						"It says nothing about whether the model has a feasible solution. "
						"If you require information on feasibility, you should set the objective to zero and reoptimize."
					  << std::endl;
			break;
		case GRB_CUTOFF:
			std::cout << "Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. "
						"No solution information is available."
					  << std::endl;
			break;
		case GRB_ITERATION_LIMIT:
			std::cout << "Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter,"
						" or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter."
					  << std::endl;
			break;
		case GRB_NODE_LIMIT:
			std::cout << "Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter."
					  << std::endl;
			break;
		case GRB_TIME_LIMIT:
			std::cout << "Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter."
					  << std::endl;
			break;
		case GRB_SOLUTION_LIMIT:
			std::cout << "Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter."
					  << std::endl;
			break;
		case GRB_INTERRUPTED:
			std::cout << "Optimization was terminated by the user."
					  << std::endl;
			break;
		case GRB_NUMERIC:
			std::cout << "Optimization was terminated due to unrecoverable numerical difficulties."
					  << std::endl;
			break;
		case GRB_SUBOPTIMAL:
			std::cout << "Unable to satisfy optimality tolerances; a sub-optimal solution is available."
					  << std::endl;
			break;
		case GRB_INPROGRESS:
			std::cout << "An asynchronous optimization call was made, but the associated optimization run is not yet complete."
					  << std::endl;
			break;
		case GRB_USER_OBJ_LIMIT:
			std::cout << "User specified an objective limit (a bound on either the best objective or the best bound), and that limit has been reached."
					  << std::endl;
			break;
		default:
			std::cout << "The solver has not been runned yet"
					  << std::endl;
	}
}

void GurobiCommon::displayOutput(bool doDisplay)
{
	model_.set(GRB_IntParam_OutputFlag, doDisplay);
}

double GurobiCommon::feasibilityTolerance() const
{
	return model_.get(GRB_DoubleParam_FeasibilityTol);
}

void GurobiCommon::feasibilityTolerance(double tol)
{
	assert(1e-2 >= tol && tol >= 1e-9);
	model_.set(GRB_DoubleParam_FeasibilityTol, tol);
}

double GurobiCommon::optimalityTolerance() const
{
	return model_.get(GRB_DoubleParam_OptimalityTol);
}

void GurobiCommon::optimalityTolerance(double tol)
{
	assert(1e-2 >= tol && tol >= 1e-9);
	model_.set(GRB_DoubleParam_OptimalityTol, tol);
}

void GurobiCommon::problem(int nrvar, int nreq, int nrineq)
{
	for(int i = 0; i < nrvar_; ++i)
	{
		model_.remove(*(vars_+i));
	}

	for(int i = 0; i < nreq_; ++i)
	{
		model_.remove(*(eqconstr_+i));
	}

	for(int i = 0; i < nrineq_; ++i)
	{
		model_.remove(*(ineqconstr_+i));
	}

	eqvars_.clear();
	ineqvars_.clear();
	lvars_.clear();
	rvars_.clear();

	nrvar_ = nrvar;
	nreq_ = nreq;
	nrineq_ = nrineq;

	Q_.resize(nrvar, nrvar);

	C_.resize(nrvar);
	Beq_.resize(nreq);
	Bineq_.resize(nrineq);
	X_.resize(nrvar);

	vars_ = model_.addVars(nrvar, GRB_CONTINUOUS);

	for(int i=0; i < Q_.size(); ++i)
	{
		lvars_.push_back(*(vars_+i/nrvar));
		rvars_.push_back(*(vars_+i%nrvar));
	}

	eqconstr_ = model_.addConstrs(nreq);
	std::vector<char> eqsense(static_cast<size_t>(nreq), '=');
	model_.set(GRB_CharAttr_Sense, eqconstr_, eqsense.data(), nreq);

	eqvars_.reserve(static_cast<size_t>(nrvar*nreq));
	for(int i = 0; i < nrvar; ++i)
	{
		eqvars_.insert(eqvars_.end(), static_cast<size_t>(nreq), *(vars_+i));
	}

	ineqconstr_ = model_.addConstrs(nrineq);
	std::vector<char> ineqsense(static_cast<size_t>(nrineq), '<');
	model_.set(GRB_CharAttr_Sense, ineqconstr_, ineqsense.data(), nrineq);

	for(int i = 0; i < nrvar; ++i)
	{
		ineqvars_.insert(ineqvars_.end(), static_cast<size_t>(nrineq_), *(vars_+i));
	}
}

/**
	*												GurobiDense
	*/


GurobiDense::GurobiDense()
{ }


GurobiDense::GurobiDense(int nrvar, int nreq, int nrineq)
{
	problem(nrvar, nreq, nrineq);
}


void GurobiDense::problem(int nrvar, int nreq, int nrineq)
{
	GurobiCommon::problem(nrvar, nreq, nrineq);
}

void GurobiDense::updateConstr(GRBConstr* constrs, const std::vector<GRBVar>& vars,
	const Eigen::MatrixXd& A, const Eigen::VectorXd& b, int len)
{
	assert(A.rows() == len);
	assert(b.rows() == len);

	if (len > 0)
	{
		for(int i = 0; i < nrvar_; ++i)
		{
			model_.chgCoeffs(constrs, vars.data()+len*i, A.col(i).data(), static_cast<int>(A.rows()));
		}
	}

	model_.set(GRB_DoubleAttr_RHS, constrs, b.data(), len);
}


bool GurobiDense::solve(const MatrixXd& Q, const VectorXd& C,
	const MatrixXd& Aeq, const VectorXd& Beq,
	const MatrixXd& Aineq, const VectorXd& Bineq,
	const VectorXd& XL, const VectorXd& XU)
{
	//Objective
	GRBQuadExpr qexpr;
	qexpr.addTerms(Q.data(), rvars_.data(), lvars_.data(), static_cast<int>(Q.size()));

	GRBLinExpr lexpr;
	lexpr.addTerms(C.data(), vars_, nrvar_);
	model_.setObjective(0.5*qexpr+lexpr);

	model_.set(GRB_DoubleAttr_LB, vars_, XL.data(), nrvar_);
	model_.set(GRB_DoubleAttr_UB, vars_, XU.data(), nrvar_);

	//Update eq and ineq, column by column
	updateConstr(eqconstr_, eqvars_, Aeq, Beq, nreq_);
	updateConstr(ineqconstr_, ineqvars_, Aineq, Bineq, nrineq_);

	model_.optimize();

	fail_ = model_.get(GRB_IntAttr_Status);
	bool success = fail_ == GRB_OPTIMAL;
	iter_ = model_.get(GRB_IntAttr_BarIterCount);
	double* result = model_.get(GRB_DoubleAttr_X, vars_, nrvar_);
	X_ = Map<VectorXd>(result, nrvar_);

	return success;
}


/**
	*												GurobiSparse
	*/


GurobiSparse::GurobiSparse()
{ }


GurobiSparse::GurobiSparse(int nrvar, int nreq, int nrineq)
{
  problem(nrvar, nreq, nrineq);
}


void GurobiSparse::problem(int nrvar, int nreq, int nrineq)
{
	GurobiCommon::problem(nrvar, nreq, nrineq);
}

void GurobiSparse::updateConstr(GRBConstr* constrs, const std::vector<GRBVar>& vars,
			const Eigen::SparseMatrix<double>& A,
			const Eigen::SparseVector<double>& b, int len)
{
	if(len > 0)
	{
		//Update constrs
		std::vector<double> zeros(static_cast<size_t>(len), 0.0);
		for(int k = 0; k < A.outerSize(); ++k)
		{
			model_.chgCoeffs(constrs, vars.data()+len*k, zeros.data(), len);
			for (SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
			{
				model_.chgCoeff(*(constrs+it.row()), *(vars_+it.col()), it.value());
			}
		}

		//Update RHSes
		for(SparseVector<double>::InnerIterator it(b); it; ++it)
		{
			(constrs+it.row())->set(GRB_DoubleAttr_RHS, it.value());
		}
	}
}


bool GurobiSparse::solve(const SparseMatrix<double>& Q, const SparseVector<double>& C,
	const SparseMatrix<double>& Aeq, const SparseVector<double>& Beq,
	const SparseMatrix<double>& Aineq, const SparseVector<double>& Bineq,
	const VectorXd& XL, const VectorXd& XU)
{

	//Objective: quadratic terms
	GRBQuadExpr qexpr;
	for(int k = 0; k<Q.outerSize(); ++k)
	{
		for (SparseMatrix<double>::InnerIterator it(Q,k); it; ++it)
		{
			qexpr.addTerm(0.5*it.value(), *(vars_+it.row()), *(vars_+it.col()));
		}
	}

	//Objective: linear terms
	for (SparseVector<double>::InnerIterator it(C); it; ++it)
	{
		qexpr.addTerm(it.value(), *(vars_+it.row()));
	}

	model_.setObjective(qexpr);

	//Bounds
	model_.set(GRB_DoubleAttr_LB, vars_, XL.data(), nrvar_);
	model_.set(GRB_DoubleAttr_UB, vars_, XU.data(), nrvar_);

	//Update eq
	updateConstr(eqconstr_, eqvars_, Aeq, Beq, nreq_);
	updateConstr(ineqconstr_, ineqvars_, Aineq, Bineq, nrineq_);

	model_.optimize();

	fail_ = model_.get(GRB_IntAttr_Status);
	bool success = fail_ == GRB_OPTIMAL;
	double* result = model_.get(GRB_DoubleAttr_X, vars_, nrvar_);
	iter_ = model_.get(GRB_IntAttr_BarIterCount);
	X_ = Map<VectorXd>(result, nrvar_);

	return success;
}

} // namespace Eigen
