#pragma once
#ifndef ____MosekQuadProgram__
#define ____MosekQuadProgram__

#include "QuadProg.h"

#include <memory>
#include <map>
#include <mosek.h>

struct MosekVarBounds;
struct MosekObjectiveMatrix;
struct MosekObjectiveVector;
struct MosekConstraintMatrix;
struct MosekConstraintVector;

class QuadProgMosek : public QuadProg {
private:

	int numVars, numCons, numIneqs, numEqs;
	std::shared_ptr<MosekVarBounds> bounds;

	std::shared_ptr<MosekObjectiveMatrix> objectiveMat;
	std::shared_ptr<MosekObjectiveVector> objectiveVec;
	double objectiveConstant;

	std::shared_ptr<MosekConstraintMatrix> inequalityMat;
	std::shared_ptr<MosekConstraintVector> inequalityVec;

	std::shared_ptr<MosekConstraintMatrix> equalityMat;
	std::shared_ptr<MosekConstraintVector> equalityVec;

	std::shared_ptr<Eigen::VectorXd> lowerVariableBound;
	std::shared_ptr<Eigen::VectorXd> upperVariableBound;

	std::map<MSKiparame, MSKint32t> paramsInt;
	std::map<MSKdparame, MSKrealt> paramsDouble;

public:

	QuadProgMosek();
	virtual ~QuadProgMosek();

	virtual void setNumberOfVariables(int numVars);
	virtual void setNumberOfInequalities(int numIneq);
	virtual void setNumberOfEqualities(int numEq);

	virtual void setObjectiveMatrix(const Eigen::SparseMatrix<double> & mat);
	virtual void setObjectiveVector(const Eigen::VectorXd & vector);
	virtual void setObjectiveConstant(double constant);

	virtual void setLowerVariableBound(const Eigen::VectorXd & bounds);
	virtual void setUpperVariableBound(const Eigen::VectorXd & bounds);

	virtual void setInequalityMatrix(const Eigen::SparseMatrix<double> & mat);
	virtual void setInequalityVector(const Eigen::VectorXd & vector);

	virtual void setEqualityMatrix(const Eigen::SparseMatrix<double> & mat);
	virtual void setEqualityVector(const Eigen::VectorXd & vector);

	virtual bool solve();

	virtual Eigen::VectorXd getPrimalSolution();
	virtual Eigen::VectorXd getDualInequality();
	virtual Eigen::VectorXd getDualEquality();
	virtual Eigen::VectorXd getDualLower();
	virtual Eigen::VectorXd getDualUpper();

	bool setParamInt(int name, int value);
	bool setParamDouble(int name, double value);
};

#endif // ____MosekQuadProgram__