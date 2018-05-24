#include "QuadProgMosek.h"

#ifdef _MEX_
#include "mex.h"
#endif

struct MosekVarBoundData {
	MSKint32t variableIndex;
	MSKboundkeye key;
	double lowerBound, upperBound;

	MosekVarBoundData(MSKint32t variableIndex, MSKboundkeye key, double lowerBound, double upperBound)
		: variableIndex(variableIndex), key(key), lowerBound(lowerBound), upperBound(upperBound) {}
};

struct MosekVarBounds {
	std::vector<MSKboundkeye> keys; // bound key
	std::vector<double>lowerBounds; // lower bound
	std::vector<double>upperBounds; // upper bound

	MosekVarBounds() {}

	void addBounds(MSKboundkeye key, double lowerBound, double upperBound) {
		keys.push_back(key);
		lowerBounds.push_back(lowerBound);
		upperBounds.push_back(upperBound);
	}

	MosekVarBoundData getBoundData(MSKint32t variableIndex) {
		return MosekVarBoundData(variableIndex, keys[variableIndex],
			lowerBounds[variableIndex], upperBounds[variableIndex]);
	}
};

struct MosekObjectiveMatrixData {
	MSKint32t numValues, *qsubi, *qsubj;
	double * qval;

	MosekObjectiveMatrixData(MSKint32t numValues, MSKint32t * qsubi, MSKint32t * qsubj, double * qval)
		: numValues(numValues), qsubi(qsubi), qsubj(qsubj), qval(qval) {}
};

struct MosekObjectiveMatrix {
	std::vector<MSKint32t>  rowIndices;
	std::vector<MSKint32t>  colIndices;
	std::vector<double>     values;

	MosekObjectiveMatrix() {}

	static std::shared_ptr<MosekObjectiveMatrix>
		fromSparseMatrix(const Eigen::SparseMatrix<double> & mat) {
		auto Q = std::shared_ptr<MosekObjectiveMatrix>(new MosekObjectiveMatrix());
		for (int row = 0; row < mat.outerSize(); ++row) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(mat, row); it; ++it) {
				if (it.row() >= it.col()) {
					Q->addTriplet(it.row(), it.col(), it.value());
				}
			}
		}
		return Q;
	}

	static std::shared_ptr<MosekObjectiveMatrix>
		fromDenseMatrix(const Eigen::MatrixXd & mat) {
		auto Q = std::shared_ptr<MosekObjectiveMatrix>(new MosekObjectiveMatrix());
		for (int row = 0; row < mat.rows(); ++row) {
			/// We only take the lower triangle
			for (int col = 0; col < mat.cols() && row >= col; col++) {
				if (mat(row, col) != 0) {
					Q->addTriplet(row, col, mat(row, col));
				}
			}
		}
		return Q;
	}

	void addTriplet(MSKint32t row, MSKint32t col, double val) {
		rowIndices.push_back(row);
		colIndices.push_back(col);
		values.push_back(val);
	}

	MosekObjectiveMatrixData getData() {
		MSKint32t * qsubi_p = NULL, *qsubj_p = NULL;
		double * qval_p = NULL;
		if (values.size() != 0) {
			qsubi_p = &rowIndices[0];
			qsubj_p = &colIndices[0];
			qval_p = &values[0];
		}
		return MosekObjectiveMatrixData(values.size(), qsubi_p, qsubj_p, qval_p);
	}
};

typedef MosekObjectiveVector MosekObjectiveVectorData;

struct MosekObjectiveVector {
	Eigen::VectorXd c;

	MosekObjectiveVector() {}
	MosekObjectiveVectorData getData() { return *this; }

	static std::shared_ptr<MosekObjectiveVector> fromVector(const Eigen::VectorXd & vector) {
		auto c = std::shared_ptr<MosekObjectiveVector>(new MosekObjectiveVector());
		c->c = vector;
		return c;
	}
};

struct MosekConstraintVectorData {
	MSKint32t columnNum, numNonZerosInCol, *rowIndexPtr;
	double * valuePtr;

	MosekConstraintVectorData(MSKint32t columnNum, MSKint32t numNonZerosInCol, MSKint32t * rowIndexPtr,
		double * valuePtr) : columnNum(columnNum), numNonZerosInCol(numNonZerosInCol), rowIndexPtr(rowIndexPtr),
		valuePtr(valuePtr) {}
};

struct MosekConstraintVector {
	std::vector<MSKint32t> rowIndices;
	std::vector<double> values;

	MosekConstraintVector() {}

	void addValueAtRow(MSKint32t rowIndex, double value) {
		rowIndices.push_back(rowIndex);
		values.push_back(value);
	}

	bool empty() {
		return values.size() == 0;
	}

	MosekConstraintVectorData getData(MSKint32t col) {
		MSKint32t numNonZerosInCol = static_cast<MSKint32t>(values.size());
		MSKint32t * rowIndexPtr = (numNonZerosInCol > 0) ? &rowIndices[0] : NULL;
		double * valuePtr = (numNonZerosInCol > 0) ? &values[0] : NULL;
		return MosekConstraintVectorData(col, numNonZerosInCol, rowIndexPtr, valuePtr);
	}

	static std::shared_ptr<MosekConstraintVector> fromVector(const Eigen::VectorXd & vector) {
		auto Aj = std::shared_ptr<MosekConstraintVector>(new MosekConstraintVector());
		for (int i = 0; i < vector.rows(); i++) {
			Aj->addValueAtRow(i, vector(i));
		}
		return Aj;
	}
};

typedef MosekConstraintMatrix MosekConstraintMatixData;

struct MosekConstraintMatrix {
	std::vector< std::pair<MSKint32t, std::shared_ptr<MosekConstraintVector> > > columns;

	MosekConstraintMatrix() {}

	MosekConstraintMatixData getData() {
		return *this;
	}

	static std::shared_ptr<MosekConstraintMatrix>
		fromSparseMatrix(const Eigen::SparseMatrix<double> & mat) {
		auto A = std::shared_ptr<MosekConstraintMatrix>(new MosekConstraintMatrix());
		auto A_map = std::map< int, std::shared_ptr<MosekConstraintVector> >();

		for (int row = 0; row < mat.outerSize(); row++) {
			auto Ak = std::shared_ptr<MosekConstraintVector>(new MosekConstraintVector());
			for (Eigen::SparseMatrix<double>::InnerIterator it(mat, row); it; ++it) {
				if (A_map.count(it.col()) == 0) {
					A_map[row] = std::shared_ptr<MosekConstraintVector>(new MosekConstraintVector());
				}
				A_map[row]->addValueAtRow(static_cast<MSKint32t>(it.row()), it.value());
			}
		}

		for (auto itr = A_map.begin(); itr != A_map.end(); itr++) {
			A->columns.push_back(std::make_pair(static_cast<MSKint32t>(itr->first), itr->second));
		}

		return A;
	}

	static std::shared_ptr<MosekConstraintMatrix>
		fromDenseMatrix(const Eigen::MatrixXd & mat) {
		auto A = std::shared_ptr<MosekConstraintMatrix>(new MosekConstraintMatrix());
		for (int col = 0; col < mat.cols(); col++) {
			auto Ak = std::shared_ptr<MosekConstraintVector>(new MosekConstraintVector());
			for (int row = 0; row < mat.rows(); row++) {
				if (mat(row, col) != 0) {
					Ak->addValueAtRow(static_cast<MSKint32t>(row), mat(row, col));
				}
			}
			if (!Ak->empty()) {
				A->columns.push_back(std::make_pair(static_cast<MSKint32t>(col), Ak));
			}
		}
		return A;
	}
};

///
///
///
///
static MSKenv_t __mosek_env = NULL;
static MSKenv_t __getMosekEnv() {
	return __mosek_env;
}

/// Sets up the environment, if needed, and returns if it was successful.
static MSKrescodee __setupMosekEnvIfNeeded() {
	MSKrescodee result = MSK_RES_OK;
	if (__mosek_env == NULL) {
		result = MSK_makeenv(&__mosek_env, NULL);
		if (result == MSK_RES_OK) {
			result = MSK_putlicensepath(__mosek_env, kMosekLicensePath);
		}
		assert(result == MSK_RES_OK);
	}
	return result;
}

static MSKtask_t __mosek_task = NULL;
static MSKtask_t __getMosekTask() {
	return __mosek_task;
}
static void __setMosekTask(const MSKtask_t &task) {
	__mosek_task = task;
}

/// Used as logging output for our quadratic program
static void MSKAPI __mosekLog(void *handle, MSKCONST char str[]) {
#ifdef _MEX_
	mexPrintf("%s\n", str); // NOT WORKING!
#else
	printf("%s\n", str);
	fflush(stdout);
#endif
}

///
///
///
///

QuadProgMosek::QuadProgMosek() : QuadProg() {
	objectiveConstant = 0.0;
	numVars = 0;
	numCons = 0;
	numIneqs = 0;
	numEqs = 0;
}

QuadProgMosek::~QuadProgMosek() {}

///
///
///
void QuadProgMosek::setNumberOfVariables(int numVars) {
	this->numVars = numVars;
}

void QuadProgMosek::setLowerVariableBound(const Eigen::VectorXd & bounds) {
	lowerVariableBound = std::shared_ptr<Eigen::VectorXd>(new Eigen::VectorXd(bounds.rows()));
	*lowerVariableBound = bounds;
}

void QuadProgMosek::setUpperVariableBound(const Eigen::VectorXd & bounds) {
	upperVariableBound = std::shared_ptr<Eigen::VectorXd>(new Eigen::VectorXd(bounds.rows()));
	*upperVariableBound = bounds;
}

///
///
///

void QuadProgMosek::setObjectiveMatrix(const Eigen::SparseMatrix<double> & mat) {
	objectiveMat = MosekObjectiveMatrix::fromSparseMatrix(mat);
}

void QuadProgMosek::setObjectiveVector(const Eigen::VectorXd & vector) {
	objectiveVec = MosekObjectiveVector::fromVector(vector);
}

void QuadProgMosek::setObjectiveConstant(double constant) {
	objectiveConstant = constant;
}

///
///
///
///

void QuadProgMosek::setNumberOfInequalities(int numIneqs) {
	this->numIneqs = numIneqs;
	numCons = numIneqs + numEqs;
}

void QuadProgMosek::setInequalityMatrix(const Eigen::SparseMatrix<double> & mat) {
	assert(mat.rows() == numIneqs);
	inequalityMat = MosekConstraintMatrix::fromSparseMatrix(mat);
}

void QuadProgMosek::setInequalityVector(const Eigen::VectorXd & vector) {
	assert(vector.rows() == numIneqs);
	inequalityVec = MosekConstraintVector::fromVector(vector);
}

///
///
///
///

void QuadProgMosek::setNumberOfEqualities(int numEqs) {
	this->numEqs = numEqs;
	numCons = numIneqs + numEqs;
}

void QuadProgMosek::setEqualityMatrix(const Eigen::SparseMatrix<double> & mat) {
	assert(mat.rows() == numEqs);
	equalityMat = MosekConstraintMatrix::fromSparseMatrix(mat);
}

void QuadProgMosek::setEqualityVector(const Eigen::VectorXd & vector) {
	assert(vector.rows() == numEqs);
	equalityVec = MosekConstraintVector::fromVector(vector);
}

///
///
///
///

struct DoNextTask {
	bool doIt;
	std::function<MSKrescodee(void)> nextOp;
	static void DoNoThing() {}

	DoNextTask(bool doIt = true) : doIt(doIt) {}

	DoNextTask doNext(std::function<MSKrescodee(void)> task) {
		bool doNext = false;
		if (doIt) {
			MSKrescodee result = task();
			if (result != MSK_RES_OK) {
				char symname[MSK_MAX_STR_LEN];
				char desc[MSK_MAX_STR_LEN];
				printf("An error occurred while optimizing.\n");
				MSK_getcodedesc(result, symname, desc);
				printf("Error %s - \"%s\"\n", symname, desc);
			}
			doNext = result == MSK_RES_OK;
		}
		return DoNextTask(doNext);
	}
};

bool QuadProgMosek::solve() {

	const MSKint32t kNumVars = static_cast<MSKint32t>(numVars);
	const MSKint32t kNumCons = static_cast<MSKint32t>(numCons);

	MSKrescodee   r = __setupMosekEnvIfNeeded();
	MSKenv_t      env = __getMosekEnv();
	MSKtask_t     task = __getMosekTask();

	DoNextTask(true).doNext([&]() {
		MSKrescodee result = MSK_RES_OK;
		if (task != NULL) {
			result = MSK_deletetask(&task);
		}
		return result;
	}).doNext([&]() {
		MSKrescodee result = MSK_maketask(env, kNumCons, kNumVars, &task);
		__setMosekTask(task);
		return result;
	}).doNext([&]() {
		// Log to stdout, file, or none
		//return MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, __mosekLog);
		return MSK_linkfiletotaskstream(task, MSK_STREAM_LOG, "mosek.log", 0);
		//return true;
	}).doNext([&]() {
		// Set parameters
		MSKrescodee result = MSK_RES_OK;
		for (auto it = paramsInt.begin(); it != paramsInt.end() && result == MSK_RES_OK; ++it) {
			result = MSK_putintparam(task, it->first, it->second);
		}
		for (auto it = paramsDouble.begin(); it != paramsDouble.end() && result == MSK_RES_OK; ++it) {
			result = MSK_putdouparam(task, it->first, it->second);
		}
		return result;
	}).doNext([&]() {

		/* Append 'NUMCON' empty constraints.
		The constraints will initially have no bounds. */
		return MSK_appendcons(task, kNumCons);
	}).doNext([&]() {

		/* Append 'NUMVAR' variables.
		The variables will initially be fixed at zero (x=0). */
		return MSK_appendvars(task, kNumVars);
	}).doNext([&]() {

		/* Optionally add a constant term to the objective. */
		return MSK_putcfix(task, this->objectiveConstant);
	}).doNext([&]() {
		MSKrescodee result = MSK_RES_OK;

		if (objectiveVec) {
			auto cData = objectiveVec->getData();
			for (MSKint32t j = 0; j < cData.c.size() && result == MSK_RES_OK; ++j) {
				/* Set the linear term c_j in the objective.*/
				result = MSK_putcj(task, j, cData.c[j]);
			}
		}
		return result;
	}).doNext([&]() {
		// First, set all variables to be FREE. Lower and upper bounds will be added next.
		MSKrescodee result = MSK_RES_OK;
		for (MSKint32t j = 0; j < kNumVars && result == MSK_RES_OK; ++j) {
			result = MSK_putvarbound(task, j, MSK_BK_FR, -MSK_INFINITY, MSK_INFINITY);
		}
		return result;
	}).doNext([&]() {
		// Set lower/upper bound
		MSKrescodee result = MSK_RES_OK;
		const double infinity = std::numeric_limits<double>::infinity();
		//                           [ (l, u),   (-inf, u),  (l, inf), (-inf, inf) ]
		const MSKboundkeye keys[4] = { MSK_BK_FX, MSK_BK_UP, MSK_BK_LO, MSK_BK_FR };

		if ((lowerVariableBound == 0) != (upperVariableBound == 0)) {
			printf("Both lower and upper bounds must be set.\n");
			return MSK_RES_ERR_UNKNOWN;
		}

		if (lowerVariableBound && upperVariableBound) {
			assert(lowerVariableBound->size() == upperVariableBound->size());
			for (MSKint32t j = 0; j < kNumVars && result == MSK_RES_OK; ++j) {
				MSKboundkeye key;
				double lb = (*lowerVariableBound)(j);
				double ub = (*upperVariableBound)(j);
				const double infinity = std::numeric_limits<double>::infinity();
				// [ (l, u), (-inf, u), (l, inf), (-inf, inf) ]
				const MSKboundkeye keys[4] = { MSK_BK_FX,MSK_BK_UP,MSK_BK_LO,MSK_BK_FR };
				bool lb_is_inf = false, ub_is_inf = false;

				if (lb == -infinity) {
					lb = -MSK_INFINITY;
					lb_is_inf = true;
				}

				if (ub == infinity) {
					ub = +MSK_INFINITY;
					ub_is_inf = true;
				}

				key = keys[2 * (ub_is_inf ? 1 : 0) + (lb_is_inf ? 1 : 0)];
				result = MSK_putvarbound(task, j, key, lb, ub);
			}
		}

		return result;
	}).doNext([&]() {
		// Inequality constraint
		MSKrescodee result = MSK_RES_OK;

		if (inequalityMat) {
			auto aData = inequalityMat->getData();
			for (auto itr = aData.columns.begin(); itr != aData.columns.end() && result == MSK_RES_OK; itr++) {
				auto col = itr->second->getData(itr->first);
				result = MSK_putacol(
					task,
					col.columnNum,
					col.numNonZerosInCol,
					col.rowIndexPtr,
					col.valuePtr);
			}
		}

		if (inequalityVec) {
			auto bData = inequalityVec->getData(0);
			for (MSKint32t i = 0; i < bData.numNonZerosInCol && result == MSK_RES_OK; ++i) {
				result = MSK_putconbound(
					task,
					bData.rowIndexPtr[i],  /* Index of constraint.*/
					MSK_BK_UP,             /* Bound key.*/
					-MSK_INFINITY,         /* Numerical value of lower bound.*/
					bData.valuePtr[i]);    /* Numerical value of upper bound.*/
			}
		}

		return result;
	}).doNext([&]() {
		// Equality constraint
		MSKrescodee result = MSK_RES_OK;

		if (equalityMat) {
			auto aData = equalityMat->getData();
			for (auto itr = aData.columns.begin(); itr != aData.columns.end() && result == MSK_RES_OK; itr++) {
				auto col = itr->second->getData(itr->first);
				// Equalities come after inequalities, so add numIneq to the row indices
				for (int k = 0; k < col.numNonZerosInCol && result == MSK_RES_OK; ++k) {
					result = MSK_putaij(
						task,
						col.rowIndexPtr[k] + numIneqs,
						col.columnNum,
						col.valuePtr[k]);
				}
			}
		}

		if (equalityVec) {
			auto bData = equalityVec->getData(0);
			for (MSKint32t i = 0; i < bData.numNonZerosInCol && result == MSK_RES_OK; ++i) {
				// Equalities come after inequalities, so add numIneq to the row indices
				result = MSK_putconbound(
					task,
					bData.rowIndexPtr[i] + numIneqs,  /* Index of constraint.*/
					MSK_BK_FX,                        /* Bound key.*/
					bData.valuePtr[i],                /* Numerical value of lower bound.*/
					bData.valuePtr[i]);               /* Numerical value of upper bound.*/
			}
		}

		return result;
	}).doNext([&]() {
		/* Input the Q for the objective. */

		auto Qdata = this->objectiveMat->getData();
		return MSK_putqobj(task, Qdata.numValues, Qdata.qsubi, Qdata.qsubj, Qdata.qval);
	}).doNext([&]() {
		MSKrescodee trmcode;
		/* Run optimizer */

		return MSK_optimizetrm(task, &trmcode);
	});

	MSKsolstae solsta;
	MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
	return solsta == MSK_SOL_STA_OPTIMAL || solsta == MSK_SOL_STA_NEAR_OPTIMAL;
}

Eigen::VectorXd QuadProgMosek::getPrimalSolution() {

	Eigen::VectorXd x(numVars);

	MSKtask_t task = __getMosekTask();
	if (task == NULL) {
		return x;
	}

	MSK_getxx(task,
		MSK_SOL_ITR,    /* Request the interior solution. */
		&x[0]);

	return x;
}

Eigen::VectorXd QuadProgMosek::getDualInequality() {

	Eigen::VectorXd y(numIneqs);

	MSKtask_t task = __getMosekTask();
	if (task == NULL) {
		return y;
	}

	DoNextTask(true).doNext([&]() {
		MSKrescodee result = MSK_RES_OK;
		for (int k = 0; k < numIneqs && result == MSK_RES_OK; ++k) {
			result = MSK_getyslice(
				task,
				MSK_SOL_ITR,
				k,
				k + 1,
				&y[k]);
		}
		return result;
	});

	y = -y; // To match matlab's output

	return y;
}

Eigen::VectorXd QuadProgMosek::getDualEquality() {

	Eigen::VectorXd y(numEqs);

	MSKtask_t task = __getMosekTask();
	if (task == NULL) {
		return y;
	}

	DoNextTask(true).doNext([&]() {
		MSKrescodee result = MSK_RES_OK;
		for (int k = 0; k < numEqs && result == MSK_RES_OK; ++k) {
			result = MSK_getyslice(
				task,
				MSK_SOL_ITR,
				k + numIneqs,
				k + numIneqs + 1,
				&y[k]);
		}
		return result;
	});

	y = -y; // To match matlab's output

	return y;
}

Eigen::VectorXd QuadProgMosek::getDualLower() {

	Eigen::VectorXd y(numVars);

	MSKtask_t task = __getMosekTask();
	if (task == NULL) {
		return y;
	}

	DoNextTask(true).doNext([&]() {
		MSKrescodee result = MSK_RES_OK;
		for (int k = 0; k < numVars && result == MSK_RES_OK; ++k) {
			result = MSK_getslxslice(
				task,
				MSK_SOL_ITR,
				k,
				k + 1,
				&y[k]);
		}
		return result;
	});

	return y;
}

Eigen::VectorXd QuadProgMosek::getDualUpper() {

	Eigen::VectorXd y(numVars);

	MSKtask_t task = __getMosekTask();
	if (task == NULL) {
		return y;
	}

	DoNextTask(true).doNext([&]() {
		MSKrescodee result = MSK_RES_OK;
		for (int k = 0; k < numVars && result == MSK_RES_OK; ++k) {
			result = MSK_getsuxslice(
				task,
				MSK_SOL_ITR,
				k,
				k + 1,
				&y[k]);
		}
		return result;
	});

	return y;
}

bool QuadProgMosek::setParamInt(int name, int value) {
	paramsInt.insert(std::pair<MSKiparame, MSKint32t>((MSKiparame)name, (MSKint32t)value));
	return true;
}

bool QuadProgMosek::setParamDouble(int name, double value) {
	paramsDouble.insert(std::pair<MSKdparame, MSKrealt>((MSKdparame)name, (MSKrealt)value));
	return true;
}
