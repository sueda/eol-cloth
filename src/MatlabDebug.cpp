#include "MatlabDebug.h"

#include <fstream>
#include <iomanip>

using namespace std;
using namespace Eigen;

//template <typename Derived>
//
//void sparse_to_file_as_dense(const Eigen::EigenBase<Derived>& mat, string var_name)
//{
//	MatrixXd dMat;
//	dMat = MatrixXd(mat);
//
//	ofstream ofs;
//	ofs.open("test.m", ofstream::out | ofstream::app);
//
//	ofs << var_name;
//	ofs << " = [";
//	ofs << dMat;
//	ofs << "];\n\n";
//
//	ofs.close();
//}

// This is messy
void sparse_to_file_as_sparse_m(const SparseMatrix<double>& mat, string var_name)
{
	MatrixXd dMat;
	dMat = MatrixXd(mat);

	ofstream ofs;
	ofs.open("test.m", ofstream::out | ofstream::app);

	ofs << "i = [ ";

	for (int i = 0; i < dMat.rows(); i++) {
		for (int j = 0; j < dMat.cols(); j++) {
			if (dMat(i, j) != 0.0) {
				ofs << i + 1 << " ";
			}
		}
	}

	ofs << "]';\n";

	ofs << "j = [ ";

	for (int i = 0; i < dMat.rows(); i++) {
		for (int j = 0; j < dMat.cols(); j++) {
			if (dMat(i, j) != 0.0) {
				ofs << j + 1 << " ";
			}
		}
	}

	ofs << "]';\n";

	ofs << "v = [ ";

	for (int i = 0; i < dMat.rows(); i++) {
		for (int j = 0; j < dMat.cols(); j++) {
			if (dMat(i, j) != 0.0) {
				ofs << setprecision(16) << dMat(i,j) << " ";
			}
		}
	}

	ofs << "]';\n";

	ofs << var_name;
	ofs << " = sparse(i,j,v,";
	ofs << dMat.rows() << "," << dMat.cols() << ");\n\n";

	ofs.close();
}

void sparse_to_file_as_dense(const SparseMatrix<double>& mat, string var_name)
{
	MatrixXd dMat;
	dMat = MatrixXd(mat);

	ofstream ofs;
	ofs.open("test.m", ofstream::out | ofstream::app);

	ofs << var_name;
	ofs << " = [";
	ofs << dMat;
	ofs << "];\n\n";

	ofs.close();
}

void mat_to_file(const Eigen::MatrixXd& mat, string var_name)
{
	ofstream ofs;
	ofs.open("test.m", ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = [";
	ofs << mat;
	ofs << "];\n\n";

	ofs.close();
}

void mat_to_file(const Eigen::Matrix4d& mat, string var_name)
{
	ofstream ofs;
	ofs.open("test.m", ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = [";
	ofs << mat;
	ofs << "];\n\n";

	ofs.close();
}

void mat_to_file(const Eigen::MatrixXi& mat, string var_name)
{
	ofstream ofs;
	ofs.open("test.m", ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = [";
	ofs << mat;
	ofs << "];\n\n";

	ofs.close();
}

void vec_to_file(const Eigen::VectorXd& vec, std::string var_name)
{
	ofstream ofs;
	ofs.open("test.m", ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = [\n";
	ofs << vec;
	ofs << "\n];\n\n";

	ofs.close();
}

void vec_to_file(const Eigen::VectorXi& vec, std::string var_name)
{
	ofstream ofs;
	ofs.open("test.m", ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = [";
	ofs << vec;
	ofs << "];\n\n";

	ofs.close();
}

void vec_to_file(const Eigen::Vector3d& vec, std::string var_name)
{
	ofstream ofs;
	ofs.open("test.m", ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = [";
	ofs << vec;
	ofs << "];\n\n";

	ofs.close();
}

void double_to_file(double d, string var_name) 
{
	ofstream ofs;
	ofs.open("test.m", ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = ";
	ofs << d;
	ofs << ";\n\n";

	ofs.close();
}

void EOL_outputter(int EOLE, int EOLV, int verts, int faces)
{
	ofstream ofs;
	ofs.open("timings.csv", ofstream::out | ofstream::app);

	ofs << EOLV << ", " << EOLE << ", " << verts << ", " << faces << ", ";

	ofs.close();
}