#pragma once
#ifndef __MATLABDEBUG__
#define __MATLABDEBUG__

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Sparse>
#include <Eigen/Dense>

//template <typename Derived>

//void sparse_to_file_as_dense(const Eigen::EigenBase<Derived>& mat, std::string var_name);

void sparse_to_file_as_sparse_m(const Eigen::SparseMatrix<double>& mat, std::string var_name);

void sparse_to_file_as_dense(const Eigen::SparseMatrix<double>& mat, std::string var_name);

void mat_to_file(const Eigen::MatrixXd& mat, std::string var_name);

void mat_to_file(const Eigen::Matrix4d& mat, std::string var_name);

void mat_to_file(const Eigen::MatrixXi& mat, std::string var_name);

void vec_to_file(const Eigen::VectorXd& vec, std::string var_name);

void vec_to_file(const Eigen::VectorXi& vec, std::string var_name);

void vec_to_file(const Eigen::Vector3d& vec, std::string var_name);

void double_to_file(double d, std::string var_name);

void EOL_outputter(int EOLE, int EOLV, int verts, int faces);

//#include "MatlabDebugINL.cpp"

#endif
