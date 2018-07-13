#pragma once
#ifndef __matlabOutputs__
#define __matlabOutputs__

#include "external\ArcSim\mesh.hpp"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen\Dense>
#include <Eigen\Sparse>

#include <string>

void mat_s2s_file(const Eigen::SparseMatrix<double>& mat, const std::string &var_name, const std::string &file_name, const bool& overwrite);

void vec_to_file(const Eigen::VectorXd& vec, const std::string &var_name, const std::string &file_name, const bool &overwrite);

void mesh2m(const Mesh& mesh, const std::string &file_name, const bool &overwrite);

#endif