#include "matlabOutputs.h"

#include <fstream>
#include <iomanip>
#include <Eigen\Dense>

using namespace std;
using namespace Eigen;

void mat_to_file(const Eigen::MatrixXd& mat, const string &var_name, const string &file_name, const bool &overwrite)
{
	ofstream ofs;
	if(overwrite) ofs.open(file_name, ofstream::out | ofstream::trunc);
	else ofs.open(file_name, ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = [";
	ofs << mat;
	ofs << "];\n\n";

	ofs.close();
}

void mat_to_file(const Eigen::MatrixXi& mat, const string &var_name, const string &file_name, const bool &overwrite)
{
	ofstream ofs;
	if (overwrite) ofs.open(file_name, ofstream::out | ofstream::trunc);
	else ofs.open(file_name, ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = [";
	ofs << mat;
	ofs << "];\n\n";

	ofs.close();
}

void vec_to_file(const Eigen::VectorXd& vec, const string &var_name, const string &file_name, const bool &overwrite)
{
	ofstream ofs;
	if (overwrite) ofs.open(file_name, ofstream::out | ofstream::trunc);
	else ofs.open(file_name, ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = [\n";
	ofs << vec;
	ofs << "\n];\n\n";

	ofs.close();
}

void vec_to_file(const Eigen::VectorXi& vec, const string &var_name, const string &file_name, const bool &overwrite)
{
	ofstream ofs;
	if (overwrite) ofs.open(file_name, ofstream::out | ofstream::trunc);
	else ofs.open(file_name, ofstream::out | ofstream::app);

	ofs << setprecision(20);

	ofs << var_name;
	ofs << " = [\n";
	ofs << vec;
	ofs << "\n];\n\n";

	ofs.close();
}

void mesh2m(const Mesh& mesh, const string &file_name, const bool &overwrite)
{
	MatrixXd x_X(mesh.nodes.size(), 5);
	MatrixXi faces2(3, mesh.faces.size());
	VectorXi isEOL(mesh.nodes.size());
	for (int i = 0; i < mesh.nodes.size(); i++) {

		isEOL(i) = 0;
		x_X(i, 0) = mesh.nodes[i]->x[0];
		x_X(i, 1) = mesh.nodes[i]->x[1];
		x_X(i, 2) = mesh.nodes[i]->x[2];
		x_X(i, 3) = mesh.nodes[i]->verts[0]->u[0];
		x_X(i, 4) = mesh.nodes[i]->verts[0]->u[1];
	}

	for (int i = 0; i < mesh.faces.size(); i++) {
		faces2.col(i) = Vector3i(mesh.faces[i]->v[0]->node->index, mesh.faces[i]->v[1]->node->index, mesh.faces[i]->v[2]->node->index);
	}

	mat_to_file(x_X, "x_X", file_name, overwrite);
	vec_to_file(isEOL, "isEol", file_name, false);
	VectorXi vvv(3);
	vvv << 1, 1, 1;
	mat_to_file(faces2.colwise() += vvv, "faces", file_name, false);
}