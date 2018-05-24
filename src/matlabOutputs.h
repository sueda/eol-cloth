#pragma once
#ifndef __matlabOutputs__
#define __matlabOutputs__

#include "external\ArcSim\mesh.hpp"

#include <string>

void mesh2m(const Mesh& mesh, const std::string &file_name, const bool &overwrite);

#endif