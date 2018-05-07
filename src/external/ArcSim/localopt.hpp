#ifndef LOCALOPT_HPP
#define LOCALOPT_HPP

#include "geometry.hpp"
#include "constraint.hpp"

template<Space s>
void local_opt(std::vector<Node*>& nodes, std::vector<Face*>& faces, std::vector<Edge*>& edges,
               const std::vector<Constraint*>& cons);

#endif
