/*
Copyright ©2013 The Regents of the University of California
(Regents). All Rights Reserved. Permission to use, copy, modify, and
distribute this software and its documentation for educational,
research, and not-for-profit purposes, without fee and without a
signed licensing agreement, is hereby granted, provided that the
above copyright notice, this paragraph and the following two
paragraphs appear in all copies, modifications, and
distributions. Contact The Office of Technology Licensing, UC
Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
(510) 643-7201, for commercial licensing opportunities.

IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT,
INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS
IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#ifndef IO_HPP
#define IO_HPP

#include "mesh.hpp"
#include "util.hpp"
//#include <zlib.h>

static double angle(const Vec3 &x0, const Vec3 &x1, const Vec3 &x2);
std::vector<Face*> triangulateARC(const std::vector<Vert*> &verts);

void triangle_to_obj(const std::string &infile, const std::string &outfile);

void load_obj(Mesh &mesh, const std::string &filename);
void load_objs(std::vector<Mesh*> &meshes, const std::string &prefix);

void load_cloth_obj(Mesh &mesh, const std::string &filename);


#endif
