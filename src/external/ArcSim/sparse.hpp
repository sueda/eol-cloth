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

#ifndef SPARSE_HPP
#define SPARSE_HPP

#include "util.hpp"
#include <fstream>
#include <utility>
#include <vector>

inline size_t find_index (int i, const std::vector<int> &indices) {
    for (size_t ii = 0; ii < indices.size(); ii++)
        if (indices[ii] == i)
            return ii;
    return indices.size();
}

template <typename T> void insert_index (int i, int j,
                                         std::vector<int> &indices,
                                         std::vector<T> &entries) {
    indices.insert(indices.begin() + j, i);
    entries.insert(entries.begin() + j, T(0));
}

template <typename T> struct SpVec {
    std::vector<int> indices;
    std::vector<T> entries;
    T operator[] (int i) const {
        size_t j = find_index(i, indices);
        if (j >= indices.size() || indices[j] != i)
            return T(0);
        else
            return entries[j];
    }
    T &operator[] (int i) {// inserts entry as side-effect
        size_t j = find_index(i, indices);
        if (j >= indices.size() || indices[j] != i)
            insert_index((int)i, (int)j, indices, entries);
        return entries[j];
    }
};

template <typename T>
std::ostream &operator<< (std::ostream &out, const SpVec<T> &v) {
    out << "[";
    for (int k = 0; k < v.indices.size(); k++)
        out << (k==0 ? "" : ", ") << v.indices[k] << ": " << v.entries[k];
    out << "]";
    return out;
}

template <typename T> struct SpMat {
    int m, n;
    std::vector< SpVec<T> > rows;
    SpMat (): m(0), n(0), rows() {}
    explicit SpMat (int m, int n): m(m), n(n), rows(m) {}
    T operator() (int i, int j) const {
        return rows[i][j];
    }
    T &operator() (int i, int j) {// inserts entry as side-effect
        return rows[i][j];
    }
};

template <typename T>
std::ostream &operator<< (std::ostream &out, const SpMat<T> &A) {
    out << "[";
    for (int i = 0; i < A.m; i++) {
        const SpVec<T> &row = A.rows[i];
        for (int jj = 0; jj < row.indices.size(); jj++) {
            int j = row.indices[jj];
            const T &aij = row.entries[jj];
            out << (i==0 && jj==0 ? "" : ", ") << "(" << i << "," << j
                << "): " << aij;
        }
    }
    out << "]";
    return out;
}

#endif
