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

#ifndef WINPORT_HPP
#define WINPORT_HPP
// MS Windows bindings, etc

#if defined(_WIN32) && !defined(__CYGWIN__)

#pragma warning(disable:4018) // signed/unsigned mismatch
#pragma warning(disable:4244) // conversion from 'double' to 'float', possible loss of data
#pragma warning(disable:4996) // this function or variable may be unsafe
#pragma warning(disable:4251) // class needs to have dll-interface to be used by clients
#pragma warning(disable:4800) // forcing value to bool 'true' or 'false'
#pragma warning(disable:161)  // unrecognized #pragma

#define _USE_MATH_DEFINES // just to have M_PI
#include <cmath>
#include <math.h>
#include <algorithm>

#include <windows.h>
#undef min
#undef max
#include <stdio.h>
#include <iostream>
#define snprintf _snprintf

//#include <boost/math/special_functions/fpclassify.hpp> 
//template <class T> inline bool isfinite(const T& number) { return boost::math::isfinite(number); }
//template <class T> inline bool   finite(const T& number) { return boost::math::isfinite(number); }

inline double sqrt(int n) { return sqrt(double(n)); }

template <class T> inline T log2(const T& number) { return log(number)/log(T(2)); }

extern std::ostream cdbg;

#endif

#endif
