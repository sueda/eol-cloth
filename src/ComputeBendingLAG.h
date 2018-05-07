#pragma once
#ifndef __ComputeBendingLAG__
#define __ComputeBendingLAG__

void ComputeBendingLAG(
	const double *x0, // [input 3x1] World position of vertex 0
	const double *x1, // [input 3x1] World position of vertex 1
	const double *x2, // [input 3x1] World position of vertex 2
	const double *x3, // [input 3x1] World position of vertex 3
	const double *X0, // [input 2x1] Material position of vertex 0
	const double *X1, // [input 2x1] Material position of vertex 1
	const double *X2, // [input 2x1] Material position of vertex 2
	const double *X3, // [input 2x1] Material position of vertex 3
	double beta,      // [input 1x1] Bending stiffness 
	double *W,        // [output 1x1] Bending potential energy
	double *f,        // [output 12x1] Bending force vector
	double *K)        // [output 12x12] Bending stiffness matrix
	;

#endif