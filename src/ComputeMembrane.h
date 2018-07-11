#pragma once
#ifndef __ComputeMembrane__
#define __ComputeMembrane__

void ComputeMembrane(
	const double *xa, // [input 3x1] World position of vertex A
	const double *xb, // [input 3x1] World position of vertex B
	const double *xc, // [input 3x1] World position of vertex C
	const double *Xa, // [input 2x1] Material position of vertex A
	const double *Xb, // [input 2x1] Material position of vertex B
	const double *Xc, // [input 2x1] Material position of vertex C
	double e,         // [input 1x1] Young's modulus -pascals
	double nu,        // [input 1x1] Poisson's ratio 0.3 or 0.2
	const double *P,  // [input 2x3] Projection matrix (column major)
	const double *Q,  // [input 2x2] Polar decomposition of \bar{F} (column major)
	double *W,        // [output 1x1] Membrane potential energy
	double *f,        // [output 9x1] Membrane force vector
	double *K)        // [output 9x9] Membrane stiffness matrix
	;

#endif