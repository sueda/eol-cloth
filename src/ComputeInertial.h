#pragma once
#ifndef __ComputeInertial__
#define __ComputeInertial__

void ComputeInertial(
	const double *xa, // [input 3x1] World position of vertex A
	const double *xb, // [input 3x1] World position of vertex B
	const double *xc, // [input 3x1] World position of vertex C
	const double *Xa, // [input 2x1] Material position of vertex A
	const double *Xb, // [input 2x1] Material position of vertex B
	const double *Xc, // [input 2x1] Material position of vertex C
	const double *g,  // [input 3x1] 3D gravity vector
	double rho,       // [input 1x1] Density (mass per area)
	double *W,        // [output 1x1] Gravitational potential energy
	double *f,        // [output 9x1] Gravity force vector
	double *M)        // [output 9x9] Inertia matrix
	;

#endif