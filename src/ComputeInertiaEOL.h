#pragma once
#ifndef __ComputeInertialEOL__
#define __ComputeInertialEOL__

void ComputeInertiaEOL(
	const double *xa,
	const double *xb,
	const double *xc,
	const double *Xa,
	const double *Xb,
	const double *Xc,
	const double *g,
	double rho,
	double *W,
	double *f,
	double *M)
	;

#endif