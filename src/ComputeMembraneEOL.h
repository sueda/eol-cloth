#pragma once
#ifndef __ComputeMembraneEOL__
#define __ComputeMembraneEOL__

void ComputeMembraneEOL(
	const double *xa,
	const double *xb,
	const double *xc,
	const double *Xa,
	const double *Xb,
	const double *Xc,
	double e,
	double nu,
	const double *P,
	const double *Q,
	double *W,
	double *f,
	double *K)
	;

#endif