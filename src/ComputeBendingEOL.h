#pragma once
#ifndef __ComputeBendingEOL__
#define __ComputeBendingEOL__

void ComputeBendingEOL(
	const double *x0,
	const double *x1,
	const double *x2,
	const double *x3,
	const double *X0,
	const double *X1,
	const double *X2,
	const double *X3,
	double beta,
	double *W,
	double *f,
	double *K)
	;

#endif