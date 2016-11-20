#include "cblas.h"
#include "routines.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>

// a [N*N] * x [N*NRHS] = b [N*NRHS]
// solution out in b
// order is col-based by default (because of Fortran)

//http://physics.oregonstate.edu/~landaur/nacphy/lapack/routines/dgesv.html
void dgesv(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int n = info[0]->Int32Value();
	int nrhs = info[1]->Int32Value();
	double *a = reinterpret_cast<double*>(GET_CONTENTS(info[2].As<v8::Float64Array>()));
	int lda = n;
	double *b = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	//ipiv size = (N)
	int *ipiv = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	int ldb = n;
	int res;
	int matrix_layout = LAPACK_ROW_MAJOR;
	if (matrix_layout == LAPACK_ROW_MAJOR) {
		//http://stackoverflow.com/questions/30212743/using-lapacke-zgetrs-with-lapack-row-major-causes-illegal-memory-access
		ldb = nrhs;
	}
	res = LAPACKE_dgesv( matrix_layout, n, nrhs, a, lda, ipiv, b, ldb );
	info.GetReturnValue().Set(res);
}

void sgesv(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int n = info[0]->Int32Value();
	int nrhs = info[1]->Int32Value();
	float *a = reinterpret_cast<float*>(GET_CONTENTS(info[2].As<v8::Float32Array>()));
	int lda = n;
	float *b = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	//ipiv size = (N)
	int *ipiv = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	int ldb = n;
	int res;
	int matrix_layout = LAPACK_ROW_MAJOR;
	if (matrix_layout == LAPACK_ROW_MAJOR) {
		ldb = nrhs;
	}
	res = LAPACKE_sgesv( matrix_layout, n, nrhs, a, lda, ipiv, b, ldb );
	info.GetReturnValue().Set(res);
}

