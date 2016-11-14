#include "cblas.h"
#include "routines.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include <lapacke.h>
using std::swap;

#if !defined(LAPACK_ROW_MAJOR)
	#define LAPACK_ROW_MAJOR               101
	#define LAPACK_COL_MAJOR               102
#endif

extern int dgesv_(int matrix_layout, int n, int nrhs, double *a, int lda, int* ipiv, double *b, int ldb);
extern int sgesv_(int matrix_layout, int n, int nrhs, float *a, int lda, int* ipiv, float *b, int ldb);

// a [N*N] * x [N*NRHS] = b [N*NRHS]
// solution out in b
// order is col-based by default (because of Fortran)

void dgesv(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int n = info[0]->Int32Value();
	int nrhs = info[1]->Int32Value();
	double *a = reinterpret_cast<double*>(GET_CONTENTS(info[2].As<v8::Float64Array>()));
	int lda = n;
	int ipiv[n];
	double *b = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	int ldb = n;
	int matrix_layout = LAPACK_ROW_MAJOR;
	if (matrix_layout == LAPACK_ROW_MAJOR) {
		//http://stackoverflow.com/questions/30212743/using-lapacke-zgetrs-with-lapack-row-major-causes-illegal-memory-access
		ldb = nrhs;
	}
	//int res = LAPACKE_dgesv( matrix_layout, n, nrhs, a, lda, ipiv, b, ldb );
	int res = dgesv_( matrix_layout, n, nrhs, a, lda, ipiv, b, ldb );
	info.GetReturnValue().Set(res);
}

void sgesv(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int n = info[0]->Int32Value();
	int nrhs = info[1]->Int32Value();
	float *a = reinterpret_cast<float*>(GET_CONTENTS(info[2].As<v8::Float32Array>()));
	int lda = n;
	int ipiv[n];
	float *b = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	int ldb = n;
	int matrix_layout = LAPACK_ROW_MAJOR;
	if (matrix_layout == LAPACK_ROW_MAJOR) {
		ldb = nrhs;
	}
	//int res = LAPACKE_sgesv( matrix_layout, n, nrhs, a, lda, ipiv, b, ldb );
	int res = sgesv_( matrix_layout, n, nrhs, a, lda, ipiv, b, ldb );
	info.GetReturnValue().Set(res);
}
