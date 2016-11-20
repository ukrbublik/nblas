#include "cblas.h"
#include "routines.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>

//http://www.netlib.org/clapack/clapack-3.2.1-CMAKE/SRC/VARIANTS/lu/CR/dgetrf.c
void dgetrf(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int m = info[0]->Int32Value();
	int n = info[1]->Int32Value();
	double *a = reinterpret_cast<double*>(GET_CONTENTS(info[2].As<v8::Float64Array>()));
	int lda = m;
	int *ipiv = reinterpret_cast<int*>(GET_CONTENTS(info[3].As<v8::Int32Array>()));
	//ipiv size = (min(M,N))
	int res;
	int matrix_layout = LAPACK_ROW_MAJOR;
	if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
	}
	res = LAPACKE_dgetrf( matrix_layout, m, n, a, lda, ipiv );
	info.GetReturnValue().Set(res);
}

void sgetrf(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int m = info[0]->Int32Value();
	int n = info[1]->Int32Value();
	float *a = reinterpret_cast<float*>(GET_CONTENTS(info[2].As<v8::Float32Array>()));
	int lda = m;
	int *ipiv = reinterpret_cast<int*>(GET_CONTENTS(info[3].As<v8::Int32Array>()));
	//ipiv size = (min(M,N))
	int res;
	int matrix_layout = LAPACK_ROW_MAJOR;
	if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
	}
	res = LAPACKE_sgetrf( matrix_layout, m, n, a, lda, ipiv );
	info.GetReturnValue().Set(res);
}
