#include "blas_sparse.h"
#include "routines.h"

/**
*
*   Sparse matrix times dense vector.  This routine computes
*   y = alpha * op(A) * x + beta * y, where op(A) is either A,
*   the transpose of A, and x and y are dense vectors.
*
*   The routine returns 0 if sucessful, -1 otherwise. If A
*   is an invalid handle, the routine returns -1 and does not
*   modify any of its parameters.
*
*
*   @param trans    compute operation with the transpose of
*   A, if trans is set to blas_transpose.  Otherwise, it computes
*   with A untransposed.   (Note: trans set to blas_conj_trans has
*   the same effect as blas_trans for real valued matrices.)
*
*   @param alpha scalar multiplier of A
*   @param A     (input) sparse matrix handle
*   @param x     (input) dense vector,
*   @param incx  increment between elements of x
*   @param beta  scalar multiplier of y
*   @param y     (input/output) dense vector, one output contains
*                   alpha * op(A) * x + beta * y.
*   @param incy  increment (stride) between elemens of y
*
*/
// A [m * n] * x [ n * 1 ] = y [ m * 1 ]
void dusmv(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const enum blas_trans_type transA = static_cast<blas_trans_type>(info[0]->Uint32Value());
	const double alpha = info[1]->NumberValue();
	const int A = info[2]->Int32Value();
	const double *x = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	const int incx = info[4]->Int32Value();
	double *y = reinterpret_cast<double*>(GET_CONTENTS(info[5].As<v8::Float64Array>()));
	int incy = info[6]->Int32Value();
	int res = BLAS_dusmv(transA, alpha, A, x, incx, y, incy);
	info.GetReturnValue().Set(res);
}
void susmv(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const enum blas_trans_type transA = static_cast<blas_trans_type>(info[0]->Uint32Value());
	const double alpha = info[1]->NumberValue();
	const int A = info[2]->Int32Value();
	const float *x = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	const int incx = info[4]->Int32Value();
	float *y = reinterpret_cast<float*>(GET_CONTENTS(info[5].As<v8::Float32Array>()));
	int incy = info[6]->Int32Value();
	int res = BLAS_susmv(transA, alpha, A, x, incx, y, incy);
	info.GetReturnValue().Set(res);
}


/**
*
*   Sparse matrix triangular solve.  Given a sparse triangular
*   matrix A and scalar alpha, this routine computes
*   x = alpha * op(A)^{-1}x.
*
*   The routine returns 0 if sucessful, -1 otherwise.
*   The matrix handle A must be refer to an upper or lower triangular
*   matrix, otherwise the routine returns -1,  and does not
*   modify any of its parameters.
*
*   @param trans    computes operation with the transpose of
*   A, if trans is set to blas_transpose.  Otherwise, it computes
*   with A untransposed.   (Note: trans set to blas_conj_trans has
*   the same effect as blas_trans for real valued matrices.)
* 
*   @param alpha scalar multiplier of A.  
*
*   @param T     (input) handle representing a triangular sparse matrix.  
*       Note that T must have been created with the blas_triangular property 
*       set, as well as either the blas_upper or blas_lower property to 
*       denote upper or lower triangular structure.
*
*   @param x     (input/output) on input, it contains the original
*                   right-hand-side of linear equation to be solved.
*               On output, it contains the solution.
*   @param incx  stride between entries of x
*
*/
void dussv(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const enum blas_trans_type transA = static_cast<blas_trans_type>(info[0]->Uint32Value());
	const double alpha = info[1]->NumberValue();
	const int A = info[2]->Int32Value();
	double *x = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	const int incx = info[4]->Int32Value();
	int res = BLAS_dussv(transA, alpha, A, x, incx);
	info.GetReturnValue().Set(res);
}
void sussv(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const enum blas_trans_type transA = static_cast<blas_trans_type>(info[0]->Uint32Value());
	const double alpha = info[1]->NumberValue();
	const int A = info[2]->Int32Value();
	float *x = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	const int incx = info[4]->Int32Value();
	int res = BLAS_sussv(transA, alpha, A, x, incx);
	info.GetReturnValue().Set(res);
}