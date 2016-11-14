#include "blas_sparse.h"
#include "routines.h"


/**

   Sparse matrix triangular solve.  Given a sparse triangular
   matrix A, dense rectangular matrix B and scalar alpha, this routine 
   computes
   B = alpha * op(T)^{-1}B.

   The routine returns 0 if sucessful, -1 otherwise.
   The matrix handle A must be refer to an upper or lower triangular
   matrix, otherwise the routine returns -1,  and does not
   modify any of its parameters.

   @param trans    computes operation with the transpose of
   A, if trans is set to blas_transpose.  Otherwise, it computes
   with A untransposed.   (Note: trans set to blas_conj_trans has
   the same effect as blas_trans for real valued matrices.)
 
   @param alpha scalar multiplier of A.  

   @param A     (input) handle representing a triangular sparse matrix.  
       Note that A must have been created with the blas_triangular property 
       set, as well as either the blas_upper or blas_lower property to 
       denote upper or lower triangular structure.

   @param B     (input/output) on input, it contains the original
                   right-hand-side of linear equation to be solved.
               On output, it contains the solution.
   @param ldB	stride between outer indices of B. (i.e. between rows, if 
   					B is column ordered, or between colums if 
					B is row-ordered.)

*/
void dussm(const v8::FunctionCallbackInfo<v8::Value>& info) {
  enum blas_order_type order = blas_rowmajor;
  const enum blas_trans_type transt = static_cast<blas_trans_type>(info[0]->Uint32Value());
  const int nrhs = info[1]->Int32Value();
  const double alpha = info[2]->NumberValue();
  const int A = info[3]->Int32Value();
  double *B = reinterpret_cast<double*>(GET_CONTENTS(info[4].As<v8::Float64Array>()));
  //int m = BLAS_usgp(A, blas_num_rows);
  int n = BLAS_usgp(A, blas_num_cols);
  int ldb = n;
  if (order == blas_rowmajor) {
    ldb = nrhs;
  }
  int res = BLAS_dussm(order, transt, nrhs, alpha, A, B, ldb);
  info.GetReturnValue().Set(res);
}
void sussm(const v8::FunctionCallbackInfo<v8::Value>& info) {
  enum blas_order_type order = blas_rowmajor;
  const enum blas_trans_type transt = static_cast<blas_trans_type>(info[0]->Uint32Value());
  const int nrhs = info[1]->Int32Value();
  const float alpha = info[2]->NumberValue();
  const int A = info[3]->Int32Value();
  float *B = reinterpret_cast<float*>(GET_CONTENTS(info[4].As<v8::Float32Array>()));
  //int m = BLAS_usgp(A, blas_num_rows);
  int n = BLAS_usgp(A, blas_num_cols);
  int ldb = n;
  if (order == blas_rowmajor) {
    ldb = nrhs;
  }
  int res = BLAS_sussm(order, transt, nrhs, alpha, A, B, ldb);
  info.GetReturnValue().Set(res);
}


// a [M*N] * b [N*NRHS] + c = c [M*NRHS]
// B and C are dense matrices, A is sparse matrix, out in C
void dusmm(const v8::FunctionCallbackInfo<v8::Value>& info) {
  enum blas_order_type order = blas_rowmajor;
  const enum blas_trans_type transa = static_cast<blas_trans_type>(info[0]->Uint32Value());
  const int nrhs = info[1]->Int32Value();
  const double alpha = info[2]->NumberValue();
  const int A = info[3]->Int32Value();
  double *B = reinterpret_cast<double*>(GET_CONTENTS(info[4].As<v8::Float64Array>()));
  double *C = reinterpret_cast<double*>(GET_CONTENTS(info[5].As<v8::Float64Array>()));
  int m = BLAS_usgp(A, blas_num_rows);
  int n = BLAS_usgp(A, blas_num_cols);
  int ldb = n;
  int ldc = m;
  if (order == blas_rowmajor) {
    ldb = nrhs;
    ldc = nrhs;
  }
  int res = BLAS_dusmm(order, transa, nrhs, alpha, A, B, ldb, C, ldc);
  info.GetReturnValue().Set(res);
}
void susmm(const v8::FunctionCallbackInfo<v8::Value>& info) {
  enum blas_order_type order = blas_rowmajor;
  const enum blas_trans_type transa = static_cast<blas_trans_type>(info[0]->Uint32Value());
  const int nrhs = info[1]->Int32Value();
  const float alpha = info[2]->NumberValue();
  const int A = info[3]->Int32Value();
  float *B = reinterpret_cast<float*>(GET_CONTENTS(info[4].As<v8::Float32Array>()));
  float *C = reinterpret_cast<float*>(GET_CONTENTS(info[5].As<v8::Float32Array>()));
  int m = BLAS_usgp(A, blas_num_rows);
  int n = BLAS_usgp(A, blas_num_cols);
  int ldb = n;
  int ldc = m;
  if (order == blas_rowmajor) {
    ldb = nrhs;
    ldc = nrhs;
  }
  int res = BLAS_susmm(order, transa, nrhs, alpha, A, B, ldb, C, ldc);
  info.GetReturnValue().Set(res);
}
