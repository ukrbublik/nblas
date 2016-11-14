#include "blas_sparse.h"
#include "routines.h"


// Creation Routines
void duscr_begin(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int m = info[0]->Int32Value();
	const int n = info[1]->Int32Value();
	blas_sparse_matrix spm = BLAS_duscr_begin(m, n);
	info.GetReturnValue().Set(spm);
}
void suscr_begin(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int m = info[0]->Int32Value();
	const int n = info[1]->Int32Value();
	blas_sparse_matrix spm = BLAS_suscr_begin(m, n);
	info.GetReturnValue().Set(spm);
}

void duscr_block_begin(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int Mb = info[0]->Int32Value();
	const int Nb = info[1]->Int32Value();
	const int k = info[2]->Int32Value();
	const int l = info[3]->Int32Value();
	blas_sparse_matrix A = BLAS_duscr_block_begin(Mb, Nb, k, l);
	info.GetReturnValue().Set(A);
}
void suscr_block_begin(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int Mb = info[0]->Int32Value();
	const int Nb = info[1]->Int32Value();
	const int k = info[2]->Int32Value();
	const int l = info[3]->Int32Value();
	blas_sparse_matrix A = BLAS_suscr_block_begin(Mb, Nb, k, l);
	info.GetReturnValue().Set(A);
}

void duscr_variable_block_begin(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int Mb = info[0]->Int32Value();
	const int Nb = info[1]->Int32Value();
	const int *k = reinterpret_cast<int*>(GET_CONTENTS(info[2].As<v8::Int32Array>()));
	const int *l = reinterpret_cast<int*>(GET_CONTENTS(info[3].As<v8::Int32Array>()));
	blas_sparse_matrix A = BLAS_duscr_variable_block_begin(Mb, Nb, k, l);
	info.GetReturnValue().Set(A);
}
void suscr_variable_block_begin(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int Mb = info[0]->Int32Value();
	const int Nb = info[1]->Int32Value();
	const int *k = reinterpret_cast<int*>(GET_CONTENTS(info[2].As<v8::Int32Array>()));
	const int *l = reinterpret_cast<int*>(GET_CONTENTS(info[3].As<v8::Int32Array>()));
	blas_sparse_matrix A = BLAS_suscr_variable_block_begin(Mb, Nb, k, l);
	info.GetReturnValue().Set(A);
}


// Insertion Routines
void duscr_insert_entry(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const double val = info[1]->NumberValue();
	const int i = info[2]->Int32Value();
	const int j = info[3]->Int32Value();
	int res = BLAS_duscr_insert_entry(A, val, i, j);
	info.GetReturnValue().Set(res);
}
void suscr_insert_entry(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const float val = info[1]->NumberValue();
	const int i = info[2]->Int32Value();
	const int j = info[3]->Int32Value();
	int res = BLAS_suscr_insert_entry(A, val, i, j);
	info.GetReturnValue().Set(res);
}

void duscr_insert_entries(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const int nz = info[1]->Int32Value();
	const double *val = reinterpret_cast<double*>(GET_CONTENTS(info[2].As<v8::Float64Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[3].As<v8::Int32Array>()));
	const int *jndx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	int res = BLAS_duscr_insert_entries(A, nz, val, indx, jndx);
	info.GetReturnValue().Set(res);
}
void suscr_insert_entries(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const int nz = info[1]->Int32Value();
	const float *val = reinterpret_cast<float*>(GET_CONTENTS(info[2].As<v8::Float32Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[3].As<v8::Int32Array>()));
	const int *jndx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	int res = BLAS_suscr_insert_entries(A, nz, val, indx, jndx);
	info.GetReturnValue().Set(res);
}

void duscr_insert_col(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const int j = info[1]->Int32Value();
	const int nz = info[2]->Int32Value();
	const double *val = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	int res = BLAS_duscr_insert_col(A, j, nz, val, indx);
	info.GetReturnValue().Set(res);
}
void suscr_insert_col(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const int j = info[1]->Int32Value();
	const int nz = info[2]->Int32Value();
	const float *val = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	int res = BLAS_suscr_insert_col(A, j, nz, val, indx);
	info.GetReturnValue().Set(res);
}

void duscr_insert_row(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const int i = info[1]->Int32Value();
	const int nz = info[2]->Int32Value();
	const double *val = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	const int *jndx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	int res = BLAS_duscr_insert_row(A, i, nz, val, jndx);
	info.GetReturnValue().Set(res);
}
void suscr_insert_row(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const int i = info[1]->Int32Value();
	const int nz = info[2]->Int32Value();
	const float *val = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	const int *jndx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	int res = BLAS_suscr_insert_row(A, i, nz, val, jndx);
	info.GetReturnValue().Set(res);
}

void duscr_insert_clique(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const int k = info[1]->Int32Value();
	const int l = info[2]->Int32Value();
	const double *val = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	const int row_stride = info[4]->Int32Value();
	const int col_stride = info[5]->Int32Value();
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[6].As<v8::Int32Array>()));
	const int *jndx = reinterpret_cast<int*>(GET_CONTENTS(info[7].As<v8::Int32Array>()));
	int res = BLAS_duscr_insert_clique(A, k, l, val, row_stride, col_stride, indx, jndx);
	info.GetReturnValue().Set(res);
}
void suscr_insert_clique(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const int k = info[1]->Int32Value();
	const int l = info[2]->Int32Value();
	const float *val = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	const int row_stride = info[4]->Int32Value();
	const int col_stride = info[5]->Int32Value();
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[6].As<v8::Int32Array>()));
	const int *jndx = reinterpret_cast<int*>(GET_CONTENTS(info[7].As<v8::Int32Array>()));
	int res = BLAS_suscr_insert_clique(A, k, l, val, row_stride, col_stride, indx, jndx);
	info.GetReturnValue().Set(res);
}

void duscr_insert_block(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const double *val = reinterpret_cast<double*>(GET_CONTENTS(info[1].As<v8::Float64Array>()));
	const int row_stride = info[2]->Int32Value();
	const int col_stride = info[3]->Int32Value();
	const int i = info[4]->Int32Value();
	const int j = info[5]->Int32Value();
	int res = BLAS_duscr_insert_block(A, val, row_stride, col_stride, i, j);
	info.GetReturnValue().Set(res);
}
void suscr_insert_block(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const float *val = reinterpret_cast<float*>(GET_CONTENTS(info[1].As<v8::Float32Array>()));
	const int row_stride = info[2]->Int32Value();
	const int col_stride = info[3]->Int32Value();
	const int i = info[4]->Int32Value();
	const int j = info[5]->Int32Value();
	int res = BLAS_suscr_insert_block(A, val, row_stride, col_stride, i, j);
	info.GetReturnValue().Set(res);
}


// Completion of Construction Routines
void duscr_end(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	int res = BLAS_duscr_end(A);
	info.GetReturnValue().Set(res);
}
void suscr_end(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	int res = BLAS_suscr_end(A);
	info.GetReturnValue().Set(res);
}

// Matrix Property Routines
void _usgp(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const int pname = info[1]->Int32Value();
	int res = BLAS_usgp(A, pname);
	info.GetReturnValue().Set(res);
}
void _ussp(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	const int pname = info[1]->Int32Value();
	int res = BLAS_ussp(A, pname);
	info.GetReturnValue().Set(res);
}

// Destruction Routine
void _usds(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int A = info[0]->Int32Value();
	int res = BLAS_usds(A);
	info.GetReturnValue().Set(res);
}
