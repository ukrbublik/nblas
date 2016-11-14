#include "blas_sparse.h"
#include "routines.h"


// https://software.intel.com/ru-ru/node/468508
void dusdot(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int nz = info[0]->Int32Value();
	const double *x = reinterpret_cast<double*>(GET_CONTENTS(info[1].As<v8::Float64Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[2].As<v8::Int32Array>()));
	const double *y = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	const int incy = info[4]->Int32Value();
	double res = 0;
	BLAS_dusdot(blas_no_conj, nz, x, indx, y, incy, &res, blas_zero_base);
	info.GetReturnValue().Set(res);
}
void susdot(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int nz = info[0]->Int32Value();
	const float *x = reinterpret_cast<float*>(GET_CONTENTS(info[1].As<v8::Float32Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[2].As<v8::Int32Array>()));
	const float *y = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	const int incy = info[4]->Int32Value();
	float res = 0;
	BLAS_susdot(blas_no_conj, nz, x, indx, y, incy, &res, blas_zero_base);
	info.GetReturnValue().Set(res);
}

//https://software.intel.com/ru-ru/node/468508
void dusaxpy(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int nz = info[0]->Int32Value();
	const double alpha = info[1]->NumberValue();
	const double *x = reinterpret_cast<double*>(GET_CONTENTS(info[2].As<v8::Float64Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[3].As<v8::Int32Array>()));
	double *y = reinterpret_cast<double*>(GET_CONTENTS(info[4].As<v8::Float64Array>()));
	const int incy = info[5]->Int32Value();
	BLAS_dusaxpy(nz, alpha, x, indx, y, incy, blas_zero_base);
}
void susaxpy(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int nz = info[0]->Int32Value();
	const float alpha = info[1]->NumberValue();
	const float *x = reinterpret_cast<float*>(GET_CONTENTS(info[2].As<v8::Float32Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[3].As<v8::Int32Array>()));
	float *y = reinterpret_cast<float*>(GET_CONTENTS(info[4].As<v8::Float32Array>()));
	const int incy = info[5]->Int32Value();
	BLAS_susaxpy(nz, alpha, x, indx, y, incy, blas_zero_base);
}

//https://software.intel.com/ru-ru/node/468516
void dusga(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int nz = info[0]->Int32Value();
	const double *y = reinterpret_cast<double*>(GET_CONTENTS(info[1].As<v8::Float64Array>()));
	const int incy = info[2]->Int32Value();
	double *x = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	BLAS_dusga(nz, y, incy, x, indx, blas_zero_base);
}
void susga(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int nz = info[0]->Int32Value();
	const float *y = reinterpret_cast<float*>(GET_CONTENTS(info[1].As<v8::Float32Array>()));
	const int incy = info[2]->Int32Value();
	float *x = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	BLAS_susga(nz, y, incy, x, indx, blas_zero_base);
}

//https://software.intel.com/ru-ru/node/468518
void dusgz(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int nz = info[0]->Int32Value();
	double *y = reinterpret_cast<double*>(GET_CONTENTS(info[1].As<v8::Float64Array>()));
	const int incy = info[2]->Int32Value();
	double *x = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	BLAS_dusgz(nz, y, incy, x, indx, blas_zero_base);
}
void susgz(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int nz = info[0]->Int32Value();
	float *y = reinterpret_cast<float*>(GET_CONTENTS(info[1].As<v8::Float32Array>()));
	const int incy = info[2]->Int32Value();
	float *x = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	BLAS_susgz(nz, y, incy, x, indx, blas_zero_base);
}

//https://software.intel.com/ru-ru/node/468522
void dussc(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int nz = info[0]->Int32Value();
	const double *x = reinterpret_cast<double*>(GET_CONTENTS(info[1].As<v8::Float64Array>()));
	double *y = reinterpret_cast<double*>(GET_CONTENTS(info[2].As<v8::Float64Array>()));
	const int incy = info[3]->Int32Value();
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	BLAS_dussc(nz, x, y, incy, indx, blas_zero_base);
}
void sussc(const v8::FunctionCallbackInfo<v8::Value>& info) {
	const int nz = info[0]->Int32Value();
	const float *x = reinterpret_cast<float*>(GET_CONTENTS(info[1].As<v8::Float32Array>()));
	float *y = reinterpret_cast<float*>(GET_CONTENTS(info[2].As<v8::Float32Array>()));
	const int incy = info[3]->Int32Value();
	const int *indx = reinterpret_cast<int*>(GET_CONTENTS(info[4].As<v8::Int32Array>()));
	BLAS_sussc(nz, x, y, incy, indx, blas_zero_base);
}
