#include "routines.h"
#include <string.h>

// transpose a -> b
void dTrTo(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int r = info[0]->Int32Value();
	int c = info[1]->Int32Value();
	double *a = reinterpret_cast<double*>(GET_CONTENTS(info[2].As<v8::Float64Array>()));
	double *b = reinterpret_cast<double*>(GET_CONTENTS(info[3].As<v8::Float64Array>()));
	int i, j;
	for (i = 0 ; i < r ; i++) {
		for (j = 0 ; j < c ; j++) {
			b[j * r + i] = a[i * c + j];
		}
	}
}

void sTrTo(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int r = info[0]->Int32Value();
	int c = info[1]->Int32Value();
	float *a = reinterpret_cast<float*>(GET_CONTENTS(info[2].As<v8::Float32Array>()));
	float *b = reinterpret_cast<float*>(GET_CONTENTS(info[3].As<v8::Float32Array>()));
	int i, j;
	for (i = 0 ; i < r ; i++) {
		for (j = 0 ; j < c ; j++) {
			b[j * r + i] = a[i * c + j];
		}
	}
}

// transpose a -> a  ("ip" => "in place")
void dTrIp(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int r = info[0]->Int32Value();
	int c = info[1]->Int32Value();
	double *a = reinterpret_cast<double*>(GET_CONTENTS(info[2].As<v8::Float64Array>()));
	double* b = new double[r*c];
	memcpy(b, a, sizeof(double)*r*c);
	int i, j;
	for (i = 0 ; i < r ; i++) {
		for (j = 0 ; j < c ; j++) {
			a[j * r + i] = b[i * c + j];
		}
	}
	delete b;
}

void sTrIp(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int r = info[0]->Int32Value();
	int c = info[1]->Int32Value();
	float *a = reinterpret_cast<float*>(GET_CONTENTS(info[2].As<v8::Float32Array>()));
	float* b = new float[r*c];
	memcpy(b, a, sizeof(float)*r*c);
	int i, j;
	for (i = 0 ; i < r ; i++) {
		for (j = 0 ; j < c ; j++) {
			a[j * r + i] = b[i * c + j];
		}
	}
	delete b;
}


// memcpy
void BufCopy(const v8::FunctionCallbackInfo<v8::Value>& info) {
	unsigned char *dst = reinterpret_cast<unsigned char*>(GET_CONTENTS(info[0].As<v8::Int8Array>()));
	int dstOffsetBytes = info[1]->Int32Value();
	unsigned char *src = reinterpret_cast<unsigned char*>(GET_CONTENTS(info[2].As<v8::Int8Array>()));
	int srcOffsetBytes = info[3]->Int32Value();
	int bytes = info[4]->Int32Value();
	memcpy((void*)(dst + dstOffsetBytes), (void*)(src + srcOffsetBytes), bytes);
}

// memset
void BufSet(const v8::FunctionCallbackInfo<v8::Value>& info) {
	unsigned char *dst = reinterpret_cast<unsigned char*>(GET_CONTENTS(info[0].As<v8::Int8Array>()));
	int dstOffsetBytes = info[1]->Int32Value();
	int val = info[2]->Int32Value(); //int, but interpreted as an unsigned char
	int bytes = info[3]->Int32Value();
	memset((void*)(dst + dstOffsetBytes), val, bytes);
}


// diagonal matrix
void dMatrixDiagonal(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int r = info[0]->Int32Value();
	int c = info[1]->Int32Value();
	double *a = reinterpret_cast<double*>(GET_CONTENTS(info[2].As<v8::Float64Array>()));
	double val = info[3]->NumberValue();

	memset(a, 0, sizeof(double) * r * c);
	int i, m = (r < c ? r : c);
	for(i = 0 ; i < m ; i++) {
		a[i * c + i] = val;
	}
}
void sMatrixDiagonal(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int r = info[0]->Int32Value();
	int c = info[1]->Int32Value();
	float *a = reinterpret_cast<float*>(GET_CONTENTS(info[2].As<v8::Float32Array>()));
	float val = (float) info[3]->NumberValue();

	memset(a, 0, sizeof(float) * r * c);
	int i, m = (r < c ? r : c);
	for(i = 0 ; i < m ; i++) {
		a[i * c + i] = val;
	}
}

// matrix of ones or other number
void dMatrixFill(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int r = info[0]->Int32Value();
	int c = info[1]->Int32Value();
	double *a = reinterpret_cast<double*>(GET_CONTENTS(info[2].As<v8::Float64Array>()));
	double val = info[3]->NumberValue();

	int len = r * c, i;
	for(i = 0 ; i < len ; i++) {
		a[i] = val;
	}
}
void sMatrixFill(const v8::FunctionCallbackInfo<v8::Value>& info) {
	int r = info[0]->Int32Value();
	int c = info[1]->Int32Value();
	float *a = reinterpret_cast<float*>(GET_CONTENTS(info[2].As<v8::Float32Array>()));
	float val = (float) info[3]->NumberValue();

	int len = r * c, i;
	for(i = 0 ; i < len ; i++) {
		a[i] = val;
	}
}