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
