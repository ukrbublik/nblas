#include "cblas.h"
#include "routines.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <lapacke.h>

//http://physics.oregonstate.edu/~landaur/nacphy/lapack/routines/dgesvx.html
void dgesvx(const v8::FunctionCallbackInfo<v8::Value>& info) {
	//sizes:
	//A - (LDA,N)
	//X - (LDX,NRHS)
	//B - (LDB,NRHS)
	//af - (LDAF,N)
	//ipiv - (N)
	//r,c - (N)
	//ferr,berr - (NRHS)
	//LDA,LDX,LDB,LDAF >= N
	const enum LAPACK_FACT ifact = static_cast<LAPACK_FACT>(info[0]->Uint32Value());
	char fact = (ifact == LAPACK_FACT_F ? 'F' : (ifact == LAPACK_FACT_E ? 'E' : 'N'));
	const enum LAPACK_TRANS itrans = static_cast<LAPACK_TRANS>(info[1]->Uint32Value());
	char trans = (itrans == LAPACK_TRANS_T ? 'T' : (itrans == LAPACK_TRANS_C ? 'C' : 'N'));
	int n = info[2]->Int32Value();
	int nrhs = info[3]->Int32Value();
	double *a = reinterpret_cast<double*>(GET_CONTENTS(info[4].As<v8::Float64Array>()));
	int lda = n;
	double *af = reinterpret_cast<double*>(GET_CONTENTS(info[5].As<v8::Float64Array>()));
	int ldaf = n;
	int *ipiv = reinterpret_cast<int*>(GET_CONTENTS(info[6].As<v8::Int32Array>()));
	const enum LAPACK_EQUED iequed = static_cast<LAPACK_EQUED>(info[7]->Uint32Value());
	char equed = (iequed == LAPACK_EQUED_R ? 'R' : (iequed == LAPACK_EQUED_C ? 'C' 
		: (iequed == LAPACK_EQUED_B ? 'B' : 'N')));
	double *r = reinterpret_cast<double*>(GET_CONTENTS(info[8].As<v8::Float64Array>()));
	double *c = reinterpret_cast<double*>(GET_CONTENTS(info[9].As<v8::Float64Array>()));
	double *b = reinterpret_cast<double*>(GET_CONTENTS(info[10].As<v8::Float64Array>()));
	int ldb = n;
	double *x = reinterpret_cast<double*>(GET_CONTENTS(info[11].As<v8::Float64Array>()));
	int ldx = n;
	double rcond;
	double ferr[nrhs];
	double berr[nrhs];
	double rpivot;

	int res;
	int matrix_layout = LAPACK_ROW_MAJOR;
	if (matrix_layout == LAPACK_ROW_MAJOR) {
		ldb = nrhs;
		ldx = nrhs;
	}
	res = LAPACKE_dgesvx( matrix_layout, fact, trans, n, nrhs, a, lda, af, ldaf, ipiv,
		 &equed, r, c, b, ldb, x, ldx, &rcond, ferr, berr, &rpivot );
	info.GetReturnValue().Set(res);
}

void sgesvx(const v8::FunctionCallbackInfo<v8::Value>& info) {
	//sizes:
	//A - (LDA,N)
	//X - (LDX,NRHS)
	//B - (LDB,NRHS)
	//af - (LDAF,N)
	//ipiv - (N)
	//r,c - (N)
	//ferr,berr - (NRHS)
	//LDA,LDX,LDB,LDAF >= N
	const enum LAPACK_FACT ifact = static_cast<LAPACK_FACT>(info[0]->Uint32Value());
	char fact = (ifact == LAPACK_FACT_F ? 'F' : (ifact == LAPACK_FACT_E ? 'E' : 'N'));
	const enum LAPACK_TRANS itrans = static_cast<LAPACK_TRANS>(info[1]->Uint32Value());
	char trans = (itrans == LAPACK_TRANS_T ? 'T' : (itrans == LAPACK_TRANS_C ? 'C' : 'N'));
	int n = info[2]->Int32Value();
	int nrhs = info[3]->Int32Value();
	double *a = reinterpret_cast<double*>(GET_CONTENTS(info[4].As<v8::Float64Array>()));
	int lda = n;
	double *af = reinterpret_cast<double*>(GET_CONTENTS(info[5].As<v8::Float64Array>()));
	int ldaf = n;
	int *ipiv = reinterpret_cast<int*>(GET_CONTENTS(info[6].As<v8::Int32Array>()));
	const enum LAPACK_EQUED iequed = static_cast<LAPACK_EQUED>(info[7]->Uint32Value());
	char equed = (iequed == LAPACK_EQUED_R ? 'R' : (iequed == LAPACK_EQUED_C ? 'C' 
		: (iequed == LAPACK_EQUED_B ? 'B' : 'N')));
	double *r = reinterpret_cast<double*>(GET_CONTENTS(info[8].As<v8::Float64Array>()));
	double *c = reinterpret_cast<double*>(GET_CONTENTS(info[9].As<v8::Float64Array>()));
	double *b = reinterpret_cast<double*>(GET_CONTENTS(info[10].As<v8::Float64Array>()));
	int ldb = n;
	double *x = reinterpret_cast<double*>(GET_CONTENTS(info[11].As<v8::Float64Array>()));
	int ldx = n;
	double rcond;
	double ferr[nrhs];
	double berr[nrhs];
	double rpivot;

	int res;
	int matrix_layout = LAPACK_ROW_MAJOR;
	if (matrix_layout == LAPACK_ROW_MAJOR) {
		ldb = nrhs;
		ldx = nrhs;
	}
	res = LAPACKE_dgesvx( matrix_layout, fact, trans, n, nrhs, a, lda, af, ldaf, ipiv,
		 &equed, r, c, b, ldb, x, ldx, &rcond, ferr, berr, &rpivot );
	info.GetReturnValue().Set(res);
}
