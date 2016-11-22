#ifndef ROUTINES_H
#define ROUTINES_H

#include <node.h>

#define GET_CONTENTS(view) \
(static_cast<unsigned char*>(view->Buffer()->GetContents().Data()) + view->ByteOffset())

#define MIN(a, b) (a < b ? a : b)
#define MAX(a, b) (a > b ? a : b)

// For LAPACK Expert Driver Routines
enum LAPACK_FACT {
	LAPACK_FACT_N = 0, // The matrix A will	be copied to AF and factored
	LAPACK_FACT_F = 1, // Use already factored form of A from AF and IPIV as input
	LAPACK_FACT_E = 2 // The matrix A will be equilibrated if necessary, then copied to AF	and factored
};
enum LAPACK_TRANS {
	LAPACK_TRANS_N = 0, // No transpose
	LAPACK_TRANS_T = 1, // Transpose A**T
	LAPACK_TRANS_C = 2 // Transpose A**H
};
enum LAPACK_EQUED {
	LAPACK_EQUED_N = 0, // No equilibration	(always	true if	FACT = 'N')
	LAPACK_EQUED_R = 1, // Row equilibration, i.e., A has been premultiplied by diag(R)
	LAPACK_EQUED_C = 2, // Column equilibration, i.e.,	A has been postmultiplied	by diag(C)
	LAPACK_EQUED_B = 3 // Both row and column equilibration,	i.e., A has	been replaced by diag(R) * A * diag(C)
};

// BLAS Level 1
void dasum(const v8::FunctionCallbackInfo<v8::Value>& info);
void sasum(const v8::FunctionCallbackInfo<v8::Value>& info);
void daxpy(const v8::FunctionCallbackInfo<v8::Value>& info);
void saxpy(const v8::FunctionCallbackInfo<v8::Value>& info);
void dcopy(const v8::FunctionCallbackInfo<v8::Value>& info);
void scopy(const v8::FunctionCallbackInfo<v8::Value>& info);
void ddot(const v8::FunctionCallbackInfo<v8::Value>& info);
void sdot(const v8::FunctionCallbackInfo<v8::Value>& info);
void dsdot(const v8::FunctionCallbackInfo<v8::Value>& info);
void sdsdot(const v8::FunctionCallbackInfo<v8::Value>& info);
void dnrm2(const v8::FunctionCallbackInfo<v8::Value>& info);
void snrm2(const v8::FunctionCallbackInfo<v8::Value>& info);
void drot(const v8::FunctionCallbackInfo<v8::Value>& info);
void srot(const v8::FunctionCallbackInfo<v8::Value>& info);
void drotg(const v8::FunctionCallbackInfo<v8::Value>& info);
void srotg(const v8::FunctionCallbackInfo<v8::Value>& info);
void drotm(const v8::FunctionCallbackInfo<v8::Value>& info);
void srotm(const v8::FunctionCallbackInfo<v8::Value>& info);
void drotmg(const v8::FunctionCallbackInfo<v8::Value>& info);
void srotmg(const v8::FunctionCallbackInfo<v8::Value>& info);
void dscal(const v8::FunctionCallbackInfo<v8::Value>& info);
void sscal(const v8::FunctionCallbackInfo<v8::Value>& info);
void dswap(const v8::FunctionCallbackInfo<v8::Value>& info);
void sswap(const v8::FunctionCallbackInfo<v8::Value>& info);
void idamax(const v8::FunctionCallbackInfo<v8::Value>& info);
void isamax(const v8::FunctionCallbackInfo<v8::Value>& info);

// BLAS Level 2
void dgbmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void sgbmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dgemv(const v8::FunctionCallbackInfo<v8::Value>& info);
void sgemv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dger(const v8::FunctionCallbackInfo<v8::Value>& info);
void sger(const v8::FunctionCallbackInfo<v8::Value>& info);
void dsbmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void ssbmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dspmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void sspmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dspr(const v8::FunctionCallbackInfo<v8::Value>& info);
void sspr(const v8::FunctionCallbackInfo<v8::Value>& info);
void dspr2(const v8::FunctionCallbackInfo<v8::Value>& info);
void sspr2(const v8::FunctionCallbackInfo<v8::Value>& info);
void dsymv(const v8::FunctionCallbackInfo<v8::Value>& info);
void ssymv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dsyr(const v8::FunctionCallbackInfo<v8::Value>& info);
void ssyr(const v8::FunctionCallbackInfo<v8::Value>& info);
void dsyr2(const v8::FunctionCallbackInfo<v8::Value>& info);
void ssyr2(const v8::FunctionCallbackInfo<v8::Value>& info);
void dtbmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void stbmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dtbsv(const v8::FunctionCallbackInfo<v8::Value>& info);
void stbsv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dtpmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void stpmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dtpsv(const v8::FunctionCallbackInfo<v8::Value>& info);
void stpsv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dtrmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void strmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dtrsv(const v8::FunctionCallbackInfo<v8::Value>& info);
void strsv(const v8::FunctionCallbackInfo<v8::Value>& info);

// BLAS Level 3
void dgemm(const v8::FunctionCallbackInfo<v8::Value>& info);
void sgemm(const v8::FunctionCallbackInfo<v8::Value>& info);
void dsymm(const v8::FunctionCallbackInfo<v8::Value>& info);
void ssymm(const v8::FunctionCallbackInfo<v8::Value>& info);
void dsyrk(const v8::FunctionCallbackInfo<v8::Value>& info);
void ssyrk(const v8::FunctionCallbackInfo<v8::Value>& info);
void dsyr2k(const v8::FunctionCallbackInfo<v8::Value>& info);
void ssyr2k(const v8::FunctionCallbackInfo<v8::Value>& info);
void dtrmm(const v8::FunctionCallbackInfo<v8::Value>& info);
void strmm(const v8::FunctionCallbackInfo<v8::Value>& info);
void dtrsm(const v8::FunctionCallbackInfo<v8::Value>& info);
void strsm(const v8::FunctionCallbackInfo<v8::Value>& info);

// LAPACK
void dgesv(const v8::FunctionCallbackInfo<v8::Value>& info);
void sgesv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dgesvx(const v8::FunctionCallbackInfo<v8::Value>& info);
void sgesvx(const v8::FunctionCallbackInfo<v8::Value>& info);
void dgetrf(const v8::FunctionCallbackInfo<v8::Value>& info);
void sgetrf(const v8::FunctionCallbackInfo<v8::Value>& info);

// SPBLAS Creation
void duscr_begin(const v8::FunctionCallbackInfo<v8::Value>& info);
void suscr_begin(const v8::FunctionCallbackInfo<v8::Value>& info);
void duscr_block_begin(const v8::FunctionCallbackInfo<v8::Value>& info);
void suscr_block_begin(const v8::FunctionCallbackInfo<v8::Value>& info);
void duscr_variable_block_begin(const v8::FunctionCallbackInfo<v8::Value>& info);
void suscr_variable_block_begin(const v8::FunctionCallbackInfo<v8::Value>& info);

void duscr_insert_entry(const v8::FunctionCallbackInfo<v8::Value>& info);
void suscr_insert_entry(const v8::FunctionCallbackInfo<v8::Value>& info);
void duscr_insert_entries(const v8::FunctionCallbackInfo<v8::Value>& info);
void suscr_insert_entries(const v8::FunctionCallbackInfo<v8::Value>& info);
void duscr_insert_col(const v8::FunctionCallbackInfo<v8::Value>& info);
void suscr_insert_col(const v8::FunctionCallbackInfo<v8::Value>& info);
void duscr_insert_row(const v8::FunctionCallbackInfo<v8::Value>& info);
void suscr_insert_row(const v8::FunctionCallbackInfo<v8::Value>& info);
void duscr_insert_clique(const v8::FunctionCallbackInfo<v8::Value>& info);
void suscr_insert_clique(const v8::FunctionCallbackInfo<v8::Value>& info);
void duscr_insert_block(const v8::FunctionCallbackInfo<v8::Value>& info);
void suscr_insert_block(const v8::FunctionCallbackInfo<v8::Value>& info);

void duscr_end(const v8::FunctionCallbackInfo<v8::Value>& info);
void suscr_end(const v8::FunctionCallbackInfo<v8::Value>& info);
void _usgp(const v8::FunctionCallbackInfo<v8::Value>& info);
void _ussp(const v8::FunctionCallbackInfo<v8::Value>& info);
void _usds(const v8::FunctionCallbackInfo<v8::Value>& info);

// SPBLAS Level 1
void dusdot(const v8::FunctionCallbackInfo<v8::Value>& info);
void susdot(const v8::FunctionCallbackInfo<v8::Value>& info);
void dusaxpy(const v8::FunctionCallbackInfo<v8::Value>& info);
void susaxpy(const v8::FunctionCallbackInfo<v8::Value>& info);
void dusga(const v8::FunctionCallbackInfo<v8::Value>& info);
void susga(const v8::FunctionCallbackInfo<v8::Value>& info);
void dusgz(const v8::FunctionCallbackInfo<v8::Value>& info);
void susgz(const v8::FunctionCallbackInfo<v8::Value>& info);
void dussc(const v8::FunctionCallbackInfo<v8::Value>& info);
void sussc(const v8::FunctionCallbackInfo<v8::Value>& info);

// SPBLAS Level 2
void dusmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void susmv(const v8::FunctionCallbackInfo<v8::Value>& info);
void dussv(const v8::FunctionCallbackInfo<v8::Value>& info);
void sussv(const v8::FunctionCallbackInfo<v8::Value>& info);

// SPBLAS Level 3
void dussm(const v8::FunctionCallbackInfo<v8::Value>& info);
void sussm(const v8::FunctionCallbackInfo<v8::Value>& info);
void dusmm(const v8::FunctionCallbackInfo<v8::Value>& info);
void susmm(const v8::FunctionCallbackInfo<v8::Value>& info);

// Other
void dTrTo(const v8::FunctionCallbackInfo<v8::Value>& info);
void sTrTo(const v8::FunctionCallbackInfo<v8::Value>& info);
void dTrIp(const v8::FunctionCallbackInfo<v8::Value>& info);
void sTrIp(const v8::FunctionCallbackInfo<v8::Value>& info);
void BufCopy(const v8::FunctionCallbackInfo<v8::Value>& info);

#endif
