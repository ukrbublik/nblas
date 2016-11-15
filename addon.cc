#include <node.h>
#include "routines.h"

void Init(v8::Local<v8::Object> exports) {
  // BLAS Level 1
  NODE_SET_METHOD(exports, "dasum", dasum);
  NODE_SET_METHOD(exports, "sasum", sasum);
  NODE_SET_METHOD(exports, "daxpy", daxpy);
  NODE_SET_METHOD(exports, "saxpy", saxpy);
  NODE_SET_METHOD(exports, "dcopy", dcopy);
  NODE_SET_METHOD(exports, "scopy", scopy);
  NODE_SET_METHOD(exports, "ddot", ddot);
  NODE_SET_METHOD(exports, "sdot", sdot);
  NODE_SET_METHOD(exports, "dsdot", dsdot);
  NODE_SET_METHOD(exports, "sdsdot", sdsdot);
  NODE_SET_METHOD(exports, "dnrm2", dnrm2);
  NODE_SET_METHOD(exports, "snrm2", snrm2);
  NODE_SET_METHOD(exports, "drot", drot);
  NODE_SET_METHOD(exports, "srot", srot);
  NODE_SET_METHOD(exports, "drotg", drotg);
  NODE_SET_METHOD(exports, "srotg", srotg);
  NODE_SET_METHOD(exports, "drotm", drotm);
  NODE_SET_METHOD(exports, "srotm", srotm);
  NODE_SET_METHOD(exports, "drotmg", drotmg);
  NODE_SET_METHOD(exports, "srotmg", srotmg);
  NODE_SET_METHOD(exports, "dscal", dscal);
  NODE_SET_METHOD(exports, "sscal", sscal);
  NODE_SET_METHOD(exports, "dswap", dswap);
  NODE_SET_METHOD(exports, "sswap", sswap);
  NODE_SET_METHOD(exports, "idamax", idamax);
  NODE_SET_METHOD(exports, "isamax", isamax);

  // BLAS Level 2
  NODE_SET_METHOD(exports, "dgbmv", dgbmv);
  NODE_SET_METHOD(exports, "sgbmv", sgbmv);
  NODE_SET_METHOD(exports, "dgemv", dgemv);
  NODE_SET_METHOD(exports, "sgemv", sgemv);
  NODE_SET_METHOD(exports, "dger", dger);
  NODE_SET_METHOD(exports, "sger", sger);
  NODE_SET_METHOD(exports, "dsbmv", dsbmv);
  NODE_SET_METHOD(exports, "ssbmv", ssbmv);
  NODE_SET_METHOD(exports, "dspmv", dspmv);
  NODE_SET_METHOD(exports, "sspmv", sspmv);
  NODE_SET_METHOD(exports, "dspr", dspr);
  NODE_SET_METHOD(exports, "sspr", sspr);
  NODE_SET_METHOD(exports, "dspr2", dspr2);
  NODE_SET_METHOD(exports, "sspr2", sspr2);
  NODE_SET_METHOD(exports, "dsymv", dsymv);
  NODE_SET_METHOD(exports, "ssymv", ssymv);
  NODE_SET_METHOD(exports, "dsyr", dsyr);
  NODE_SET_METHOD(exports, "ssyr", ssyr);
  NODE_SET_METHOD(exports, "dsyr2", dsyr2);
  NODE_SET_METHOD(exports, "ssyr2", ssyr2);
  NODE_SET_METHOD(exports, "dtbmv", dtbmv);
  NODE_SET_METHOD(exports, "stbmv", stbmv);
  NODE_SET_METHOD(exports, "dtbsv", dtbsv);
  NODE_SET_METHOD(exports, "stbsv", stbsv);
  NODE_SET_METHOD(exports, "dtpmv", dtpmv);
  NODE_SET_METHOD(exports, "stpmv", stpmv);
  NODE_SET_METHOD(exports, "dtpsv", dtpsv);
  NODE_SET_METHOD(exports, "stpsv", stpsv);
  NODE_SET_METHOD(exports, "dtrmv", dtrmv);
  NODE_SET_METHOD(exports, "strmv", strmv);
  NODE_SET_METHOD(exports, "dtrsv", dtrsv);
  NODE_SET_METHOD(exports, "strsv", strsv);

  // BLAS Level 3
  NODE_SET_METHOD(exports, "dgemm", dgemm);
  NODE_SET_METHOD(exports, "sgemm", sgemm);
  NODE_SET_METHOD(exports, "dsymm", dsymm);
  NODE_SET_METHOD(exports, "ssymm", ssymm);
  NODE_SET_METHOD(exports, "dsyrk", dsyrk);
  NODE_SET_METHOD(exports, "ssyrk", ssyrk);
  NODE_SET_METHOD(exports, "dsyr2k", dsyr2k);
  NODE_SET_METHOD(exports, "ssyr2k", ssyr2k);
  NODE_SET_METHOD(exports, "dtrmm", dtrmm);
  NODE_SET_METHOD(exports, "strmm", strmm);
  NODE_SET_METHOD(exports, "dtrsm", dtrsm);
  NODE_SET_METHOD(exports, "strsm", strsm);

  // LAPACK
  NODE_SET_METHOD(exports, "dgesv", dgesv);
  NODE_SET_METHOD(exports, "sgesv", sgesv);

  // SPBLAS Creation
  NODE_SET_METHOD(exports, "duscr_begin", duscr_begin);
  NODE_SET_METHOD(exports, "suscr_begin", suscr_begin);
  NODE_SET_METHOD(exports, "duscr_block_begin", duscr_block_begin);
  NODE_SET_METHOD(exports, "suscr_block_begin", suscr_block_begin);
  NODE_SET_METHOD(exports, "duscr_variable_block_begin", duscr_variable_block_begin);
  NODE_SET_METHOD(exports, "suscr_variable_block_begin", suscr_variable_block_begin);

  NODE_SET_METHOD(exports, "duscr_insert_entry", duscr_insert_entry);
  NODE_SET_METHOD(exports, "suscr_insert_entry", suscr_insert_entry);
  NODE_SET_METHOD(exports, "duscr_insert_entries", duscr_insert_entries);
  NODE_SET_METHOD(exports, "suscr_insert_entries", suscr_insert_entries);
  NODE_SET_METHOD(exports, "duscr_insert_col", duscr_insert_col);
  NODE_SET_METHOD(exports, "suscr_insert_col", suscr_insert_col);
  NODE_SET_METHOD(exports, "duscr_insert_row", duscr_insert_row);
  NODE_SET_METHOD(exports, "suscr_insert_row", suscr_insert_row);
  NODE_SET_METHOD(exports, "duscr_insert_clique", duscr_insert_clique);
  NODE_SET_METHOD(exports, "suscr_insert_clique", suscr_insert_clique);
  NODE_SET_METHOD(exports, "duscr_insert_block", duscr_insert_block);
  NODE_SET_METHOD(exports, "suscr_insert_block", suscr_insert_block);

  NODE_SET_METHOD(exports, "duscr_end", duscr_end);
  NODE_SET_METHOD(exports, "suscr_end", suscr_end);
  NODE_SET_METHOD(exports, "_usgp", _usgp);
  NODE_SET_METHOD(exports, "_ussp", _ussp);
  NODE_SET_METHOD(exports, "_usds", _usds);

  // SPBLAS Level 1
  NODE_SET_METHOD(exports, "dusdot", dusdot);
  NODE_SET_METHOD(exports, "susdot", susdot);
  NODE_SET_METHOD(exports, "dusaxpy", dusaxpy);
  NODE_SET_METHOD(exports, "susaxpy", susaxpy);
  NODE_SET_METHOD(exports, "dusga", dusga);
  NODE_SET_METHOD(exports, "susga", susga);
  NODE_SET_METHOD(exports, "dusgz", dusgz);
  NODE_SET_METHOD(exports, "susgz", susgz);
  NODE_SET_METHOD(exports, "dussc", dussc);
  NODE_SET_METHOD(exports, "sussc", sussc);

  // SPBLAS Level 2
  NODE_SET_METHOD(exports, "dusmv", dusmv);
  NODE_SET_METHOD(exports, "susmv", susmv);
  NODE_SET_METHOD(exports, "dussv", dussv);
  NODE_SET_METHOD(exports, "sussv", sussv);
  
  // SPBLAS Level 3
  NODE_SET_METHOD(exports, "dussm", dussm);
  NODE_SET_METHOD(exports, "sussm", sussm);
  NODE_SET_METHOD(exports, "dusmm", dusmm);
  NODE_SET_METHOD(exports, "susmm", susmm);

  // Other
  NODE_SET_METHOD(exports, "dTrTo", dTrTo);
  NODE_SET_METHOD(exports, "cTrTo", cTrTo);
  NODE_SET_METHOD(exports, "dTrIp", dTrIp);
  NODE_SET_METHOD(exports, "cTrIp", cTrIp);
}

NODE_MODULE(addon, Init)
