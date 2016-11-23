(function () {
  'use strict';

  var nblas = require('./build/Release/addon');
  var assert = require('assert');

  /**
   * @doc
   * BLAS:
   *  https://software.intel.com/ru-ru/node/468390
   *  https://software.intel.com/ru-ru/node/468426
   *  https://software.intel.com/ru-ru/node/468478
   * LAPACK:
   *  http://physics.oregonstate.edu/~landaur/nacphy/lapack/simple.html
   *  https://software.intel.com/ru-ru/node/468874
   * SPBLAS:
   *  http://math.nist.gov/spblas/
   *  http://www.cerfacs.fr/algor/reports/2001/TR_PA_01_24.pdf
   *  http://www.netlib.org/blas/blast-forum/chapter3.pdf
   *
   *  see also this alt libs
   *  http://librsb.sourceforge.net/
   *  https://github.com/xianyi/OpenBLAS
   **/

  // from enums declared in functions/cblas.h
  nblas.NoTrans = 111;
  nblas.Trans = 112;
  nblas.ConjTrans = 113;

  nblas.Upper = 121;
  nblas.Lower = 122;

  nblas.NonUnit = 131;
  nblas.Unit = 132;

  nblas.Left = 141;
  nblas.Right = 142;

  // from enums declared in lib/blas_enum.h
  nblas.SymmetryType = {
    blas_general          : 231,
    blas_symmetric        : 232,
    blas_hermitian        : 233,
    blas_triangular       : 234,
    blas_lower_triangular : 235,
    blas_upper_triangular : 236,
    blas_lower_symmetric  : 237,
    blas_upper_symmetric  : 238,
    blas_lower_hermitian  : 239,
    blas_upper_hermitian  : 240
  };
  nblas.FieldType = {
    blas_complex          : 241, //not supported here
    blas_real             : 242, //not supported here
    blas_double_precision : 243,
    blas_single_precision : 244
  };
  nblas.SizeType = {
    blas_num_rows      : 251,
    blas_num_cols      : 252,
    blas_num_nonzeros  : 253
  };
  nblas.HandleType = {
    blas_invalid_handle : 261,
    blas_new_handle     : 262,
    blas_open_handle    : 263,
    blas_valid_handle   : 264
  };

  // for lapack
  nblas.Lapack = {
    Fact: {
      N: 0,
      F: 1,
      E: 2,
    },
    Trans: {
      N: 0,
      T: 1,
      C: 2,
    },
    Equed: {
      N: 0,
      R: 1,
      C: 2,
      B: 3,
    },
  };


  // enforce strict type checking
  function typeCheck(array) {
    if (array.constructor === Float64Array)
      return true;
    else if (array.constructor === Float32Array)
      return false;

    throw new Error('invalid type!');
  }

  // BLAS Level 1 Routines and Functions
  // https://software.intel.com/ru-ru/node/468390
  nblas.asum = function (x) {
    return typeCheck(x) ?
      nblas.dasum(x.length, x, 1) :
      nblas.sasum(x.length, x, 1);
  };

  nblas.axpy = function (x, y, alpha) {
    alpha = alpha || 1.0;
    return typeCheck(x) ?
      nblas.daxpy(x.length, alpha, x, 1, y, 1) :
      nblas.saxpy(x.length, alpha, x, 1, y, 1);
  };

  nblas.copy = function (x, y) {
    return typeCheck(x) ?
      nblas.dcopy(x.length, x, 1, y, 1) :
      nblas.scopy(x.length, x, 1, y, 1);
  };

  nblas.dot = function (x, y) {
    return typeCheck(x) ?
      nblas.ddot(x.length, x, 1, y, 1) :
      nblas.sdot(x.length, x, 1, y, 1);
  };

  nblas.nrm2 = function (x) {
    return typeCheck(x) ?
      nblas.dnrm2(x.length, x, 1) :
      nblas.snrm2(x.length, x, 1);
  };

  nblas.rot = function (x, y, c, s) {
    return typeCheck(x) ?
      nblas.drot(x.length, x, 1, y, 1, c, s) :
      nblas.srot(x.length, x, 1, y, 1, c, s);
  };

  nblas.rotg = function (x, y, c, s) {
    return typeCheck(x) ?
      nblas.drotg(x, y, c, s) :
      nblas.srotg(x, y, c, s);
  };

  nblas.rotm = function (x, y, param) {
    return typeCheck(x) ?
      nblas.drotm(x.length, x, 1, y, 1, param) :
      nblas.srotm(x.length, x, 1, y, 1, param);
  };

  nblas.rotmg = function (d1, d2, x1, y1, param) {
    return typeCheck(x) ?
      nblas.drotmg(d1, d2, x1, y1, param) :
      nblas.srotmg(d1, d2, x1, y1, param);
  };

  nblas.scal = function (x, alpha) {
    return typeCheck(x) ?
      nblas.dscal(x.length, alpha, x, 1) :
      nblas.sscal(x.length, alpha, x, 1);
  };

  nblas.swap = function (x, y) {
    return typeCheck(x) ?
      nblas.dswap(x.length, x, 1, y, 1) :
      nblas.sswap(x.length, x, 1, y, 1);
  };

  nblas.iamax = function (x) {
    return typeCheck(x) ?
      nblas.idamax(x.length, x, 1) :
      nblas.isamax(x.length, x, 1);
  };

  nblas.iamin = function (x) {
    return typeCheck(x) ?
      nblas.idamin(x.length, x, 1) :
      nblas.isamin(x.length, x, 1);
  };

  // BLAS Level 2 Routines
  // https://software.intel.com/ru-ru/node/468426
  nblas.gbmv = function (a, x, y, kl, ku, alpha, beta, trans) {
    trans = trans || nblas.NoTrans;
    kl = kl || 0;
    ku = ku || 0;
    alpha = alpha || 1.0;
    beta = beta || 0.0;
    return typeCheck(a) ?
      nblas.dgbmv(trans, x.length, y.length, kl, ku, alpha, a, x.length, x, 1, beta, y, 1) :
      nblas.sgbmv(trans, x.length, y.length, kl, ku, alpha, a, x.length, x, 1, beta, y, 1);
  };

  nblas.gemv = function (a, x, y, alpha, beta, trans) {
    trans = trans || nblas.NoTrans;
    alpha = alpha || 1.0;
    beta = beta || 0.0;
    return typeCheck(a) ?
      nblas.dgemv(trans, x.length, y.length, alpha, a, x.length, x, 1, beta, y, 1) :
      nblas.sgemv(trans, x.length, y.length, alpha, a, x.length, x, 1, beta, y, 1);
  };

  nblas.ger = function (a, x, y, alpha) {
    alpha = alpha || 1.0;
    return typeCheck(a) ?
      nblas.dger(x.length, y.length, alpha, x, 1, y, 1, a, x.length) :
      nblas.sger(x.length, y.length, alpha, x, 1, y, 1, a, x.length);
  };

  nblas.sbmv = function (a, x, y, k, uplo, alpha, beta) {
    uplo = uplo || nblas.Upper;
    k = k || 0;
    alpha = alpha || 1.0;
    beta = beta || 0.0;
    return typeCheck(a) ?
      nblas.dsbmv(uplo, x.length, k, alpha, a, x.length, x, 1, beta, y, 1) :
      nblas.ssbmv(uplo, x.length, k, alpha, a, x.length, x, 1, beta, y, 1);
  };

  nblas.spmv = function (ap, x, y, uplo, alpha, beta) {
    uplo = uplo || nblas.Upper;
    alpha = alpha || 1.0;
    beta = beta || 0.0;
    return typeCheck(ap) ?
      nblas.dspmv(uplo, x.length, alpha, ap, x, 1, beta, y, 1) :
      nblas.sspmv(uplo, x.length, alpha, ap, x, 1, beta, y, 1);
  };

  nblas.spr = function (ap, x, uplo, alpha) {
    uplo = uplo || nblas.Upper;
    alpha = alpha || 1.0;
    return typeCheck(ap) ?
      nblas.dspr(uplo, x.length, alpha, x, 1, ap) :
      nblas.sspr(uplo, x.length, alpha, x, 1, ap);
  };

  nblas.spr2 = function (ap, x, y, uplo, alpha) {
    uplo = uplo || nblas.Upper;
    alpha = alpha || 1.0;
    return typeCheck(ap) ?
      nblas.dspr2(uplo, x.length, alpha, x, 1, y, 1, ap) :
      nblas.sspr2(uplo, x.length, alpha, x, 1, y, 1, ap);
  };

  nblas.symv = function (a, x, y, uplo, alpha, beta) {
    uplo = uplo || nblas.Upper;
    alpha = alpha || 1.0;
    beta = beta || 0.0;
    return typeCheck(a) ?
      nblas.dsymv(uplo, x.length, alpha, a, x.length, x, 1, beta, y, 1) :
      nblas.ssymv(uplo, x.length, alpha, a, x.length, x, 1, beta, y, 1);
  };

  nblas.syr = function (a, x, uplo, alpha) {
    uplo = uplo || nblas.Upper;
    alpha = alpha || 1.0;
    return typeCheck(a) ?
      nblas.dsyr(uplo, x.length, alpha, x, 1, a, x.length) :
      nblas.ssyr(uplo, x.length, alpha, x, 1, a, x.length);
  };

  nblas.syr2 = function (a, x, y, uplo, alpha) {
    uplo = uplo || nblas.Upper;
    alpha = alpha || 1.0;
    return typeCheck(a) ?
      nblas.dsyr2(uplo, x.length, alpha, x, 1, y, 1, a, x.length) :
      nblas.ssyr2(uplo, x.length, alpha, x, 1, y, 1, a, x.length);
  };

  nblas.tbmv = function (a, x, y, uplo, trans, diag) {
    uplo = uplo || nblas.Upper;
    trans = trans || nblas.NoTrans;
    diag = diag || nblas.NonUnit;
    return typeCheck(a) ?
      nblas.dtbmv(uplo, trans, diag, x.length, 0, a, x.length, x, 1) :
      nblas.stbmv(uplo, trans, diag, x.length, 0, a, x.length, x, 1);
  };

  nblas.tbsv = function (a, x, uplo, trans, diag) {
    uplo = uplo || nblas.Upper;
    trans = trans || nblas.NoTrans;
    diag = diag || nblas.NonUnit;
    return typeCheck(a) ?
      nblas.dtbsv(uplo, trans, diag, x.length, 0, a, x.length, x, 1) :
      nblas.stbsv(uplo, trans, diag, x.length, 0, a, x.length, x, 1);
  };

  nblas.tpmv = function (ap, x, uplo, trans, diag) {
    uplo = uplo || nblas.Upper;
    trans = trans || nblas.NoTrans;
    diag = diag || nblas.NonUnit;
    return typeCheck(ap) ?
      nblas.dtpmv(uplo, trans, diag, x.length, ap, x, 1) :
      nblas.stpmv(uplo, trans, diag, x.length, ap, x, 1);
  };

  nblas.tpsv = function (ap, x, uplo, trans, diag) {
    uplo = uplo || nblas.Upper;
    trans = trans || nblas.NoTrans;
    diag = diag || nblas.NonUnit;
    return typeCheck(ap) ?
      nblas.dtpsv(uplo, trans, diag, x.length, ap, x, 1) :
      nblas.stpsv(uplo, trans, diag, x.length, ap, x, 1);
  };

  nblas.trmv = function (a, x, uplo, trans, diag) {
    uplo = uplo || nblas.Upper;
    trans = trans || nblas.NoTrans;
    diag = diag || nblas.NonUnit;
    return typeCheck(a) ?
      nblas.dtrmv(uplo, trans, diag, x.length, a, x.length, x, 1) :
      nblas.strmv(uplo, trans, diag, x.length, a, x.length, x, 1);
  };

  nblas.trsv = function (a, x, uplo, trans, diag) {
    uplo = uplo || nblas.Upper;
    trans = trans || nblas.NoTrans;
    diag = diag || nblas.NonUnit;
    return typeCheck(a) ?
      nblas.dtrsv(uplo, trans, diag, x.length, a, x.length, x, 1) :
      nblas.strsv(uplo, trans, diag, x.length, a, x.length, x, 1);
  };

  // BLAS Level 3 Routines
  // https://software.intel.com/ru-ru/node/468478
  nblas.gemm = function (a, b, c, m, n, k, transa, transb, alpha, beta) {
    transa = transa || nblas.NoTrans;
    transb = transb || nblas.NoTrans;
    alpha = alpha || 1.0;
    beta = beta || 0.0;
    var lda = (transa == nblas.NoTrans ? k : m);
    var ldb = (transa == nblas.NoTrans ? n : k);
    var ldc = n;
    return typeCheck(a) ?
      nblas.dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) :
      nblas.sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
  };

  nblas.symm = function (a, b, c, m, n, side, uplo, alpha, beta) {
    side = side || nblas.Left;
    uplo = uplo || nblas.Upper;
    alpha = alpha || 1.0;
    beta = beta || 0.0;
    return typeCheck(a) ?
      nblas.dsymm(side, uplo, m, n, alpha, a, m, b, n, beta, c, m) :
      nblas.ssymm(side, uplo, m, n, alpha, a, m, b, n, beta, c, m);
  };

  nblas.syrk = function (a, c, n, k, uplo, trans, alpha, beta) {
    uplo = uplo || nblas.Upper;
    trans = trans || nblas.NoTrans;
    alpha = alpha || 1.0;
    beta = beta || 0.0;
    var lda = (trans == nblas.NoTrans ? n : k);
    var ldc = n;
    return typeCheck(a) ?
      nblas.dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc) :
      nblas.ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);
  };

  nblas.syr2k = function (a, b, c, n, k, uplo, trans, alpha, beta) {
    uplo = uplo || nblas.Upper;
    trans = trans || nblas.NoTrans;
    alpha = alpha || 1.0;
    beta = beta || 0.0;
    var lda = (trans == nblas.NoTrans ? n : k);
    var ldb = (trans == nblas.NoTrans ? n : k);
    var ldc = n;
    return typeCheck(a) ?
      nblas.dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc) :
      nblas.ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
  };

  nblas.trmm = function (a, b, m, n, side, uplo, transa, diag, alpha) {
    side = side || nblas.Left;
    uplo = uplo || nblas.Upper;
    transa = transa || nblas.NoTrans;
    diag = diag || nblas.NonUnit;
    alpha = alpha || 1.0;
    var lda = m;
    var ldb = m;
    return typeCheck(a) ?
      nblas.dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) :
      nblas.strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
  };

  nblas.trsm = function (a, b, m, n, side, uplo, transa, diag, alpha) {
    side = side || nblas.Left;
    uplo = uplo || nblas.Upper;
    transa = transa || nblas.NoTrans;
    diag = diag || nblas.NonUnit;
    alpha = alpha || 1.0;
    var lda = m;
    var ldb = m;
    return typeCheck(a) ?
      nblas.dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) :
      nblas.strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
  };


  // LAPACK
  // http://physics.oregonstate.edu/~landaur/nacphy/lapack/simple.html
  // https://software.intel.com/ru-ru/node/468874

  // A [m,m] * x [m,n] = B [m,n]
  nblas.gesv = function (A, B, m, n, ipiv) {
    var type = A.constructor;
    if (ipiv === undefined)
      ipiv = new type(m);
    assert(A.length >= m*m);
    assert(B.length >= m*n);
    assert(ipiv.length == m);

    return typeCheck(A) ?
      nblas.dgesv(m, n, A, B, ipiv) :
      nblas.sgesv(m, n, A, B, ipiv);
  };

  nblas.gesvx = function (A, B, X, m, n, AF, IPIV, fact, trans, equed, R, C) {
    var type = A.constructor;
    fact = fact || nblas.Lapack.Fact.N;
    trans = trans || nblas.Lapack.Trans.N;
    equed = equed || nblas.Lapack.Equed.N;
    if (R === undefined)
      R = new type(m);
    if (C === undefined)
      C = new type(m);
    if (IPIV === undefined && fact == nblas.Lapack.Fact.N)
      IPIV = new Int32Array(m);
    if (AF === undefined && fact == nblas.Lapack.Fact.N)
      AF = new type(m*m);

    assert(A.length >= m*m);
    assert(B.length >= m*n);
    assert(X.length >= m*n);
    assert(AF.length >= m*m);
    assert(IPIV.length >= 1*m);
    assert(R.length >= 1*m);
    assert(C.length >= 1*m);

    return typeCheck(A) ?
      nblas.dgesvx(fact, trans, m, n, A, AF, IPIV, equed, R, C, B, X) :
      nblas.sgesvx(fact, trans, m, n, A, AF, IPIV, equed, R, C, B, X);
  };

  nblas.getrf = function (A, IPIV, m, n) {
    assert(A.length >= m*n);
    assert(IPIV.length >= Math.min(m,n));

    return typeCheck(A) ?
      nblas.dgetrf(m, n, A, IPIV) :
      nblas.sgetrf(m, n, A, IPIV);
  }


  // SPBLAS
  // http://math.nist.gov/spblas/
  // http://www.cerfacs.fr/algor/reports/2001/TR_PA_01_24.pdf
  // http://www.netlib.org/blas/blast-forum/chapter3.pdf

  // SPBLAS Creation
  // Construction
  nblas.uscr_begin = function(double, m, n) {
    return double ?
      nblas.duscr_begin(m, n) :
      nblas.suscr_begin(m, n);
  };
  // Block construction
  // Mb, Nb - blocks count, k, l = blocks size
  // M = Mb * k, N = Nb * l
  nblas.uscr_block_begin = function(double, Mb, Nb, k, l) {
    return double ?
      nblas.duscr_block_begin(Mb, Nb, k, l) :
      nblas.suscr_block_begin(Mb, Nb, k, l);
  };
  // Variable block construction
  // K - array of size Mb, L - array of size Nb
  nblas.uscr_variable_block_begin = function(double, Mb, Nb, K, L) {
    return double ?
      nblas.duscr_variable_block_begin(Mb, Nb, K, L) :
      nblas.suscr_variable_block_begin(Mb, Nb, K, L);
  };
  // Insertion
  nblas.uscr_insert_entry = function(A, val, i, j) {
    var double = nblas.usgp(A, nblas.FieldType.blas_double_precision);
    return double ?
      nblas.duscr_insert_entry(A, val, i, j) :
      nblas.suscr_insert_entry(A, val, i, j);
  };
  nblas.uscr_insert_entries = function(A, nz, vals, indx, jndx) {
    return typeCheck(vals) ?
      nblas.duscr_insert_entries(A, nz, vals, indx, jndx) :
      nblas.suscr_insert_entries(A, nz, vals, indx, jndx);
  };
  nblas.uscr_insert_col = function(A, j, nz, vals, indx) {
    return typeCheck(vals) ?
      nblas.duscr_insert_col(A, j, nz, vals, indx) :
      nblas.suscr_insert_col(A, j, nz, vals, indx);
  };
  nblas.uscr_insert_row = function(A, i, nz, vals, jndx) {
    return typeCheck(vals) ?
      nblas.duscr_insert_row(A, i, nz, vals, jndx) :
      nblas.suscr_insert_row(A, i, nz, vals, jndx);
  };
  nblas.uscr_insert_clique = function(A, k, l, vals, row_stride, col_stride, indx, jndx) {
    return typeCheck(vals) ?
      nblas.duscr_insert_clique(A, k, l, vals, row_stride, col_stride, indx, jndx) :
      nblas.suscr_insert_clique(A, k, l, vals, row_stride, col_stride, indx, jndx);
  };
  nblas.uscr_insert_block = function(A, vals, row_stride, col_stride, i, j) {
    return typeCheck(vals) ?
      nblas.duscr_insert_block(A, vals, row_stride, col_stride, i, j) :
      nblas.suscr_insert_block(A, vals, row_stride, col_stride, i, j);
  };
  // Completion of Construction Routines
  nblas.uscr_end = function(A) {
    var double = nblas.usgp(A, nblas.FieldType.blas_double_precision);
    return double ?
      nblas.duscr_end(A) :
      nblas.suscr_end(A);
  };
  // Matrix Property Routines
  nblas.usgp = function(A, pname) {
    return nblas._usgp(A, pname);
  };
  nblas.ussp = function(A, pname) {
    return nblas._ussp(A, pname);
  };
  // Destruction Routine
  nblas.usds = function(A) {
    return nblas._usds(A);
  };


  // SPBLAS Level 1
  // sparse dot product 
  nblas.usdot = function(x, indx, y) {
    var incy = 1;
    return typeCheck(x) ?
      nblas.dusdot(x.length, x, indx, y, incy) :
      nblas.susdot(x.length, x, indx, y, incy);
  };
  // sparse vector update
  nblas.usaxpy = function(x, indx, y, alpha) {
    var incy = 1;
    return typeCheck(x) ?
      nblas.dusaxpy(x.length, alpha, x, indx, y, incy) :
      nblas.susaxpy(x.length, alpha, x, indx, y, incy);
  };
  // sparse gather
  nblas.usga = function(x, indx, y) {
    var incy = 1;
    return typeCheck(x) ?
      nblas.dusga(x.length, y, incy, x, indx) :
      nblas.susga(x.length, y, incy, x, indx);
  };
  // sparse gather and zero
  nblas.usgz = function(x, indx, y) {
    var incy = 1;
    return typeCheck(x) ?
      nblas.dusgz(x.length, y, incy, x, indx) :
      nblas.susgz(x.length, y, incy, x, indx);
  };
  // sparse scatter
  nblas.ussc = function(x, indx, y) {
    var incy = 1;
    return typeCheck(x) ?
      nblas.dussc(x.length, x, y, incy, indx) :
      nblas.sussc(x.length, x, y, incy, indx);
  };


  // SPBLAS Level 2
  // sparse matrix-vector multiply
  // A [m * n] * x [ n * 1 ] = y [ m * 1 ]
  nblas.usmv = function (A, x, y, trans, alpha) {
    trans = trans || nblas.NoTrans;
    alpha = alpha || 1.0;
    var incx = 1;
    var incy = 1;
    var m = nblas.usgp(A, nblas.SizeType.blas_num_rows);
    var n = nblas.usgp(A, nblas.SizeType.blas_num_cols);
    assert(x.length >= n);
    assert(y.length >= m);
    return typeCheck(x) ?
      nblas.dusmv(trans, alpha, A, x, incx, y, incy) :
      nblas.susmv(trans, alpha, A, x, incx, y, incy);
  }
  // sparse triangular solve 
  nblas.ussv = function (A, x, trans, alpha) {
    trans = trans || nblas.NoTrans;
    alpha = alpha || 1.0;
    var incx = 1;
    var m = nblas.usgp(A, nblas.SizeType.blas_num_rows);
    var n = nblas.usgp(A, nblas.SizeType.blas_num_cols);
    assert(x.length >= n);
    return typeCheck(x) ?
      nblas.dussv(trans, alpha, A, x, incx) :
      nblas.sussv(trans, alpha, A, x, incx);
  }

  // SPBLAS Level 3
  // sparse matrix-matrix multiply
  nblas.usmm = function(A, B, C, nrhs, trans, alpha) {
    trans = trans || nblas.NoTrans;
    alpha = alpha || 1.0;
    var m = nblas.usgp(A, nblas.SizeType.blas_num_rows);
    var n = nblas.usgp(A, nblas.SizeType.blas_num_cols);
    assert(B.length >= n * nrhs);
    assert(C.length >= m * nrhs);
    return typeCheck(B) ?
      nblas.dusmm(trans, nrhs, alpha, A, B, C) :
      nblas.susmm(trans, nrhs, alpha, A, B, C);
  };
  // sparse triangular solve
  nblas.ussm = function(A, B, nrhs, trans, alpha) {
    trans = trans || nblas.NoTrans;
    alpha = alpha || 1.0;
    var m = nblas.usgp(A, nblas.SizeType.blas_num_rows);
    var n = nblas.usgp(A, nblas.SizeType.blas_num_cols);
    assert(B.length >= n * nrhs);
    return typeCheck(B) ?
      nblas.dussm(trans, nrhs, alpha, A, B) :
      nblas.sussm(trans, nrhs, alpha, A, B);
  };


  // Other
  nblas.TrTo = function(a, b, m, n) {
    return typeCheck(a) ?
      nblas.dTrTo(m, n, a, b) :
      nblas.sTrTo(m, n, a, b);
  };
  nblas.TrIp = function(a, m, n) {
    return typeCheck(a) ?
      nblas.dTrIp(m, n, a) :
      nblas.sTrIp(m, n, a);
  };
  //nblas.BufCopy(dst, src, bytes)


  module.exports = nblas;
}());
