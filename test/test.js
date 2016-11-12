(function () {
  'use strict';

  var assert = require('assert'),
      nblas = require('../addon');

  describe('?asum', function () {
    it('works for different sizes', function () {
      assert.equal(6, nblas.asum(new Float64Array([1, 2, 3])));
      assert.equal(6, nblas.asum(new Float64Array([1, 2, 0, 3])));
    });

    it('works for vectors containing negative values', function () {
      assert.equal(6, nblas.asum(new Float64Array([-1, -2, 1, 2])));
      assert.equal(13, nblas.asum(new Float64Array([-1, -2, 1, 2, 3, 4])));
    });
  });

  describe('?axpy', function () {
    it('works for a + b', function () {
      var a = new Float64Array([1, 2, 3]),
          b = new Float64Array([6, 5, 4]),
          ans = new Float64Array([7, 7, 7]);

      nblas.axpy(a, b);
      assert.deepEqual(ans, b);
    });

    it('works for (-a) + b', function () {
      var a = new Float64Array([1, 2, 3]),
          b = new Float64Array([6, 5, 4]),
          ans = new Float64Array([5, 3, 1]);

      nblas.axpy(a, b, -1);
      assert.deepEqual(ans, b);
    });

    it('works for (2a) + b', function () {
      var a = new Float64Array([-1, 3, 8, 1]),
          b = new Float64Array([-10, 2, 4, -1]),
          ans = new Float64Array([-12, 8, 20, 1]);

      nblas.axpy(a, b, 2);
      assert.deepEqual(ans, b);
    });
  });

  describe('?copy', function () {
    it('works as expected for different sizes', function () {
      var a = new Float64Array([1, 2, 3]),
          b = new Float64Array(3);

      nblas.copy(a, b);
      assert.deepEqual(a, b);

      a = new Float64Array([5, 1, 3, 8]);
      b = new Float64Array(4);
      nblas.copy(a, b);
      assert.deepEqual(a, b);
    });
  });

  describe('?dot', function () {
    it('works as expected for different sizes', function () {
      var a = new Float64Array([1, 2, 3]),
          b = new Float64Array([4, 5, 6]);

      assert.equal(32, nblas.dot(a, b));

      a = new Float64Array([-1, 3, 7, 4]);
      b = new Float64Array([2, 1, 3, 0]);

      assert.equal(22, nblas.dot(a, b));
    });
  });

  describe('?nrm2', function () {
    it('works as expected for different sizes', function () {
      var a = new Float64Array([1, 2, 3]);

      assert.equal(Math.sqrt(14).toPrecision(6), nblas.nrm2(a).toPrecision(6));

      a = new Float64Array([3, 7, 1, 0]);

      assert.equal(Math.sqrt(59).toPrecision(6), nblas.nrm2(a).toPrecision(6));
    });
  });

  describe('?rot', function () {
    // should perform plane rotation of points
  });

  describe('?scal', function () {
    it('works as expected for different sizes and negative values', function () {
      var a = new Float64Array([1, 2, 3]),
          ans = new Float64Array([2, 4, 6]);

      nblas.scal(a, 2);
      assert.deepEqual(ans, a);

      a = new Float64Array([1, -2, 3, 0]);
      ans = new Float64Array([-2, 4, -6, 0]);

      nblas.scal(a, -2);
      assert.deepEqual(ans, a);
    });
  });

  describe('?swap', function () {
    it('works as expected for different sizes', function () {
      var a = new Float64Array([1, 2, 3]),
          acopy = new Float64Array(a),
          b = new Float64Array(3),
          bcopy = new Float64Array(b);

      nblas.swap(a, b);
      assert.deepEqual(bcopy, a);
      assert.deepEqual(acopy, b);

      a = new Float64Array([3, 4, 7, -1]);
      b = new Float64Array([1, 3, -1, 0]);
      acopy = new Float64Array(a);
      bcopy = new Float64Array(b);

      nblas.swap(a, b);
      assert.deepEqual(bcopy, a);
      assert.deepEqual(acopy, b);
    });
  });

  describe('i?amax', function () {
    it('works as expected for different values', function () {
      var a = new Float64Array([1, 2, 3]);

      assert.equal(3, a[nblas.iamax(a)]);

      a = new Float64Array([-1, -2, 0, 7]);
      assert.equal(7, a[nblas.iamax(a)]);
    });
  });

  describe('?gbmv', function () {
    // computes matrix-vector product using a general band matrix
  });

  describe('?gemv', function () {
    // computes a matrix-vector product using a general matrix
  });

  describe('?ger', function () {
    // performs a rank-1 update of a general matrix
  });

  describe('?sbmv', function () {
    // computes a matrix-vector product using a symmetric band matrix
  });

  describe('?spmv', function () {
    // computes a matrix-vector product using a symmetric packed matrix
  });

  describe('?spr', function () {
    // performs a rank-1 update of a symmetric packed matrix
  });

  describe('?spr2', function () {
    // performs a rank-2 update of a symmetric packed matrix
  });

  describe('?symv', function () {
    // computes a matrix-vector product for a symmetric matrix
  });

  describe('?syr', function () {
    // performs a rank-1 update of a symmetric matrix
  });

  describe('?syr2', function () {
    // performs a rank-2 update of a symmetric matrix
  });

  describe('?tbmv', function () {
    // computes a matrix-vector product using a triangular band matrix
  });

  describe('?tbsv', function () {
    // solves a system of linear equations whose coefficients are in a triangular band matrix
  });

  describe('?tpmv', function () {
    // computes a matrix-vector product using a triangular band matrix
  });

  describe('?tpsv', function () {
    // solves a system of linear equations whose coefficients are in a triangular packed matrix
  });

  describe('?trmv', function () {
    // computes a matrix-vector product using a triangular matrix
  });

  describe('?trsv', function () {
    // solves a system of linear equations whose coefficients are in a triangular matrix
  });

  describe('?gemm', function () {
    it('works for 3x1 * 1x3', function () {
      // computes a matrix-matrix product

      // a is 1x3 matrix
      var a = new Float64Array([
        1, 2, 3
      ]);

      // b is 3x1 matrix
      var b = new Float64Array([
        2,
        3,
        4
      ]);

      // c will hold 3x3 matrix
      var c = new Float64Array(9);

      var ans = new Float64Array([
        2, 3, 4,
        4, 6, 8,
        6, 9, 12
      ]);

      nblas.gemm(a, b, c, 3, 3, 1);
      assert.deepEqual(ans, c);
    });

    it('works for 2x2 * 2x2', function () {
      // a is 2x2 matrix
      var a = new Float64Array([
        1, 2,
        3, 4
      ]);

      // b is 2x2 matrix
      var b = new Float64Array([
        5, 6,
        7, 8
      ]);

      // c will hold 2x2 matrix
      var c = new Float64Array(4);

      var ans = new Float64Array([
        19, 22,
        43, 50
      ]);

      nblas.gemm(a, b, c, 2, 2, 2);
      assert.deepEqual(ans, c);
    });
  });

  describe('?symm', function () {
    // computes a matrix-matrix product where one input matrix is symmetric
  });

  describe('?syrk', function () {
    // performs a symmetric rank-k update
  });

  describe('?syr2k', function () {
    // performs a symmetric rank-2k update
  });

  describe('?trmm', function () {
    // computes a matrix-matrix product where one input matrix is triangular
  });

  describe('?trsm', function () {
    // solves a triangular matrix equation
  });

  describe('?gesv', function () {
    // solves for X the system of linear equations A*X = B, where A is an n-by-n matrix
    it('works for 2x2 * 2x2', function () {
      // a is 2x2 matrix
      var a = new Float64Array([
        1, 2,
        3, 4
        //1, 3, 2, 4
      ]);

      // x is 2x2 matrix, shape same as b
      var x = new Float64Array([
        5, 6,
        7, 8
        //5, 7, 6, 8
      ]);

      // b is 2x2 matrix
      var b = new Float64Array([
        19, 22,
        43, 50
        //19, 43, 22, 50
      ]);

      var res = nblas.gesv(a, b, 2, 2);
      assert.equal(res, 0);
      b = b.map(el => el.toFixed(2) );
      assert.deepEqual(b, x);
    });

    it('works for 2x2 * 1x2', function () {
      var a = new Float64Array([
        3, 1,
        1, 2
      ]);

      var x = new Float64Array([
        2, 3,
      ]);

      var b = new Float64Array([
        9, 8,
      ]);

      var res = nblas.gesv(a, b, 2, 1);
      assert.equal(res, 0);
      b = b.map(el => el.toFixed(2) );
      assert.deepEqual(b, x);
    });

    it('works for 5x5 * 5x3', function () {
      var a = new Float64Array([
        6.80,   -6.05,  -0.45,   8.32,  -9.67, 
        -2.11,  -3.30,   2.58,   2.71,  -5.14, 
        5.66,    5.36,  -2.70,   4.35,  -7.26,
        5.97,   -4.44,   0.27,  -7.17,   6.08,
        8.23,    1.08,   9.04,   2.14,  -6.87
      ]);
      var b = new Float64Array([
        4.02,   -1.56,   9.81,
        6.19,    4.00,  -4.09,
        -8.22,  -8.67,  -4.57,
        -7.57,   1.75,  -8.61,
        -3.03,   2.86,   8.99
      ]);
      var x = new Float64Array([
        -0.80,  -0.39,   0.96,
        -0.70,  -0.55,   0.22,
        0.59,    0.84,   1.90,
        1.32,   -0.10,   5.36,
        0.57,    0.11,   4.04
      ]);

      var res = nblas.gesv(a, b, 5, 3);
      assert.equal(res, 0);
      b = b.map(el => el.toFixed(2) );
      assert.deepEqual(b, x);
    });

  });
}());
