#include <iostream>
#include <octave/oct.h>

// x = gauss_mod2_cc(a,b)
// b = a*x

DEFUN_DLD (gauss_mod2_cc,
           args,
           ,
           "gaussian elimination mod2\n"
           "x = gauss_mod2_cc(a,b); ## this solves a*x = b")
{
  octave_value_list retval;
  retval(1) = boolMatrix(0,0);
  retval(0) = boolMatrix(0,0);

  const int nargin(args.length());
  if (nargin != 1 && nargin != 2) {
    print_usage();
    return retval;
  }

  boolMatrix a(args(0).array_value());
  boolMatrix x(nargin == 2 ? args(1).array_value() : boolMatrix(0,0));
  if (nargin == 2 && (a.dim1() != x.dim1() || a.dim1() != a.dim2() || x.dim2() != 1))
    error("dim mismatch");

  // (1) upper-triangular form
  const int n(a.dim1());
  const int m(a.dim2());
  for (int k=0; k<m; ++k) {
    int idx(-1);
    for (int l=k; l<n; ++l) {
      if (a(l,k)) {
        idx = l;
        break;
      }
    }
    if (idx < 0)
      return retval;

    // permute k <-> idx
    if (idx != k) {
      if (nargin == 2)
        std::swap(x(k), x(idx));
      for (int l=0; l<m; ++l)
        std::swap(a(k,l), a(idx,l));
    }

    // a(k,k)==1
    for (int j=k+1; j<n; ++j) {
      if (a(j,k)) {
        if (nargin == 2)
          x(j) ^= x(k);
        for (int l=k; l<m; ++l)
          a(j,l) ^= a(k,l);
      }
    }
  }

  if (nargin == 2) {
    //  ## (2) back-substitution
    // x(n-1) = x(n-1);
    for (int i=n-2; i>=0; --i) {
      bool sum(x(i));
      for (int j=i+1; j<n; ++j)
        sum ^= (a(i,j) & x(j));
      x(i) = sum;
    }
  }

  if (nargin == 2)
    retval(1) = x;

  retval(0) = a;
  return retval;
}
