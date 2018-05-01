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
  retval(0) = boolMatrix(0,0);

  const int nargin(args.length());
  if (nargin != 2) {
    print_usage();
    return retval;
  }

  boolMatrix a(args(0).array_value());
  boolMatrix x(args(1).array_value());
  if (a.size(0) != x.size(0) || a.size(0) != a.size(1) || x.size(1) != 1)
    error("dim mismatch");

  // (1) upper-triangular form
  const int n(a.size(0));
  for (int k=0; k<n; ++k) {
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
      std::swap(x(k), x(idx));
      for (int l=0; l<n; ++l)
        std::swap(a(k,l), a(idx,l));
    }

    // a(k,k)==1
    for (int j=k+1; j<n; ++j) {
      if (a(j,k)) {
        x(j) ^= x(k);
        for (int l=k; l<n; ++l)
          a(j,l) ^= a(k,l);
      }
    }
  }

  //  ## (2) back-substitution
  // x(n-1) = x(n-1);
  for (int i=n-2; i>=0; --i) {
    bool sum(x(i));
    for (int j=i+1; j<n; ++j)
      sum ^= (a(i,j) & x(j));
    x(i) = sum;
  }

  retval(0) = x;
  return retval;
}
