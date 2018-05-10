#include <iostream>
#include <octave/oct.h>
#include <octave/builtin-defun-decls.h>

DEFUN_DLD (lfsr_gen_cc,
           args,
           ,
           "LFSR sequence generation\n"
           "s = lfsr_gen_cc(state,taps,n);")
{
  octave_value_list retval;

  const int nargin(args.length());
  if (nargin != 3) {
    print_usage();
    return retval;
  }

  const boolNDArray state(args(0).array_value());
  const boolNDArray taps(args(1).array_value());
  const int            n(args(2).int_value());

  const int m = taps.numel();
  if (m != state.numel())
      error("dim mismatch");

  boolNDArray result(dim_vector(1,n));
  for (int i=0; i<taps.numel(); ++i)
    result(i) = state(i);

  for (int i=m; i<n; ++i) {
    int sum=0;
    for (int j=0; j<m; ++j)
      sum += result(i-m+j) * taps(j);
    result(i) = sum & 1; // using AND instead of mod
  }
  retval(0) = result;
  return retval;
}
