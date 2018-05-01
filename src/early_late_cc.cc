#include <iostream>
#include <vector>

#include <octave/oct.h>
#include <octave/ov-struct.h>

DEFUN_DLD (early_late_cc,
           args,
           ,
           "early late synchronization\n"
           "AAA")
{
  octave_value_list retval;

  const int nargin(args.length());
  if (nargin != 2) {
    print_usage();
    return retval;
  }

  octave_scalar_map state(args(0).scalar_map_value());
  double  alpha   = state.contents("alpha").double_value();
  double  period  = state.contents("period").double_value();
  double  t       = state.contents("t").double_value();
  NDArray history = state.contents("history").array_value();
  int     counter = state.contents("counter").int_value();
  double  err     = state.contents("err").double_value();

  double *h = history.fortran_vec();
  int m = history.numel();

  std::vector<int> bits;

  NDArray z = args(1).array_value();
  for (int i=0, n=z.numel(); i<n; ++i) {
    if (counter == std::floor(t)) {
      double sum=0;
      for (int j=0, nj=int(period+0.5); j<nj; ++j)
        sum += h[j];
      bits.push_back(2*(sum > 0.0)-1);

      double early=0, late=0;
      for (int j=0; j<m/2; ++j) {
        early += h[j];
        late  += h[j+1+m/2];
      }
      err      = (1-alpha)*err - alpha*(std::abs(early)-std::abs(late))/m*2;
      t       += period - counter + 0.5*err;
      counter  = 0;
    }
    for (int j=0; j<m-1; ++j)
      h[j] = h[j+1];
    h[m-1] = z(i);
    counter +=1;
  }

  boolNDArray bbits(dim_vector(1,bits.size()), 0);
  for (int i=0, n=bits.size(); i<n; ++i)
    bbits(i) = bits[i]>0;
  retval(1) = bbits;

  state.contents("counter") = counter;
  state.contents("history") = history;
  state.contents("t")       = t;
  state.contents("err")     = err;
  retval(0) = state;
  return retval;
}
