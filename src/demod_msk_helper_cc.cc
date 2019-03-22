#include <iostream>
#include <octave/oct.h>
#include <octave/builtin-defun-decls.h>

DEFUN_DLD (demod_msk_helper_cc,
           args,
           ,
           "coherent MSK demodulation helper\n"
           "inputs: IQ samples z and the data carrier t\n"
           "outputs: soft decisions for the 'X' and 'Y' bit streams"
           "[bX,bY] = demod_msk_helper_cc(z,t);")
{
  octave_value_list retval;

  int const nargin(args.length());
  if (nargin != 2) {
    print_usage();
    return retval;
  }

  FloatComplexNDArray const z(args(0).float_complex_array_value());
  FloatNDArray const        t(args(1).float_array_value());

  int const n = z.numel();
  if (t.numel() != n) {
    print_usage();
    return retval;
  }

  FloatNDArray bX(dim_vector(1,n), 0.0f);
  FloatNDArray bY(dim_vector(1,n), 0.0f);

  int jX=0, jY=0;
  float ct_last=0, st_last=0;
  for (int i=0; i<n; ++i) {
    float const ct = std::cos(t(i));
    float const st = std::sin(t(i));
    bX(jX) += ct * z(i).real();
    bY(jY) += st * z(i).imag();
    jX += (std::signbit(ct) != std::signbit(ct_last));
    jY += (std::signbit(st) != std::signbit(st_last));
    ct_last = ct;
    st_last = st;
  }
  bX.resize1(jX);
  bY.resize1(jY);
  retval(1) = bY;
  retval(0) = bX;
  return retval;
}
