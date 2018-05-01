#include <iostream>
#include <vector>
#include <complex>

#include <octave/oct.h>
#include <octave/ov-struct.h>

DEFUN_DLD (pll_cc,
           args,
           ,
           "PLL\n"
           "AAA")
{
  octave_value_list retval;

  const int nargin(args.length());
  if (nargin != 2) {
    print_usage();
    return retval;
  }

  octave_scalar_map state(args(0).scalar_map_value());

  double theta = state.contents("theta").double_value();
  double f1    = state.contents("f1").double_value();
  double ts    = state.contents("ts").double_value();
  double ud    = state.contents("ud").double_value();
  double uf    = state.contents("uf").double_value();
  double fc    = state.contents("fc").double_value();
  double dwl   = state.contents("dwl").double_value();
  // double limf1 = state.contents("limf1").double_value();
  NDArray a    = state.contents("a").array_value();
  NDArray b    = state.contents("b").array_value();

  ComplexNDArray z = args(1).complex_array_value();
  const int      n = z.numel();
  dim_vector     dv(n, 1);
  NDArray        vec_theta(dv), vec_f1(dv);
  ComplexNDArray vec_carr(dv);

  for (int i=0; i<n; ++i) {

    // phase update
    theta += f1*ts;

    // phase detector
    const std::complex<double> carr(std::exp(std::complex<double>(0, theta)));
    const double ud_old = ud;
    ud = std::arg(z(i)*std::conj(carr));

    // loop filter (1 + s*tau1) / (s*tau2)
    uf = -a(1)*uf + b(0)*ud + b(1)*ud_old;

#if 0
    const double limits[2] = {
      2*M_PI*(-dwl*5),
      2*M_PI*(+dwl*5)
    };
    if (uf < limits[0])
      uf = limits[0];
    if (uf > limits[1])
      uf = limits[1];
#endif
    // nco frequency update
    f1 = 2*M_PI*fc + uf;

    vec_f1(i)    = f1;
    vec_theta(i) = theta;
    vec_carr(i)  = carr;
  }

  state.contents("theta") = theta;
  state.contents("ud")    = ud;
  state.contents("uf")    = uf;
  state.contents("f1")    = f1;

  retval(2) = octave_value(vec_f1);
  retval(1) = octave_value(vec_carr);
  retval(0) = octave_value(vec_theta);
  return retval;
}
