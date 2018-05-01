#include <iostream>
#include <vector>

#include <octave/oct.h>
#include <octave/ov-struct.h>

#include "s4285_bitstream_decoder.hpp"

class s4285_decoder_oct : public s4285::bitstream_decoder {
public:
  s4285_decoder_oct(int taps_data, int taps_symbols, int mod)
    : s4285::bitstream_decoder(taps_data, taps_symbols, mod)
    , _frame_counter(0)
    , _chunk_counter(0)
    , _oct_frame_symb()
    , _oct_taps_data()
    , _oct_taps_symbols()
    , _oct_soft_decisions()
    , _oct_bits()
    , _oct_chunk_counter() {}

  virtual ~s4285_decoder_oct() {}

  virtual void reset() {
    std::cout << "LOCK s4285_decoder_oct::reset " << _chunk_counter << std::endl;
    s4285::bitstream_decoder::reset();
    _chunk_counter += 1;
  }
  virtual void process_bits(std::vector<int>::const_iterator i0,
                            std::vector<int>::const_iterator i1) {
    ComplexNDArray a(dim_vector(1,256));
    for (int j=0; j<256; ++j)
      a(j) = symbol(j);
    _oct_frame_symb(_frame_counter) = a;

    ComplexNDArray b(dim_vector(1,n_taps_data()));
    for (int j=0; j<n_taps_data(); ++j)
      b(j) = taps_data(j);
    _oct_taps_data(_frame_counter) = b;

    ComplexNDArray c(dim_vector(1,n_taps_symbols()));
    for (int j=0; j<n_taps_symbols(); ++j)
      c(j) = taps_symbols(j);
    _oct_taps_symbols(_frame_counter) = c;

    uint8NDArray e(dim_vector(1, 128*mod2n(modulation_type())));
    for (int j=0; j<e.numel(); ++j)
      e(j) = soft_decision(j);
    _oct_soft_decisions(_frame_counter) = e;

    uint8NDArray f(dim_vector(1, std::distance(i0, i1)));
    for (int j=0; j<std::distance(i0, i1); ++j)
      f(j) = i0[j];
    _oct_bits(_frame_counter) = f;

    _oct_chunk_counter(_frame_counter) = octave_value(_chunk_counter);
    _frame_counter += 1;
  }

  octave_value to_oct() const {
    octave_map map;
    map.setfield("symb",           _oct_frame_symb);
    map.setfield("taps_data",      _oct_taps_data);
    map.setfield("taps_symb",      _oct_taps_symbols);
    map.setfield("soft_decisions", _oct_soft_decisions);
    map.setfield("bits",           _oct_bits);
    map.setfield("chunk_counter",  _oct_chunk_counter);
    return map;
  }
protected:
private:
  int                 _frame_counter;
  int                 _chunk_counter;
  octave_value_list   _oct_frame_symb;
  octave_value_list   _oct_taps_data;
  octave_value_list   _oct_taps_symbols;
  octave_value_list   _oct_soft_decisions;
  octave_value_list   _oct_bits;
  octave_value_list   _oct_chunk_counter;
} ;

DEFUN_DLD (s4285_cc,
           args,
           ,
           "STANAG4285 decoder\n"
           "s4285_cc(z,mod)")
{
  octave_value_list retval;

  const int nargin(args.length());
  if (nargin < 1) {
    print_usage();
    return retval;
  }

  const int mod = (nargin == 2 ? args(1).int_value() : 2);

  s4285_decoder_oct d(12, 4, mod); // was: 15,6,2

  const FloatComplexNDArray z = args(0).complex_array_value();
  for (int i=0, n=z.numel(); i<n; ++i)
    d.process(z(i));

  retval(0) = d.to_oct();
  return retval;
}
