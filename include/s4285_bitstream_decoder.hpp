// -*- C++ -*-

#ifndef _S4285_BITSTREAM_DECODER_HPP_
#define _S4285_BITSTREAM_DECODER_HPP_

#include <vector>
#include <array>

#include "viterbi2.hpp"
#include "s4285_deinterleaver.hpp"
#include "s4285_channel_estimator.hpp"

namespace s4285 {

class bitstream_decoder : public channel_estimator {
public:
  bitstream_decoder(int taps_data, int taps_symbols, int mod)
    : channel_estimator(taps_data, taps_symbols, mod)
    , _deinterleaver(12*mod2n(mod))
    , _soft_symbol_buffer(256)
    , _i(0)
    , _viterbi_decoder(128)
    , _bits(128) {}

  virtual ~bitstream_decoder() {}

  virtual void reset() {
    std::cout << "LOCK bitstream_decoder::reset\n";
    frame_detector::reset();
    _viterbi_decoder.reset();
    _deinterleaver.reset();
  }

  // input are soft symbols
  virtual void process_frame_data(uint8_vector_type::const_iterator i0,
                                  uint8_vector_type::const_iterator i1) {
    for (; i0!=i1; i0+=32) {
      assert(i0 != i1);
      if (_deinterleaver.insert(i0, i0+32)) {
        assert( _soft_symbol_buffer.begin()+_i <  _soft_symbol_buffer.end());
        std::copy(_deinterleaver.begin(), _deinterleaver.end(),
                  _soft_symbol_buffer.begin()+_i);
        _i += 32;
        if (_i == _soft_symbol_buffer.size()) {
          // viterbi decode
          _viterbi_decoder.reset();

          for (int i=0,n=_soft_symbol_buffer.size(); i<n; i+=2)
            _viterbi_decoder.update(i/2, _soft_symbol_buffer[i], _soft_symbol_buffer[i+1]);

          const int tt = _viterbi_decoder.chainback(_bits);
          std::cout << "TT " << tt << std::endl;

          process_bits(_bits.begin(), _bits.end());

          _i = 0;
        }
      }
    }
  }

  virtual void process_bits(std::vector<int>::const_iterator i0,
                            std::vector<int>::const_iterator i1) {}
protected:
private:
  deinterleaver<std::uint8_t> _deinterleaver;
  uint8_vector_type _soft_symbol_buffer;
  int _i;
  viterbi2<7> _viterbi_decoder;
  std::vector<int> _bits;
} ;

} // namespace s4285

#endif // _S4285_BITSTREAM_DECODER_HPP_
