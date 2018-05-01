// -*- C++ -*-

#ifndef _S4285_S4285_CHANNEL_ESTIMATOR_HPP_
#define _S4285_S4285_CHANNEL_ESTIMATOR_HPP_

#include <array>
#include <bitset>

#include "s4285_frame_detector.hpp"
#include "constellation.hpp"

namespace s4285 {

class channel_estimator : public frame_detector {
public:
  typedef std::vector<std::uint8_t> uint8_vector_type;

  channel_estimator(int taps_data, int taps_symbols, int mod)
    : frame_detector()
    , _scramble(gen_scramble())
    , _symbols(256)
    , _soft_decisions(128*mod/2)
    , _taps_data(sps()*2*taps_data+1) // 5 SPS times (past, future)
    , _taps_symbols(taps_symbols)
    , _known_symbols(taps_symbols)
    , _modulation_type(mod)
    , _c2({0,1})
    , _c4({0,1,3,2})
    , _c8({1,0,2,3,6,7,5,4})
  {
    reset();
  }
  virtual ~channel_estimator() {}

  int modulation_type() const { return _modulation_type; }

  virtual void reset() {
    std::cout << "LOCK channel_estimator::reset\n";
    std::fill_n(_taps_data.begin(),     _taps_data.size(),     0);
    std::fill_n(_taps_symbols.begin(),  _taps_symbols.size(),  0);
    _known_symbols.reset();
  }

  complex_type taps_data(int i)    const { return _taps_data[i]; }
  complex_type taps_symbols(int i) const { return _taps_symbols[i]; }
  std::size_t n_taps_data()    const { return _taps_data.size(); }
  std::size_t n_taps_symbols() const { return _taps_symbols.size(); }

  complex_type filter_sample(int i) const {
    complex_type s = 0;
    for (int j=0; j<n_taps_data(); ++j)
      s += z(sps()*i-n_taps_data()/2+j) * _taps_data[j];
    for (int j=0, n=n_taps_symbols(); j<n; ++j)
      s += _known_symbols[j] * _taps_symbols[j];
    return s;
  }
  complex_type constrain_sample(complex_type s) {
    switch (_modulation_type) {
      case 2: return _c2.map(s);
      case 4: return _c4.map(s);
      case 8: return _c8.map(s);
      default: return 0;
    }
  }
  void filter_update(int i, complex_type s, complex_type known_symbol, bool save_symbol=true) {
    const complex_type err = s - known_symbol;
// //    if (std::abs(err) > 0.5)
//     std::cout << "err " << i << " "<< std::abs(err) << " " << s << " " << known_symbol <<" " << _symbols[i] << std::endl;

    // update weights
    const float mu    = 0.05; // was: 0.04
    const float alpha = 0.0005f/mu;
    for (int j=0; j<n_taps_data(); ++j)
      _taps_data[j] -= mu*err*std::conj(z(sps()*i-n_taps_data()/2+j)) + 0.0f*_taps_data[j];

    for (int j=0; j<n_taps_symbols(); ++j)
      _taps_symbols[j] -= mu*err*std::conj(_known_symbols[j]) + mu*alpha*_taps_symbols[j];
    // save the known symbol
    if (save_symbol)
      _known_symbols.insert(known_symbol);
  }

  template<typename T>
  void shift_vector(std::vector<T>& v, int shift) {
    if (shift < 0) {
      for (int i=0, n=v.size(); i<n; ++i)
        v[i] = (i-shift < n ? v[i-shift] : 0);
    }
    if (shift > 0) {
      for (int i=v.size()-1; i>0; --i)
        v[i] = (i-shift > 0 ? v[i-shift] : 0);
    }
  }
  void shift_filter(int shift) {
    shift_vector(_taps_data,     -shift);
    // shift_vector(_taps_symbols,  shift/sps());
    _known_symbols.shift(-shift/sps());
  }

  complex_type symbol(int i) const { return _symbols[i]; }
  std::uint8_t soft_decision(int i) const { return _soft_decisions[i]; }

  virtual bool process_frame(int shift) {
    if (shift)
      std::cout << "process_frame shift=" << shift << std::endl;
    shift_filter(shift);
    const int n_iter_known=1;

    int num_mismatch_synch = 0;
    for (int i=0, j=0; i<256; ++i) {
      complex_type s = 0;

      if (i < 80) {              // synchronization sequence
        for (int l=0; l<n_iter_known; ++l) {
          s = filter_sample(i);
          num_mismatch_synch += ((std::real(s) > 0) != (synch(i) > 0));
          filter_update(i, s, synch(i), l==n_iter_known-1);
        }
        _symbols[i] = s;
      } else {
        if (((i-80)%48) >= 32) { // reference symbols
          for (int l=0; l<n_iter_known; ++l) {
            s = filter_sample(i);
            filter_update(i, s, _scramble[i-80], l==n_iter_known-1);
          }
        } else {                 // data symbols
          s = filter_sample(i);
          filter_update(i, s, _scramble[i-80]*constrain_sample(std::conj(_scramble[i-80])*s));
          assert(j<128*modulation_type()/2);
          switch (_modulation_type) {
            case 2: {
              _soft_decisions[j++] = _c2.soft_decision(0);
              break;
            }
            case 4: {
              _soft_decisions[j++] = _c4.soft_decision(1);
              _soft_decisions[j++] = _c4.soft_decision(0);
              break;
            }
            case 8: {
              _soft_decisions[j++] = _c8.soft_decision(2);
              _soft_decisions[j++] = _c8.soft_decision(1);
              _soft_decisions[j++] = _c8.soft_decision(0);
              // here comes the erasure
              _soft_decisions[j++] = 127;
              break;
            }
          }
        }
        s *= conj(_scramble[i-80]);
        _symbols[i] = s;

        std::cout << "__symb " << i << " " << _symbols[i] << std::endl;
        if (i == 255) {
          assert(j == 128*modulation_type()/2);
          process_frame_data(_soft_decisions.begin(), _soft_decisions.end());
        }
      }
    }
    if (num_mismatch_synch)
      std::cout << "LOCK NUM_MISMATCH_SYNC " << frame_num() << " " << num_mismatch_synch << std::endl;
    return (num_mismatch_synch < 30);
  }

  virtual void process_frame_data(uint8_vector_type::const_iterator i0,
                                  uint8_vector_type::const_iterator i1) {}

protected:
  static complex_vector_type gen_scramble() {
    bool state[9] = { 1,1,1,1,1,1,1,1,1 };
  // bool taps[9] = { 0,0,0,0,1,0,0,0,1 };
    complex_vector_type scramble(176);
    for (int i=0; i<176; ++i) {
      scramble[i] = std::exp(complex_type(0.0f, float(2*M_PI*0.125*(4*state[6]+2*state[7]+state[8]))));
      for (int j=0; j<3; ++j) {
        const bool b = state[4] ^ state[8];
        for (int k=8; k>0; --k)
          state[k] = state[k-1];
        state[0] = b;
      }
    }
    return scramble;
  }
private:
  complex_vector_type _scramble;
  complex_vector_type _symbols;
  uint8_vector_type   _soft_decisions;
  complex_vector_type _taps_data;
  complex_vector_type _taps_symbols;
  complex_circular_buffer_type _known_symbols;
  int                 _modulation_type; // 2,4,8
  constellation<1>    _c2;
  constellation<2>    _c4;
  constellation<3>    _c8;

} ;

} // namespace s4285

#endif // _S4285_CHANNEL_ESTIMATOR_HPP_
