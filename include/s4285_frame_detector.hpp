// -*- C++ -*-

#ifndef _S4285_FRAME_DETECTOR_HPP_
#define _S4285_FRAME_DETECTOR_HPP_

#include <cassert>
#include <complex>
#include <numeric>

#include "circular_buffer.hpp"

namespace s4285 {

class frame_detector {
public:
  typedef std::complex<float> complex_type;

  typedef std::vector<complex_type> complex_vector_type;
  typedef std::vector<float>           real_vector_type;

  typedef circular_buffer<complex_type> complex_circular_buffer_type;
  typedef circular_buffer<float>           real_circular_buffer_type;

  const int FRAME_LENGTH     = 256;
  const int SYNCH_LENGTH     =  80;
  const int SCRAMBLE_LENGTH  = 176;
  const int CIRCULAR_BUFFER_LENGTH = FRAME_LENGTH+SYNCH_LENGTH;

  enum class state {
    UNLOCKED,
    LOCKED
  };

  class counter {
  public:
    counter(int length)
      : _frame_length(length)
      , _frame_counter(0)
      , _sample_counter(0) {}
    void reset() {
      _frame_counter = _sample_counter = 0;
    }
    std::size_t frame_num() const { return _frame_counter; }
    std::size_t sample_num() const { return _sample_counter; }
    counter& operator++() {
      const bool b(_sample_counter+1 == _frame_length);
      _frame_counter +=  b;
      _sample_counter = (b ? 0 : _sample_counter + 1);
      return *this;
    }
    operator std::size_t() const { return _frame_length*_frame_counter + _sample_counter; }
  protected:
  private:
    std::size_t _frame_length;   // length of frame
    std::size_t _frame_counter;  // number of frame
    std::size_t _sample_counter; // nuber of sample in the current frame
  } ;

  frame_detector(int samples_per_symbol=5)
    : _sps(samples_per_symbol)
    , _state(state::UNLOCKED)
    , _unlock_counter(0)
    , _synch(make_synch())
    , _z(samples_per_symbol * CIRCULAR_BUFFER_LENGTH)
    , _z_synch_corr(samples_per_symbol * (31+20))
    , _abs_z_synch_corr_filtered(samples_per_symbol * (31+20))
    , _last_locked(false)
    , _counter(samples_per_symbol*FRAME_LENGTH)
    , _last_frame_num(-1)
    , _frame_offset(-1)
    , _frame_start(-1)
    , _df(0)
    , _phase(0)
//    , _a{1,-1} // PLL loop filter with 0.05/fs bandwidth
    , _b{0.338187046465954, -0.288839024460507}
    , _ud(0) {}

  ~frame_detector() {}

  const int sps() const { return _sps; }

  float synch(int i) const { return _synch[i]; }

  int frame_num() const { return _counter.frame_num(); }

//  bool is_locked() const { return _frame_offset >= 0; }

  // returns a sample w.r.t. start of the frame
  complex_type z(int i) const { return _z.at(i + _frame_offset); }

  virtual void reset() {}

  void process(complex_type s) {
    insert(s);
    synch_corr();
    filter_abs_synch_corr();

    bool is_locked = false;
    switch (_state) {
      case state::UNLOCKED: {
        is_locked = detect_synch_single(10.0f);
        break;
      }
      case state::LOCKED: {
        is_locked = ((_frame_start == ((_counter.sample_num() + 50) % 1280)) &&
                     (_counter.frame_num() != _last_frame_num));

        if (not _last_locked && is_locked) {
          std::cout << "LOCK TEST: " << _unlock_counter << std::endl;
          is_locked = _unlock_counter < 20;
          if (!is_locked) {
            _state = state::UNLOCKED;
            _unlock_counter = 0;
            reset();
          }
        }
        break;
      }
    }

    if (is_locked && not _last_locked) {
      const int   idx_max = find_max_corr();
      const float delta_f = compute_frequency_offset(_z_synch_corr.begin() + idx_max - 4*sps(),
                                                     _z_synch_corr.begin() + idx_max + 4*sps());
      update_pll(delta_f);

      const int frame_start_old = _frame_start;
      _frame_start  = quantize_index_sps(_counter.sample_num() + idx_max + 1);
      _frame_offset = quantize_index_sps(idx_max + 1);

      std::cout << "IS_LOCKED: " << _counter.frame_num() << " " << _counter.sample_num() << " | "
                << _frame_start << " " << _frame_offset <<  " __ " << _frame_start - frame_start_old << " idx_max=" << idx_max
                << " pos= " << _z.pos() << " ___ "
                << std::size_t(_counter)
                << std::endl;

      int shift = (_state == state::LOCKED ? _frame_start - frame_start_old : 0);
      shift += (shift < 0)*sps()*FRAME_LENGTH;
      _unlock_counter += (process_frame(shift) == false);
      _state = state::LOCKED;
      _last_frame_num = _counter.frame_num();
    }
    _last_locked = is_locked;
  }

  virtual bool process_frame(int shift) {
    return false;
  }

protected:
  int quantize_index_sps(int idx) {
    return (idx % (sps()*FRAME_LENGTH));
    const int q = 8*sps();
    return (q*int(float(idx)/q) % (sps()*FRAME_LENGTH));
  }
  void insert(complex_type sample) {
    _phase  = mod_2pi(_phase + _df);
    sample *= std::exp(complex_type(0.0f, _phase));
    _z.insert(sample);
    ++_counter;
  }

  void synch_corr() { // correlation with part of synchronization sequence
    if (not _z.is_filled())
      return;

    complex_type sum = 0;
    for (int i=0; i<31*sps(); ++i)
      sum += _z[i + (7+31+20)*sps()] * synch(7 + i/sps());

    _z_synch_corr.insert(sum);
  }

  void filter_abs_synch_corr() { // moving average filter for the correlation
    if (not _z_synch_corr.is_filled())
      return;

    const int n_window = 4*sps();
    const float f = std::accumulate(_z_synch_corr.end() - n_window,
                                    _z_synch_corr.end(), 0.0f,
                                    [](float a, complex_type b) { return a+std::abs(b); })/n_window;
    _abs_z_synch_corr_filtered.insert(f);
  }

  bool detect_synch_single(float threshold = 6.0f, int offset=0) const {
    if (not _abs_z_synch_corr_filtered.is_filled())
      return false;

    // test for synchronization: expected are two peaks separated by 31*sps() samples
    const float test[3] = {
      _abs_z_synch_corr_filtered[offset+ 10      *sps()],  // first
      _abs_z_synch_corr_filtered[offset+(10+15.5)*sps()],  // mid
      _abs_z_synch_corr_filtered[offset+(10+31)  *sps()]   // last
    };
    const bool result = (test[0] > threshold*test[1] &&
                         test[2] > threshold*test[1] &&
                         std::abs((test[0]-test[2])/(test[0]+test[2])) < 1);
    return result;
  }
  bool detect_synch(float threshold, int range) const { // detect in +- range*sps()
    bool success = false;
    int i0 = 0;
    for (int i=-range*sps(); i<=+range*sps() && !success; ++i) {
      i0      = i;
      success = detect_synch_single(threshold, i);
    }
    if (success)
      std::cout << "DET_SYNC: " << i0 << std::endl;
    return success;
  }
  int find_max_corr() const { // find maximum in 1st half
    auto imax = std::max_element(_z_synch_corr.begin(),
                                 _z_synch_corr.begin()+_z_synch_corr.size()/2,
                                 [](complex_type s1, complex_type s2) { return std::norm(s1) < std::norm(s2); });
    return std::distance(_z_synch_corr.begin(), imax);
  }

  float compute_frequency_offset(complex_vector_type::const_iterator i0,
                                 complex_vector_type::const_iterator i1) const {
    float sum_w = 0, sum_wx = 0;
    for (auto i=i0; i!=i1; ++i) {
      const complex_type xy = std::conj(i[0]) * i[31*sps()];
      const float w = std::norm(xy);
      sum_w  += w;
      sum_wx += w * std::arg(xy);
    }
    return sum_wx/sum_w/31/sps(); // frequency offset in rad/sec
  }

  void update_pll(float delta_f) {
    // carrier tracking PLL update
    if (_df == 0) { // init
      _ud = _df = -delta_f;
      return;
    }
    const float ud_old = _ud;
    _ud  = -delta_f;
    _df +=_b[0]*_ud + _b[1]*ud_old;
  }

  static float mod_2pi(float p) {
    p += 2*M_PI*(p < -M_PI);
    p -= 2*M_PI*(p >  M_PI);
    return p;
  }

  static real_vector_type make_synch() {
    std::array<bool, 5> state{ 1,1,0,1,0 };
  // bool taps[5] = { 0,0,1,0,1 };
    real_vector_type synch(80);
    for (int i=0; i<80; ++i) {
      synch[i] = -2*state[4]+1; // +1.0f -> '0', -1.0f -> '1'
      const bool b = state[2] ^ state[4];
      for (int j=4; j>0; --j)
        state[j] = state[j-1];
      state[0] = b;
    }
    return synch;
  }

private:
  const int _sps;  // samples per symbol
  state   _state;
  int     _unlock_counter;
  const real_vector_type _synch; // synchronization sequence

  // circular buffers
  complex_circular_buffer_type _z;
  complex_circular_buffer_type _z_synch_corr;
  real_circular_buffer_type    _abs_z_synch_corr_filtered;

  bool    _last_locked;
  counter _counter;        // counter modulo FRAME_LENGTH*sps()
  int     _last_frame_num; //
  int     _frame_offset;   // offset of data frame in _z;
  int     _frame_start;    // w.r.t. start of frame

  // PLL for carrier tracking
  float _df;     // frequency offset in radians per sample
  float _phase;  // accumulated phase for frequency correction
//  const float _a[2];
  const float _b[2];
  float _ud;
} ;

} // namespace s4285

#endif // _S4285_FRAME_DETECTOR_HPP_
