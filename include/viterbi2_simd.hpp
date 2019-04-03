// -*- C++ -*-

#ifndef _VITERBI2_SIMD_HPP_
#define _VITERBI2_SIMD_HPP_

#include <cassert>
#include <array>
#include <algorithm>
#include <bitset>
#include <vector>

#include <arm_neon.h>

//#ifdef _AARCH64_NEON_H_

// AARCh64 SIMD "half butterfly" soft-decision viterbi decoder
template<int N>
class viterbi2 {
public:
  enum { M = 1<<N };

  typedef std::vector<uint16_t>   vec_type;
  typedef std::array<uint16_t, M> arr_type;

  viterbi2(int len)
    : _poly{0x6d, 0x4f}
    , _decisions(len<<N)
    , _metric()
    , _bits0()
    , _bits1()
    , _prev()
    , _last_max_metric(0){
    make_tables();
    reset();
  }

  void reset() {
    std::fill_n(_metric.begin(), _metric.size(), 0);
    _last_max_metric = 0;
  }

  vec_type::iterator decision(int i) { return _decisions.begin()+(i<<N); }
  vec_type::const_iterator decision(int i) const { return _decisions.begin()+(i<<N); }

  void update(int j, std::uint8_t sym0, std::uint8_t sym1) {
    arr_type new_metric;
    auto jdec = decision(j);

    uint8_t const *pb0    = &_bits0[0];
    uint8_t const *pb1    = &_bits1[0];
    uint16_t *pnew_metric = &new_metric[0];
    uint16_t *pmetric0    = &_metric[0];
    uint16_t *pmetric1    = &_metric[0]+(1<<(N-1));
    uint16_t *pdec        = &jdec[0];

    uint16x8_t mmin = vdupq_n_u16(65535);    
    uint8x8_t  s0   = vdup_n_u8(sym0 ^ 255);
    uint8x8_t  s1   = vdup_n_u8(sym1 ^ 255);
    for (int i=0; i<M/8; ++i) {
      uint8x8_t b0 = vld1_u8(pb0);
      uint8x8_t b1 = vld1_u8(pb1);
      pb0 += 8;
      pb1 += 8;
      // new_metric[i]  = (_bits0[i] ^ s0) + (_bits1[1] ^ s1);
      uint16x8_t m = vaddl_u8(veor_u8(b0, s0), veor_u8(b1, s1));

      // previous metrics:
      uint16x4_t pm0 = vld1_u16(pmetric0);
      uint16x4_t pm1 = vld1_u16(pmetric1);
      pmetric0 += 4;
      pmetric1 += 4;

      // jdec[i]        = _metric[_prev[i][0]] < _metric[_prev[i][1]];
      uint16x4_t dec = vclt_u16(pm0, pm1);
      vst1q_u16(pdec, make_4to8(dec));
      pdec += 8;

      // new_metric[i] += _metric[_prev[i][jdec[i]]];
      uint16x4_t delta_new_metric = vbsl_u16(dec, pm1, pm0);
      
      m = vaddq_u16(m, make_4to8(delta_new_metric));
      vst1q_u16(pnew_metric, m);
      pnew_metric += 8;

      mmin = vminq_u16(mmin, m);
    }
    // avoid path metric overflow    
    uint16_t const imin = vminvq_u16(mmin);
    if (imin > (1<<15)) {
      _last_max_metric -= imin;
      for (int i=0; i<M; ++i)
	new_metric[i] -= imin;
    } 
    std::copy(new_metric.begin(), new_metric.end(), _metric.begin());
  }

  float chainback(vec_type& v) {
    assert(v.size() == (_decisions.size()>>N));

    auto imax = std::max_element(_metric.begin(), _metric.end());
    int  idx_max = std::distance(_metric.begin(), imax);
    int  max_metric = *imax;
    for (int k=_decisions.size()>>N; k!=0; --k) {
      v[k-1] = idx_max&1;
      //idx_max = _prev[idx_max][decision(k-1)[idx_max] != 0];
      idx_max = (idx_max>>1) + (decision(k-1)[idx_max] != 0)*(1<<(N-1));
    }
    float const quality = float(max_metric -  _last_max_metric)/255.0f;
    _last_max_metric = max_metric;
    return quality;
  }
protected:
  static uint16x8_t make_4to8(uint16x4_t x) {
    // (0,1,2,3) x (0,1,2,3) -> (0,0,1,1,2,2,3,3)
    uint16x4x2_t tmp = vzip_u16(x, x);
    return vcombine_u16(tmp.val[0], tmp.val[1]);
  }
  
  void make_tables() {
    for (int i=0, n=M; i<n; ++i) {
      const std::bitset<N> b[2] = { i&_poly[0], i&_poly[1] };
      _bits0[i] = 255*(b[0].count()%2);
      _bits1[i] = 255*(b[1].count()%2);
    }
    for (int i=0; i<M; ++i) {
      _prev[i][0] = (i>>1);
      _prev[i][1] = _prev[i][0] + (1<<(N-1));
    }
  }
private:
  const int _poly[2];
  vec_type  _decisions;
  arr_type  _metric;
  uint8_t   _bits0[M];
  uint8_t   _bits1[M];
  int       _prev[M][2];
  int       _last_max_metric;
} ;

//#endif // _AARCH64_NEON_H_
#endif // _VITERBI2_SIMD_HPP_
