// -*- C++ -*-

#ifndef _VITERBI2_HPP_
#define _VITERBI2_HPP_

#include <cassert>
#include <algorithm>
#include <array>
#include <bitset>
#include <vector>

// soft-decision viterbi decoder
// based on Phil Karn's libfec
template<std::size_t N>
class viterbi2 {
public:
  enum { M = 1<<(N-1) };

  typedef std::vector<int>   vec_type;
  typedef std::array<int, M> arr_type;

  viterbi2(int len)
    : _poly{0x6d, 0x4f}
    , _decisions(len<<(N-1))
    , _metric()
    , _bits()
    , _prev()
    , _last_max_metric(0)
  {
    make_tables();
    reset();
  }

  void reset() {
    std::fill_n(_metric.begin(), _metric.size(), 0);
    _last_max_metric = 0;
  }

  vec_type::iterator decision(int i) { return _decisions.begin()+(i<<(N-1)); }
  vec_type::const_iterator decision(int i) const { return _decisions.begin()+(i<<(N-1)); }

  void update(int j, std::uint8_t sym0, std::uint8_t sym1) {
    const int s0 = sym0 ^ 255;
    const int s1 = sym1 ^ 255;

    arr_type new_metric;
    auto jdec = decision(j);

    int mmin[2] = {65535, 65535};
    for (int i=0; i<M; i+=2) {
      const int p0 = _prev[i][0];
      const int p1 = _prev[i][1];

      const int m0[2] = {
        _metric[p0] + (_bits[p0][0][0] ^ s0) + (_bits[p0][0][1] ^ s1),
        _metric[p1] + (_bits[p1][0][0] ^ s0) + (_bits[p1][0][1] ^ s1)
      };
      const int m1[2] = {
        _metric[p0] + (_bits[p0][1][0] ^ s0) + (_bits[p0][1][1] ^ s1),
        _metric[p1] + (_bits[p1][1][0] ^ s0) + (_bits[p1][1][1] ^ s1)
      };

      jdec[i  ] = m0[0] < m0[1];
      jdec[i+1] = m1[0] < m1[1];

      new_metric[i  ] = m0[jdec[i  ]];
      new_metric[i+1] = m1[jdec[i+1]];
      mmin[0] = std::min(mmin[0], new_metric[i]);
      mmin[1] = std::min(mmin[1], new_metric[i+1]);
    }
    // avoid path metric overflow
    int const imin = std::min(mmin[0], mmin[1]);
    if (imin > (1<<15)) {
      _last_max_metric -= imin;
      for (int i=0; i<M; ++i)
	new_metric[i] -= imin;
    }
    std::copy(new_metric.begin(), new_metric.end(), _metric.begin());
  }

  float chainback(vec_type& v) {
    assert(v.size() == (_decisions.size()>>(N-1)));

    auto imax = std::max_element(_metric.begin(), _metric.end());
    int idx_max = std::distance(_metric.begin(), imax);
    int max_metric = *imax;
    for (int k=_decisions.size()>>(N-1); k!=0; --k) {
      v[k-1] = idx_max&1;
      //      idx_max = _prev[idx_max][decision(k-1)[idx_max]];
      idx_max = (idx_max>>1) + (decision(k-1)[idx_max] != 0)*(1<<(N-1));
    }
    float const quality = float(max_metric -  _last_max_metric)/255.0;
    _last_max_metric = max_metric;
    return quality;
  }
protected:
  void make_tables() {
    for (int i=0, n=1<<N; i<n; ++i) {
      const std::bitset<N> b[2] = { _poly[0]&i, _poly[1]&i };
      _bits[i>>1][i&1][0] = 255*(b[0].count()%2);
      _bits[i>>1][i&1][1] = 255*(b[1].count()%2);
    }
    for (int i=0; i<M; ++i) {
      _prev[i][0] = (i>>1);
      _prev[i][1] = _prev[i][0] + (1<<(N-2));
    }
  }
private:
  const int _poly[2];
  vec_type  _decisions;
  arr_type  _metric;
  int       _bits[M][2][2];
  int       _prev[M][2];
  int       _last_max_metric;
} ;

#endif // _VITERBI2_HPP_
