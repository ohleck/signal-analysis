// -*- C++ -*-

#ifndef _VITERBI2_SIMPLE_HPP_
#define _VITERBI2_SIMPLE_HPP_

#include <cassert>
#include <array>
#include <algorithm>
#include <bitset>
#include <vector>


// "half butterfly" soft-decision viterbi decoder
template<int K>
class viterbi2 {
public:
  enum { M = 1<<K };

  typedef std::vector<int>   vec_type;
  typedef std::array<int, M> arr_type;

  viterbi2(int len)
    : _poly{0x6d, 0x4f}
    , _decisions(len<<K)
    , _metric()
    , _bits()
    , _prev()
    , _last_max_metric(0) {
    make_tables();
    reset();
  }

  void reset() {
    std::fill_n(_metric.begin(), _metric.size(), 0);
    _last_max_metric = 0;
  }

  vec_type::iterator decision(int i) { return _decisions.begin()+(i<<K); }
  vec_type::const_iterator decision(int i) const { return _decisions.begin()+(i<<K); }

  void update(int j, std::uint8_t sym0, std::uint8_t sym1) {
    const int s0 = sym0 ^ 255;
    const int s1 = sym1 ^ 255;

    arr_type new_metric;
    auto jdec = decision(j);

    int mmin = 65535;

    for (int i=0; i<M; ++i) {
      new_metric[i]  = (_bits[i][0] ^ s0) + (_bits[i][1] ^ s1);
      jdec[i]        = _metric[_prev[i][0]] < _metric[_prev[i][1]];
      new_metric[i] += _metric[_prev[i][jdec[i]]];
      mmin = std::min(mmin, new_metric[i]);
    }
    // avoid path metric overflow
    if (mmin > (1<<15)) {
      _last_max_metric -= mmin;
      for (int i=0; i<M; ++i)
	new_metric[i] -= mmin;
    }
    std::copy(new_metric.begin(), new_metric.end(), _metric.begin());
  }

  float chainback(vec_type& v) {
    assert(v.size() == (_decisions.size()>>K));

    auto imax = std::max_element(_metric.begin(), _metric.end());
    int idx_max = std::distance(_metric.begin(), imax);
    int max_metric = *imax;
    for (int k=_decisions.size()>>K; k!=0; --k) {
      v[k-1] = idx_max&1;
      // idx_max = _prev[idx_max][decision(k-1)[idx_max]];
      idx_max = (idx_max>>1) + (decision(k-1)[idx_max] != 0)*(1<<(K-1));
    }
    float const quality = float(max_metric -  _last_max_metric)/255.0;
    _last_max_metric = max_metric;
    return quality;
  }
protected:
  void make_tables() {
    for (int i=0, n=M; i<n; ++i) {
      const std::bitset<K> b[2] = { i&_poly[0], i&_poly[1] };
      _bits[i][0] = 255*(b[0].count()%2);
      _bits[i][1] = 255*(b[1].count()%2);
    }
    for (int i=0; i<M; ++i) {
      _prev[i][0] = (i>>1);
      _prev[i][1] = _prev[i][0] + (1<<(K-1));
    }
  }
private:
  const int _poly[2];
  vec_type  _decisions;
  arr_type  _metric;
  int       _bits[M][2];
  int       _prev[M][2];
  int       _last_max_metric;
} ;

#endif // _VITERBI2_SIMPLE_HPP_
