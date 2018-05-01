// -*- C++ -*-

#ifndef _VITERBI2_HPP_
#define _VITERBI2_HPP_

#include <array>
#include <bitset>
#include <vector>

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
  {
    make_tables();
    reset();
  }

  void reset() {
    std::fill_n(_metric.begin(), _metric.size(), 0);
  }

  vec_type::iterator decision(int i) { return _decisions.begin()+(i<<(N-1)); }
  vec_type::const_iterator decision(int i) const { return _decisions.begin()+(i<<(N-1)); }

  void update(int j, std::uint8_t sym0, std::uint8_t sym1) {
    const int s0 = sym0 ^ 255;
    const int s1 = sym1 ^ 255;

    arr_type new_metric;
    auto jdec = decision(j);

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
    }
    std::copy(new_metric.begin(), new_metric.end(), _metric.begin());
  }

  int chainback(vec_type& v) const {
    assert(v.size() == (_decisions.size()>>(N-1)));

    auto imax = std::max_element(_metric.begin(), _metric.end());
    int  idx_max = std::distance(_metric.begin(), imax);
    for (int k=_decisions.size()>>(N-1); k!=0; --k) {
      v[k-1] = idx_max&1;
      idx_max = _prev[idx_max][decision(k-1)[idx_max]];
    }
    return *imax/255;
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
} ;

#endif // _VITERBI2_HPP_
