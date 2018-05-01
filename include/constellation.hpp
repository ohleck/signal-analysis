// -*- C++ -*-

#ifndef _CONSTELLATION_HPP_
#define _CONSTELLATION_HPP_

#include <algorithm>
#include <array>
#include <complex>
#include <cmath>
#include <cstdint>
#include <vector>

template<int M>
class constellation {
public:
  typedef std::complex<float> complex_type;
  enum { N = (1<<M) };
  constellation(const std::array<int, N>& gray_code_map, int max_time_constant=120)
    : _points()
    , _soft_decision()
    , _table()
    , _max_time_constant(max_time_constant)
    , _counter(0)
    , _sigma(0)
  {
    for (int i=0; i<N; ++i)      // _points are in order of the gray_code_map
      _points[i] = std::exp(complex_type(0.0f, float(2*M_PI/N*gray_code_map[i])));

    for (int i=0; i<1024; ++i) { // llr -> probability(1)
      const float x = -7.0f + 14.0f*float(i)/1023.0f;
      _table[i] = std::uint8_t(0.5+255/(1+std::exp(x)));
    }
  }

  float soft_decision(int i) const { return _soft_decision[i]; }

  // returns the nearest point in the constellation and computes the soft decisions
  // TODO: erasures
  complex_type map(complex_type sample) {
    // (1) compute distance of sample to all points in the constellation
    std::array<float, N> distances;
    for (int i=0; i<N; ++i)
      distances[i] = std::norm(_points[i]-sample);
    const int imin = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));

    // (2) update filtered RMS of distances
    _counter = (_counter < _max_time_constant ? 1+_counter : _counter);
    const float alpha = 1.0f/_counter;
    _sigma = _sigma*(1-alpha) + 2*distances[imin]*alpha;
    std::cout << "_sigma " << _sigma << std::endl;

    // (3) compute the soft decisions based on log-likelihood ratios
    std::cout << "_sigma " << _sigma << " PROB: ";
    for (int j=0; j<M; ++j) {
      std::array<float, N/2> d[2]; // distances to all points in the constellation where bit j is '0','1'
      typename std::array<float, N/2>::iterator dj[2] = { d[0].begin(), d[1].begin() };
      for (int i=0; i<N; ++i)
        *(dj[(i>>j)&1]++) = distances[i];

      assert(dj[0] == d[0].end());
      assert(dj[1] == d[1].end());
      // log-likelihood ratio
      const float llr = 0.5*(*std::min_element(d[1].begin(), d[1].end()) -
                             *std::min_element(d[0].begin(), d[0].end()))/(_sigma*_sigma);
      // the soft decision is the probability for the bit j to be 1 (8-bit quanization 0-255)
      _soft_decision[j] = table_lookup(llr);
      std::cout << int(_soft_decision[j]) << "," << 255*1/(1+std::exp(llr)) << "," << llr << std::endl;
    }
    std::cout << std::endl;
    return _points[imin];
  }
protected:
  std::uint8_t table_lookup(float llr) const {
    const float x = (7.0f + std::max(-7.0f, std::min(+7.0f, llr)))/14.0f; // \in [0,1]
    return _table[int(0.5 + 1023*x)];
  }
private:
  std::array<complex_type, N> _points;
  std::array<std::uint8_t, M> _soft_decision;
  std::array<std::uint8_t, 1024> _table;
  int   _max_time_constant;
  int   _counter;
  float _sigma;
} ;

#endif // _CONSTELLATION_HPP_
