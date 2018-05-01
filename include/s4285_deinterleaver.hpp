// -*- C++ -*-

#ifndef _S4285_DEINTERLEAVER_HPP_
#define _S4285_DEINTERLEAVER_HPP_

#include <vector>
#include <array>

namespace s4285 {

template<typename T>
class deinterleaver {
public:
  class buffer {
  public:
    buffer(int length=512)
      : _v(length+1)
      , _i(0) {}
    void clear() {
      std::fill_n(_v.begin(), _v.size(), 0);
      _i = 0;
    }
    void resize(int n) {
      _v.resize(n+1);
      clear();
    }
    void insert(T x) {
      if (_i >= _v.size())
        std::cout << "XXXX " << _i << " " << _v.size() << std::endl;
      assert(_i < _v.size());
      _v[_i] = x;
      _i = (_i+1 == _v.size() ? 0 : _i+1);
    }
    T get() const { return _v[_i]; }
  protected:
  private:
    std::vector<T> _v;
    std::size_t _i;
  } ;

  deinterleaver(int incr)
    : _b()
    , _increment(incr)
    , _counter(0)
    , _is_filled(false) { resize(incr); }

  bool is_filled() const {
    return _is_filled;
  }
  void resize(int incr) {
    _increment = incr;
    _counter   = 0;
    _is_filled = false;
    for (int i=0; i<32; ++i)
      _b[i].resize(_increment*(31-i));
  }
  void reset() {
    _counter   = 0;
    _is_filled = false;
    for (int i=0; i<32; ++i)
      _b[i].clear();
  }

  template<typename ITER>
  bool insert(ITER i0, ITER i1) {
    assert(std::distance(i0, i1) == 32);
    _is_filled = (not _is_filled ? _counter+1 == 31*_increment     : _is_filled);
    _counter   =                  (_counter+1 == 31*_increment ? 0 : _counter+1);
    for (int i=0; i<32; ++i)
      _b[i].insert(i0[i]);

    for (int i=0; i<32; ++i)
      _out[i] = _b[(9*i)&31].get(); // using AND instead of modulo

    return _is_filled;
  }

  typename std::array<T, 32>::const_iterator begin() const { return _out.begin(); }
  typename std::array<T, 32>::const_iterator end() const { return _out.end(); }
protected:

private:
  std::array<buffer, 32> _b;
  std::array<T, 32> _out;
  int _increment;
  int _counter;
  bool _is_filled;
} ;
} // namespace s4285

#endif // _S4285_DEINTERLEAVER_HPP_
