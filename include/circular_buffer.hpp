// -*- C++ -*-

#ifndef _CIRCULAR_BUFFER_H_
#define _CIRCULAR_BUFFER_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

template<typename T>
class circular_buffer {
public:
  typedef T value_type;
  typedef std::vector<T> vector_type;
  typedef typename vector_type::iterator iterator;
  typedef typename vector_type::iterator const_iterator;

  circular_buffer(std::size_t size)
    : _v(2*size, value_type())
    , _size(size)
    , _i(0)
    , _is_filled(false) {}

  size_t size() const  { return _size; }

  iterator       begin()       { return _v.begin() + _i; }
  const_iterator begin() const { return const_cast<vector_type&>(_v).begin() + _i; }

  iterator       end()         { return begin() + _size; }
  const_iterator end() const   { return begin() + _size; }

  bool is_filled() const { return _is_filled; }
  std::size_t pos() const { return _i; }

  // get contents in the current valid range
  T& operator[](std::size_t idx) {
    assert(_i+idx >= 0); assert(_i+idx < _v.size());
    return _v[_i+idx];
  }
  T  operator[](std::size_t idx) const {
    assert(_i+idx >= 0); assert(_i+idx < _v.size());
    return _v[_i+idx];
  }

  // get contents using modulo operation
  T at(int idx) const {
    idx += _i;
    idx += _v.size()/2 * (idx < 0);
    idx -= _v.size()/2 * (idx >= _v.size());
    return _v[idx];
  }

  void insert(value_type v) {
    _v[_i] = _v[_i+_size] = v;
    _i = (_i+1 == _size ? 0 : _i+1);
    _is_filled = (not _is_filled ? _i == 0 : _is_filled);
    // _i points to the oldest sample in the buffer
    // _v[_i .. _i+_size-1] are guaranteed to be contiguous
  }

  void shift(int di) {
    _i  = (_i+_size + di) % _size;
    _i += (_i < 0) * _size;
  }
  void reset() {
    std::fill_n(_v.begin(), _v.size(), value_type());
    _i = 0;
    _is_filled = false;
  }
  void resize(std::size_t size) {
    _v.resize(size);
    reset();
  }
protected:

private:
  vector_type _v;
  std::size_t _size;
  std::size_t _i;
  bool        _is_filled;
} ;

#endif // _CIRCULAR_BUFFER_H_
