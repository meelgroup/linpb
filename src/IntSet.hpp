/***********************************************************************
Copyright (c) 2014-2020, Jan Elffers
Copyright (c) 2019-2020, Jo Devriendt
Copyright (c) 2020, Stephan Gocht

Parts of the code were copied or adapted from MiniSat.

MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010  Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***********************************************************************/

#pragma once

#include "aux.hpp"

namespace rs {

struct IntSet {  // TODO: template to long long, int128, ...?
 private:
  std::vector<int> _index = {-1};
  std::vector<int>::iterator index = _index.begin();
  static constexpr int _unused_() { return -1; }
  const int resize_factor = 2;

  bool check() {
    for (int i = 0; i < (int)_index.size() / 2; ++i) {
      assert(index[i] == _unused_() || i == keys[index[i]]);
      assert(index[-i] == _unused_() || -i == keys[index[-i]]);
    }
    for (int i = 0; i < (int)keys.size(); ++i) assert(index[keys[i]] == i);
    return true;
  }

 public:
  std::vector<int> keys;

  IntSet() {}
  IntSet(int size, const std::vector<int>& ints) {
    resize(size);
    for (int i : ints) add(i);
  }

  void resize(int size) { aux::resizeIntMap(_index, index, size, resize_factor, _unused_()); }
  size_t size() const { return keys.size(); }
  bool isEmpty() const { return size() == 0; }

  void clear() {
    assert(check());  // TODO: disable
    for (int k : keys) index[k] = _unused_();
    keys.clear();
  }

  bool has(int key) const { return _index.size() > (unsigned int)2 * std::abs(key) && index[key] != _unused_(); }

  void add(int key) {
    if (_index.size() <= (unsigned int)2 * std::abs(key)) resize(std::abs(key));
    if (index[key] != _unused_()) return;
    assert(!aux::contains(keys, key));
    index[key] = keys.size();
    keys.push_back(key);
  }

  void remove(int key) {
    if (!has(key)) return;
    int idx = index[key];
    index[keys.back()] = idx;
    aux::swapErase(keys, idx);
    index[key] = _unused_();
    assert(!has(key));
  }

  std::ostream& operator<<(std::ostream& o) const {
    for (int k : keys)
      if (has(k)) o << k << " ";
    return o;
  }
};

}  // namespace rs
