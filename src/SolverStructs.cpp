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

#include "SolverStructs.hpp"
#include "Options.hpp"
#include "globals.hpp"

namespace rs {

void ConstraintAllocator::capacity(uint32_t min_cap) {
  if (cap >= min_cap) return;

  uint32_t prev_cap = cap;
  while (cap < min_cap) {
    // NOTE: Multiply by a factor (13/8) without causing overflow, then add 2 and make the
    // result even by clearing the least significant bit. The resulting sequence of capacities
    // is carefully chosen to hit a maximum capacity that is close to the '2^32-1' limit when
    // using 'uint32_t' as indices so that as much as possible of this space can be used.
    uint32_t delta = ((cap >> 1) + (cap >> 3) + 2) & ~1;
    cap += delta;
    if (cap <= prev_cap) throw OutOfMemoryException();
  }

  assert(cap > 0);
  memory = (uint32_t*)xrealloc(memory, sizeof(uint32_t) * cap);
}

// segment tree (fast implementation of priority queue).
void OrderHeap::resize(int newsize) {
  if (cap >= newsize) return;
  // insert elements in such order that tie breaking remains intact
  std::vector<Var> variables;
  while (!empty()) variables.push_back(removeMax());
  tree.clear();
  while (cap < newsize) cap = cap * resize_factor + 1;
  tree.resize(2 * (cap + 1), -1);
  for (Var x : variables) insert(x);
}
void OrderHeap::recalculate() {  // TODO: more efficient implementation
  // insert elements in such order that tie breaking remains intact
  std::vector<Var> variables;
  while (!empty()) variables.push_back(removeMax());
  tree.clear();
  tree.resize(2 * (cap + 1), -1);
  for (Var x : variables) insert(x);
}
void OrderHeap::percolateUp(Var x) {
  for (int at = x + cap + 1; at > 1; at >>= 1) {
    if (tree[at ^ 1] == -1 || activity[x] > activity[tree[at ^ 1]])
      tree[at >> 1] = x;
    else
      break;
  }
}
bool OrderHeap::empty() const { return tree[1] == -1; }
bool OrderHeap::inHeap(Var x) const { return tree[x + cap + 1] != -1; }
void OrderHeap::insert(Var x) {
  assert(x <= cap);
  if (inHeap(x)) return;
  tree[x + cap + 1] = x;
  percolateUp(x);
}
Var OrderHeap::removeMax() {
  Var x = tree[1];
  assert(x != -1);
  tree[x + cap + 1] = -1;
  for (int at = x + cap + 1; at > 1; at >>= 1) {
    if (tree[at ^ 1] != -1 && (tree[at] == -1 || activity[tree[at ^ 1]] > activity[tree[at]]))
      tree[at >> 1] = tree[at ^ 1];
    else
      tree[at >> 1] = tree[at];
  }
  return x;
}

}  // namespace rs
