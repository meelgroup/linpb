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

#include "ConstrSimple.hpp"
#include "ConstrExp.hpp"

namespace rs {

template <typename CF, typename DG>
CeSuper ConstrSimple<CF, DG>::toExpanded(ConstrExpPools& cePools) const {
  // TODO: make this the minimal bitwidth expanded constraint?
  CePtr<ConstrExp<CF, DG>> ce = cePools.take<CF, DG>();
  ce->addRhs(rhs);
  for (const Term<CF>& t : terms) {
    ce->addLhs(t.c, t.l);
  }
  ce->orig = orig;
  if (ce->plogger) {
    ce->proofBuffer.str(std::string());
    ce->proofBuffer << proofLine;
  }
  return ce;
}

template <typename CF, typename DG>
void ConstrSimple<CF, DG>::toNormalFormLit() {
  for (Term<CF>& t : terms) {
    if (t.c < 0) {
      rhs -= t.c;
      t.c = -t.c;
      t.l = -t.l;
    }
  }
}

template <typename CF, typename DG>
void ConstrSimple<CF, DG>::toNormalFormVar() {
  for (Term<CF>& t : terms) {
    if (t.l < 0) {
      rhs -= t.c;
      t.c = -t.c;
      t.l = -t.l;
    }
  }
}

template struct ConstrSimple<int, long long>;
template struct ConstrSimple<long long, int128>;
template struct ConstrSimple<int128, int128>;
template struct ConstrSimple<int128, int256>;
template struct ConstrSimple<bigint, bigint>;

}  // namespace rs
