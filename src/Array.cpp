#include "Array.hpp"
#include <cstdlib>

idx::idx() : _min(0), _max(0), _delta(0), _dim(1) {}
idx::idx(int _dim) : _min(0), _max(_dim-1), _delta(1), _dim(_dim) {}
idx::idx(halfint _min, halfint _max, halfint _delta) : _min(_min), _max(_max), _delta(_delta) {
  // cerr << "idx::constr - _delta = " << _delta << endl;
  // PR(_min); PR(_max);
  assert(twice(_max-_min)%twice(_delta) == 0);
  _dim = twice(_max-_min)/twice(_delta) + 1;
  /*
  cerr << "idx::constr: ";
  PL(_min); PL(_max); PL(_delta); PR(_dim);
  PL(half); PR(3*half);
  //*/
}
halfint idx::min() const { return _min; }
halfint idx::max() const { return _max; }
halfint idx::delta() const { return _delta; }
int idx::dim() const { return _dim; }
int idx::operator()(halfint h) const {
  // cerr << "idx::operator()" << endl;
  // cerr << "_min = " << _min << endl;
  // cerr << "_max = " << _max << endl;
  // cerr << "_delta = " << _delta << endl;
  // cerr << "_dim = " << _dim << endl;
  // cerr << "h = " << h << endl;
  if (!(match(h))) {
    cerr << "idx::operator()(halfint h) const" << endl;
    PL(h); PL(_min); PL(_max); PL(_delta); PR(_dim);
  }
  assert(match(h));
  return twice(h-_min)/twice(_delta);
}

bool idx::match(halfint h) const {
  return ((_min<=h) and (h<=_max) and (twice(h-_min)%twice(_delta) == 0));
}

ostream& operator<<(ostream& out, const idx& IDX) {
  out << "idx object:" << endl;
  out << "  _min =   " << IDX._min << endl;
  out << "  _max =   " << IDX._max << endl;
  out << "  _delta = " << IDX._delta << endl;
  out << "  _dim =   " << IDX._dim << endl;
  return out;
}

const idx idx::idx_0(_0,_0,_1);
const idx idx::idx_lor(_0,_3,_1);
const idx idx::idx_s1(-_1,_1,_1);
const idx idx::idx_s10(-_1,_1,_2);
const idx idx::idx_s1h(-half,half,_1);
const idx idx::idx_s3h(-3*half,3*half,_1);
const idx idx::idx_s5h(-5*half,5*half,_1);

const idx& idx_0(idx::idx_0);
const idx& idx_lor(idx::idx_lor);
const idx& idx_s1(idx::idx_s1);
const idx& idx_s10(idx::idx_s10);
const idx& idx_s1h(idx::idx_s1h);
const idx& idx_s3h(idx::idx_s3h);
const idx& idx_s5h(idx::idx_s5h);



const idx& idx_spin(halfint spin) {
  if (spin==half) {
    return idx_s1h;
  } else if (spin==3*half) {
    return idx_s3h;
  } else if (spin==5*half) {
    return idx_s5h;
  }
  cerr << "Unimplemented spin (" << spin << ") in idx_spin(halfint spin) - only intended for fermions" << endl;
  exit(0);
  // never get here:
  return idx_0;
}

