#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <vector>
#include <iostream>
#include <cstdarg>
#include <cassert>
#include "halfint.hpp"
#include "utils.hpp"

using namespace std;

class idx {
public:
  idx();
  idx(int);
  idx(halfint _min, halfint _max, halfint _delta=2*half);
  halfint min() const;
  halfint max() const;
  halfint delta() const;
  int dim() const;
  int operator()(halfint h) const;
  bool match(halfint h) const;
  friend ostream& operator<<(ostream&, const idx&);
  static const idx idx_0;
  static const idx idx_lor;
  static const idx idx_s1;
  static const idx idx_s10;
  static const idx idx_s1h;
  static const idx idx_s3h;
  static const idx idx_s5h;
private:
  halfint _min;
  halfint _max;
  halfint _delta;
  int _dim;
};

extern const idx& idx_0;     // one element array  
extern const idx& idx_lor;   // lorentz index
extern const idx& idx_s1;    // spin 1             
extern const idx& idx_s1_0;  // spin 1 massless    
extern const idx& idx_s1h;   // spin 1/2           
extern const idx& idx_s3h;   // spin 3/2           
extern const idx& idx_s5h;   // spin 5/2           

const idx& idx_spin(halfint spin);

/**
   Represent an N dimensional array of elements of type T.
*/
template <typename T, int N>
class Array {
public:
  /**
     Constructor. The array is filled with default constructed
     elements. This depends on the implementation of the vector.resize()
     method, therefore there is no initialization for primitive types
     like int or double.
  */
  /*
    This solution has been chosen because it's illegal to pass non-POD types through
    ..., and halfint and idx are non-POD  */
  Array(idx _idx0, idx _idx1=idx_0, idx _idx2=idx_0, idx _idx3=idx_0, idx _idx4=idx_0,
           idx _idx5=idx_0, idx _idx6=idx_0, idx _idx7=idx_0, idx _idx8=idx_0, idx _idx9=idx_0) {
    // cerr << " idx0= " << _idx0 << endl;
    // cerr << " idx1= " << _idx1 << endl;
    // cerr << " idx2= " << _idx2 << endl;
    _idx.resize(10);
    _idx[0] = _idx0;
    _idx[1] = _idx1;
    _idx[2] = _idx2;
    _idx[3] = _idx3;
    _idx[4] = _idx4;
    _idx[5] = _idx5;
    _idx[6] = _idx6;
    _idx[7] = _idx7;
    _idx[8] = _idx8;
    _idx[9] = _idx9;
    nelement = 1;
    for (int k(0); k<N; k++) {
      nelement *= _idx[k].dim();
      //PR(nelement);
    }
    vals.resize(nelement);


    // cerr << "Array::constr. indices: " << endl;
    // for (auto it=_idx.begin(); it!=_idx.end(); it++) {
    //   cerr << *it << endl;
    // }
    // cerr << endl;

  }

  void fill(const T& val) {
    for (int i(0); i<nelement; i++) { vals[i] = val; }
  }

  T& operator()(halfint h0, halfint h1=_0, halfint h2=_0, halfint h3=_0, halfint h4=_0,
                halfint h5=_0, halfint h6=_0, halfint h7=_0, halfint h8=_0, halfint h9=_0) {
    vector<int> i;
    // cerr << "h0 = " << h0 << endl;
    // cerr << "h1 = " << h1 << endl;
    i.resize(10);
    i[0] = _idx[0](h0);
    // cerr << "itt vagyok!" << endl;
    i[1] = _idx[1](h1);
    i[2] = _idx[2](h2);
    i[3] = _idx[3](h3);
    i[4] = _idx[4](h4);
    i[5] = _idx[5](h5);
    i[6] = _idx[6](h6);
    i[7] = _idx[7](h7);
    i[8] = _idx[8](h8);
    i[9] = _idx[9](h9);
    int ind = index(i);
    return vals[ind];
  }

  const T& operator()(halfint h0, halfint h1=_0, halfint h2=_0, halfint h3=_0, halfint h4=_0,
                halfint h5=_0, halfint h6=_0, halfint h7=_0, halfint h8=_0, halfint h9=_0) const {
    vector<int> i(10);
    i[0] = _idx[0](h0);
    i[1] = _idx[1](h1);
    i[2] = _idx[2](h2);
    i[3] = _idx[3](h3);
    i[4] = _idx[4](h4);
    i[5] = _idx[5](h5);
    i[6] = _idx[6](h6);
    i[7] = _idx[7](h7);
    i[8] = _idx[8](h8);
    i[9] = _idx[9](h9);
    int ind = index(i);
    return vals[ind];
  }

private:
  vector<idx> _idx;
  int nelement;
  vector<T> vals;

  // default construction forbidden:
  Array();

  // copy construction forbidden:
  //  Array(const Array<T,N>&);

  int index(const vector<int>& i) const {
    /*
    cerr << "Array::index" << endl;
    cerr << "( ";
    for (int k(0); k<N; k++) { cerr << i[k] << ", "; }
    cerr << ")   of   ( ";
    for (int k(0); k<N; k++) { cerr << _idx[k].dim() << ", "; }
    cerr << ")" << endl;
    //*/
    for (int k(0); k<N; k++) {
      if (!(i[k]<_idx[k].dim())) {
        PR(N);
        for (int ii(0); ii<N; ii++) {
          PL(ii); PL(i[ii]); PR(_idx[ii].dim());
        }
        PL(nelement); PL(k); PL(i[k]); PR(_idx[k].dim()); } 
      assert(i[k]<_idx[k].dim());
    }
    int ind = i[0];
    for (int k(1); k<N; k++) {
      ind *= _idx[k].dim();
      ind += i[k];
    }
    return ind;
  }

};


#endif //ARRAY_HPP
