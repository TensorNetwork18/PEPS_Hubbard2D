#ifndef _VEC_H_
#define _VEC_H_

#include "tensor_base.h"
#include "mkl_interface.h"

class Vec: public Tensor_base<double> {
  public:
    virtual ~Vec()=default;
    
    Vec() { }
    Vec(const size_t n) : Tensor_base<double>(n) { }
    
    Vec(const Vec& o) : Tensor_base<double>(o) { }
    Vec(Vec&& o) : Tensor_base<double>(std::move(o)) { }
    Vec& operator = (const Vec& o) { Tensor_base<double>::operator=(o); return *this; }
    Vec& operator = (Vec&& o) { Tensor_base<double>::operator=(o); return *this; }
       
    size_t size() const { return extent(0); }

    std::shared_ptr<Vec> clone() const { return make_shared<Vec>(size()); }
    std::shared_ptr<Vec> copy() const { return make_shared<Vec>(*this); }

    void ax_plus_y (const double a, const Vec& o) { assert(size() == o.size()); daxpy(a, o, *this); }
    Vec& operator += (const Vec& o) { this->ax_plus_y(1.0, o); return *this; }
    Vec operator + (const Vec& o) const { Vec out(*this); out += o; return out; }
    Vec& operator -= (const Vec& o) { this->ax_plus_y(-1.0, o); return *this; }
    Vec operator - (const Vec& o) const { Vec out(*this); out -= o; return out; }
    
    double dot_product(const Vec& o) const { return ddot(*this, o); }
    void scale(const double a) { dscal(a, *this); }

    friend std::ostream& operator << (std::ostream& os, const Vec&);
};

#endif
