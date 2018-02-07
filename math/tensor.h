#ifndef _TENSOR_H_
#define _TENSOR_H_

#include "tensor_base.h"
#include "mkl_interface.h"

class Tensor : public Tensor_base<double> {
  public:
    virtual ~Tensor()=default;

    Tensor() : Tensor_base<double>() { }
    template<typename... args>
    Tensor(const args&... extent) : Tensor_base<double>(extent...) { }

    Tensor(const Tensor& o) : Tensor_base<double>(o) { }
    Tensor(Tensor&& o) : Tensor_base<double>(std::move(o)) { }
    Tensor& operator = (const Tensor& o) { Tensor_base<double>::operator=(o); return *this; }
    Tensor& operator = (Tensor&& o) { Tensor_base<double>::operator=(o); return *this; }
 
    std::shared_ptr<Tensor> clone() const { return std::make_shared<Tensor>(this->extent()); }
    std::shared_ptr<Tensor> copy() const { return std::make_shared<Tensor>(*this); }

    void ax_plus_y (const double a, const Tensor& o) { assert(extent_ == o.extent_); daxpy(a, o, *this); }
    Tensor& operator += (const Tensor& o) { this->ax_plus_y(1.0, o); return *this; }
    Tensor operator + (const Tensor& o) const { Tensor out(*this); out += o; return out; }
    Tensor& operator -= (const Tensor& o) { this->ax_plus_y(-1.0, o); return *this; }
    Tensor operator - (const Tensor& o) const { Tensor out(*this); out -= o; return out; }
    bool is_equal (const Tensor& o, const double thresh = 1.0e-11) const;
   
    std::vector<std::shared_ptr<Tensor>> svd(const size_t) const;
    std::shared_ptr<Tensor> permute(const std::vector<size_t>&) const;
    void contract(const double, const Tensor&, const std::vector<size_t>&, const Tensor&, const std::vector<size_t>&);
    void scale(const double);
    double max() const;
    double max_abs() const;
    double sum_diag() const;
    friend std::ostream& operator << (std::ostream& os, const Tensor&);
};
#endif /* _TENSOR_H  */
