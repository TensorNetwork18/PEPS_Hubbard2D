#ifndef _TENSOR_BASE_H_
#define _TENSOR_BASE_H_

#include <cassert>
#include <cstddef>
#include <type_traits>
#include <memory>
#include <vector>
#include <algorithm>
#include <random>
#include <initializer_list>
#include <iostream>
#include <iomanip>



template<typename DataType>
class Tensor_base {
  protected:
    std::vector<size_t> extent_;
    std::unique_ptr<DataType []> data_;
  public:
    virtual ~Tensor_base() {
      if (!extent_.empty())
        extent_.clear();
      if (pdata()) 
        data_.reset();
    }

  //Constructor
    Tensor_base()=default;

    template<typename arg0, typename... args, typename = typename std::enable_if<std::is_convertible<arg0, size_t>::value, size_t>::type>
    Tensor_base(const arg0& extent0, const args&... extents) {
      extent_ = {extent0, extents...};
      data_.reset(new DataType[size()]);
      fill(0);
    }

    Tensor_base(const std::vector<size_t>& extent) {
      extent_ = extent;
      data_.reset(new DataType[size()]);
      fill(0);
    }

    Tensor_base(const Tensor_base& o) {
     extent_ = o.extent_;
     data_.reset(new double[o.size()]);
     std::copy_n(o.pdata(), size(), pdata());
    }

    Tensor_base(Tensor_base&& o) {
      extent_ = move(o.extent_);
      data_ = move(o.data_);
    }

  //Assignment Operators
    Tensor_base& operator = (const Tensor_base& o) {
      extent_ = o.extent_;
      data_.reset(new double[o.size()]);
      std::copy_n(o.pdata(), size(), pdata());
      return *this;    
    }
    
    Tensor_base& operator = (Tensor_base&& o) {
      if (this != &o) {
        extent_ = move(o.extent_);
        data_ = move(o.data_);   
      }
      return *this;
    }

  //Basic Properties  
    size_t rank() const {
      return extent_.size(); 
    }

    size_t extent(const size_t i) const {
      if (rank() != 0 && i < rank()) {
        return extent_[i];
      } else {
        return 0;
      }
    }       

   const std::vector<size_t>& extent() const {
      return extent_;
    }

    size_t size() const {
      if (rank() != 0) {
        size_t size = 1;
        for_each(extent().begin(), extent().end(), [&](size_t i) {size *= i;});
        return size;
      } else { 
        return 0;
      } 
    }

  //Element Accessor
    DataType* pdata() const { 
      return data_.get(); 
    }

    void fill(const DataType a) { 
      std::fill_n(pdata(), size(), a); 
    }

    void set_elements(const size_t na, const double* a) {
      assert(na == size());
      std::copy_n(a, na, pdata());
    }

    void randomize(const double lower_bound = -1.0, const double upper_bound = 1.0) {
      assert(!(lower_bound>upper_bound));
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<double> dis(lower_bound, upper_bound);
      std::for_each(pdata(), pdata()+size(), [&gen, &dis](double& val) { val = dis(gen); });
    }

    template<typename... args>
    DataType* element_ptr(const args&... indx) {
      assert(sizeof...(args) <= rank());
      return pdata() + offset(indx...);
    }
    
    template<typename... args>
    DataType& element(const args&... indx) { 
      return *element_ptr(indx...); 
    }   
  
    template<typename... args>
    DataType& operator () (const args&... indx) { 
      return element(indx...); 
    }

    template<typename... args>
    const DataType* element_ptr(const args&... indx) const {
      assert(sizeof...(args) <= rank());
      return pdata() + offset(indx...);
    } 
    
    template<typename... args>
    const DataType& element(const args&... indx) const { 
      return *element_ptr(indx...); 
    }   
  
    template<typename... args>
    const DataType& operator() (const args&... indx) const { 
      return element(indx...); 
    }

    template<typename... args>
    void reshape(const args&... extent) {
      assert(sizeof...(extent) != 0);
      std::vector<size_t> tmp = {extent...};
      size_t n = 1;
      for_each(tmp.begin(), tmp.end(), [&](size_t i) {n *= i;} );
      assert(n == size());
      extent_ = tmp;
    }

    void reshape(const std::vector<size_t>& extent) {
      assert(extent.empty() == false);
      size_t n = 1;
      for_each(extent.begin(), extent.end(), [&](size_t i) {n *= i;} );
      assert(n == size());
      extent_ = extent;
    }
  
  private:
    // map n_tuple of coordinates into one-dimensional memory reference
    template<typename arg0, typename... args>
    size_t offset(const arg0& first, const args&... rest) const {
      const size_t i = rank() - sizeof...(rest) - 1;
      return first + extent(i) * offset(rest...);
    }
    
    template<typename arg0>
    size_t offset(const arg0& first) const {
      return first;
    }   
};
#endif /* _TENSOR_BASE_H  */
