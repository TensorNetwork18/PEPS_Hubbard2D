#ifndef __SRC_MATH_DIIS_H
#define __SRC_MATH_DIIS_H

#include <type_traits>
#include <utility>
#include <memory>
#include <vector>
#include <list>
#include "matrix.h"

// std::shared_ptr<T> is assumed to be a shared_pointer of some class
// which have daxpy and ddot functions.
// T must have clone() function that returns shared_ptr<T>
// Mat should be either Matrix or ZMatrix. Needs Mat::solve function

template <class T, typename Mat = Matrix>
class DIIS {
  private:
    int ndiis_;

    std::list<std::pair<std::shared_ptr<const T>, std::shared_ptr<const T>>> data_;

    std::shared_ptr<Mat> matrix_;
    std::shared_ptr<Mat> coeff_;

  public:
    DIIS() { }
    DIIS(const int ndiis) : ndiis_(ndiis), matrix_(std::make_shared<Mat>(ndiis+1, ndiis+1)), coeff_(std::make_shared<Mat>(ndiis+1, 1)) { }

    std::shared_ptr<T> extrapolate(const std::pair<std::shared_ptr<const T>, std::shared_ptr<const T>> input) {
      std::shared_ptr<const T> v = input.first;
      std::shared_ptr<const T> e = input.second;
      data_.push_back(input);

      if (data_.size() > ndiis_) {
        data_.pop_front();
        matrix_->copy_block(0, 0, ndiis_-1, ndiis_-1, *matrix_->get_submatrix(1, 1, ndiis_-1, ndiis_-1));
      }
      const size_t cnum = data_.size();
      auto data_iter = data_.begin();

      // left hand side
      for (size_t i = 0; i != cnum - 1; ++i, ++data_iter) {
        matrix_->element(cnum-1, i) = e->dot_product(*(data_iter->second));
        matrix_->element(i, cnum-1) = matrix_->element(cnum-1, i);
      }
      matrix_->element(cnum-1, cnum-1)= e->dot_product(*e);
      for (size_t i = 0; i != cnum; ++i)
        matrix_->element(cnum, i) = matrix_->element(i, cnum) = -1.0;
      matrix_->element(cnum, cnum) = 0.0;

      // right hand side
      for (size_t i = 0; i != cnum; ++i)
        coeff_->element(i,0) = 0.0;
      coeff_->element(cnum,0) = -1.0;

      // solve the linear equation
      coeff_ = coeff_->solve(*matrix_, cnum+1);

      // return a linear combination
      std::shared_ptr<T> out = input.first->clone();
      data_iter = data_.begin();
      for (size_t i = 0; i != cnum; ++i, ++data_iter)
        out->ax_plus_y(coeff_->element(i,0), *(data_iter->first));
      return out;
    }

};
#endif

