#ifndef _PEPS_H_
#define _PEPS_H_

#include <vector>
#include <memory>
#include <string>
#include "../sparse_tensor/basis.h"
#include "../sparse_tensor/qtensor.h"

class PEPS {
  private:
    size_t nrows_;
    size_t ncols_;
    std::vector<Qtensor> data_; //row_wise
  public:
    virtual ~PEPS()=default;
    PEPS()=default;
    PEPS(const size_t, const size_t, const std::vector<Qtensor>&);
    
    PEPS(const PEPS&);
    PEPS(PEPS&&);
    PEPS& operator = (const PEPS&);
    PEPS& operator = (PEPS&&);

    std::shared_ptr<PEPS> copy() const { return std::make_shared<PEPS>(*this); }
    std::shared_ptr<PEPS> conjugate_copy() const;

    size_t length() const { return data_.size(); }
    size_t nrows() const { return nrows_; }
    size_t ncols() const { return ncols_; }
    const std::vector<Qtensor>& data() const { return data_; }
    std::vector<Qtensor>& data() { return data_; }

    const Qtensor& site(const size_t i, const size_t j) const { assert(i < nrows_ && j < ncols_); return data_[i*ncols()+j]; }
    const Qtensor& operator () (const size_t i, const size_t j) const { return site(i,j); }
    Qtensor& site(const size_t i, const size_t j) { assert(i < nrows_ && j < ncols_); return data_[i*ncols()+j]; }
    Qtensor& operator () (const size_t i, const size_t j) { return site(i,j); }
    Qtensor conjugate_site(const size_t, const size_t) const;

    void randomize(const double lbound = -1.0, const double ubound = 1.0);
    void fill(const double);
    void scale(const double);
    Qtensor rdmL(const size_t, const size_t) const;
    Qtensor rdmR(const size_t, const size_t) const;
    Qtensor rdmU(const size_t, const size_t) const;
    Qtensor rdmD(const size_t, const size_t) const;
    Qtensor rdmS(const size_t, const size_t) const;

    friend std::ostream& operator << (std::ostream& os, const PEPS& o);

    void save() const;
    void load(const std::string&);
};
#endif
