#ifndef QTENSOR_H_
#define QTENSOR_H_

#include <iostream>
#include <cassert>
#include <utility>
#include <iterator>
#include <memory>
#include <vector>
#include <string>
#include <map>
#include "../mpi/mpi_interface.h"
#include "../math/tensor.h"
#include "../math/matrix.h"
#include "../math/vec.h"
#include "../shared_tensor/shared_tensor.h"
#include "basis.h"

class Qtensor {
  //TODO decouple direction info from Basis
  private:
    Qnum charge_;
    std::vector<Basis> index_;
    mutable std::map<std::vector<Qnum>, std::shared_ptr<Tensor>> data_;

    bool disk_ = false;
    mutable bool modified_ = true;
    std::string file_;

  public:
    //virtual ~Qtensor()=default;
    virtual ~Qtensor() {
      if (file_.length() != 0) {
        struct stat buf;
        if (stat(file_.c_str(), &buf) == 0) {
          remove( file_.c_str() );
        }
      }
      index_.clear();
      data_.clear();
    }

    Qtensor(const Qnum = Qnum::zero());
    Qtensor(const Qnum, const std::vector<Basis>&, const bool = true);
    Qtensor(const std::vector<Basis>&, const bool = true);

    Qtensor(const Qtensor&);
    Qtensor(Qtensor&&) noexcept;
    Qtensor& operator = (const Qtensor&);
    Qtensor& operator = (Qtensor&&) noexcept;

    std::shared_ptr<Qtensor> copy() const { return std::make_shared<Qtensor>(*this); }

    bool disk() const { return disk_; }
    string file() const { return file_; }
    void build() { if(data_.empty()) init(); }
    size_t rank() const { return index_.size(); }
    size_t block_size() const { return data_.size(); }
    size_t total_element() const {
      load();
      const size_t out = (rank() != 0) ? std::accumulate(data_.begin(), data_.end(), 0, [](size_t value, const std::map<vector<Qnum>, std::shared_ptr<Tensor>>::value_type& pair) { return value + pair.second->size(); }) : 0;
      dump();
      return out;
    }
    size_t total_size() const { return (rank() != 0) ? std::accumulate(index().begin(), index().end(), 1, [](size_t value, const Basis& e){ return value * e.total_dim(); }) : 0; }
    size_t total_dim(const size_t o) const { return (rank()!=0) ? index(o).total_dim() : 0; }
    size_t qdim(const size_t o1, const Qnum& o2) const { return (o1 < rank()) ? index(o1).qdim(o2) : 0; }
    Basis index(const size_t o) const { assert(o < rank()); return index_[o]; }

    const Qnum charge() const { return charge_; }
    Qnum charge() { return charge_; }
    bool parity() const { return (charge_.parity() == Qnum::ODD); } //ODD->True EVEN->False
    const std::vector<Basis>& index() const { return index_; } 
    std::vector<Basis>& index() { return index_; }
    const std::map<std::vector<Qnum>, std::shared_ptr<Tensor>>& data() const { return data_; }
    std::map<std::vector<Qnum>, std::shared_ptr<Tensor>>& data() { return data_; }

    bool has_block(const std::vector<Qnum>& o) const { return data_.count(o) == 1; }
    std::map<std::vector<Qnum>, std::shared_ptr<Tensor>>::const_iterator find_block(const std::vector<Qnum>& o) const { return data_.find(o); }
    std::map<std::vector<Qnum>, std::shared_ptr<Tensor>>::iterator find_block(const std::vector<Qnum>& o) { return data_.find(o); }
    std::shared_ptr<Tensor> get_block() const;
    std::shared_ptr<Tensor> get_block(const std::vector<Qnum>& o) const;
    bool put_block(const Tensor&);
    bool put_block(const std::vector<Qnum>&, const Tensor&);
    bool insert_block(const std::vector<Qnum>&, const Tensor&);

    void randomize(const double lbound = -1.0, const double ubound = 1.0);
    void fill(const double a);

    void reverse();
    double max() const;
    friend std::ostream& operator << (std::ostream& os, const Qtensor& o);
    void broadcast(int source);
    void broadcast_pack(int source);
    void allreduce();
    void allreduce_pack();
    void load() const;
    void dump() const;
    void flush();
  private:
    void init();
    void write() const;
    void read() const;
    static bool sgn_permute(const std::vector<size_t>& fidx);

  public:
    bool operator == (const Qtensor& o) const;
    Qtensor& operator += (const Qtensor& o);
    Qtensor operator + (const Qtensor& o) const { Qtensor out(*this); out += o; return out; }
    Qtensor& operator -= (const Qtensor& o);
    Qtensor operator - (const Qtensor& o) const { Qtensor out(*this); out -= o; return out; }
    //Qtensor eigs(double&, const Qtensor&) const;
    std::vector<Qtensor> geigs(std::vector<double>&, const Qtensor&, const size_t = 1u) const;
    std::vector<Qtensor> split(const int) const;
    std::vector<Qtensor> sqrt2(std::vector<double>&) const;
    void scale(const double);
    Qtensor permute(const std::vector<size_t>&, const bool sgn = true) const;
    void contract(const double, const Qtensor&, const std::vector<size_t>&, const Qtensor&, const std::vector<size_t>&, const bool sgn = true);
    std::vector<Qtensor> svd(const size_t indx, const Basis::Type, const char*) const;
};
#endif
