#ifndef QINDEX_H_
#define QINDEX_H_

#include <iostream>
#include <cassert>
#include <bitset>
#include <vector>
#include <array>
#include "qnum.h"

template<size_t B, size_t N = 10u> /*B = 10 <-> [-500,500]; B = 8 <-> [-100, 100]*/
class Qindex {
  private:
    size_t rank_ = 0;
    std::bitset<B*N> tag_;
  public:
    virtual ~Qindex()=default;
    Qindex()=default;
    Qindex(std::vector<Qnum>::const_iterator first, const size_t range) {
      assert(tag_.none());
      insert(first, range);
    }
    Qindex(const std::vector<Qnum>& o) : Qindex(o.begin(), o.size()) { }
    Qindex(const Qindex& o) : rank_(o.rank_), tag_(o.tag_) {  }
    Qindex& operator = (const Qindex<B>& o) { rank_ = o.rank_; tag_ = o.tag_; return *this; }
    bool operator == (const Qindex<B>& o) const { assert(rank_ == o.rank_); return (tag_ == o.tag_); }
    bool operator != (const Qindex<B>& o) const { return !(*this==o); }
    bool operator < (const Qindex<B>& o) const { 
     assert(rank_ == o.rank_);
     for (size_t i = rank_*B; i-- > 0;) {
       if (tag_[i]^o.tag_[i]) return o.tag_[i];
     }
     return false;
     //return tag_.to_string() < o.tag_.to_string();
     //typedef const array<long, (B*N/64)> AsArray;
     //return *const_cast<AsArray*>(reinterpret_cast<AsArray*>(&tag_)) < *const_cast<AsArray*>(reinterpret_cast<AsArray*>(&o.tag_));

    }

    size_t rank() const { return rank_; }
    std::bitset<B*N> tag() const { return tag_; }


    void insert(std::vector<Qnum>::const_iterator first, const size_t range) {
      rank_ += range;
      assert(rank_ <= N);
      for (size_t i = 0; i != range; ++i) { 
        tag_ <<= B;
        std::bitset<B> tmp(first->u1());
        for (size_t j = 0; j != B; ++j) {
          tag_[j] = tmp[j];
        }
        ++first;
      }
      return;
    }

    Qindex slice(const size_t position, const size_t range) const {
      assert(position + range <= rank_);
      Qindex<B> out;
      out.rank_ = range;
      size_t shift = (rank_ - position) * B;
      for (size_t i = range*B; i-- > 0;)
        out.tag_[i] = tag_[--shift];
      return out;
    }

    Qnum at(const size_t position) const {
      std::bitset<B> tmp;
      size_t shift = (rank_ - position) * B;
      for (size_t i = 0; i != B; ++i) {
         tmp[B-i-1] = tag_[--shift];
      }
      struct {int x:B;} s;  //signed conversion from https://graphics.stanford.edu/~seander/bithacks.html#FixedSignExtend
      s.x = tmp.to_ulong();
      return Qnum(s.x, true); //assume everything is fermion
    }

    std::vector<Qnum> convert() const {
      std::vector<Qnum> out(rank());
      for(size_t i=0; i != rank_; ++i)
        out[i] = at(i);
      return out;
    }
};

#endif
