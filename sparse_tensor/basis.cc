#include <iostream>
#include <cassert>
#include <map>
#include <iterator>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include "../mpi/mpi_interface.h"
#include "basis.h"

using namespace std;


Basis::Basis(const Type type) : type_(type) {
}

Basis::Basis(const Type type, const size_t ldim) : type_(type) {
  if (!basis_.empty()) 
    basis_.clear();
  basis_.emplace(Qnum::zero(), ldim);
}

Basis::Basis(const Type type, const vector<Qnum>& qnums) : type_(type) {
  init(qnums); 
}

Basis::Basis(const Basis& o) : type_(o.type()), basis_(o.basis()) {
}

Basis::Basis(Basis&& o) : type_(o.type()) {
  basis_ = move(o.basis_);
}

Basis& Basis::operator = (const Basis& o) {
  if (this != &o) {
    type_ = o.type();
    basis_ = o.basis();
  }
  return *this;
}

Basis& Basis::operator = (Basis&& o) {
  if (this != &o) {
    type_ = o.type();
    basis_ = move(o.basis_);
  }
  return *this;
}

bool Basis::operator == (const Basis& o) const { 
  return (type() == o.type()) && (basis().size() == o.basis().size()) && equal(basis().begin(), basis().end(), o.basis().begin()); 
}

void Basis::assign(Type type, const size_t bdim) {
  type_ = type;
  if (!basis_.empty()) 
    basis_.clear();
  basis_.emplace(Qnum::zero(), bdim);
}

void Basis::assign(Type type, const vector<Qnum>& qnums) {
  type_ = type;
  init(qnums);
}

Basis Basis::reverse() const {
  Basis out(*this);
  switch(out.type()){
    case IN:{
      out.type_ = OUT;
      break;
    }   
    case OUT:{
      out.type_ = IN;
      break;
    }
    default:
      throw logic_error("Basis::reverse failed");
  }
  return out;
}

Basis Basis::inverse() const {
  Basis out;
  out.type_ = type();
  for_each(basis_.begin(), basis_.end(),[&](const map<Qnum, size_t>::value_type& e){
    out.basis_.emplace(-e.first, e.second);
  });
  return out;
}

Basis Basis::fuse(const Basis& o) const {
  assert(type() == o.type());
  Basis out;
  out.type_ = type();
  for (const auto& i : basis()){
    for (const auto& j : o.basis()){
      Qnum qn = i.first + j.first;
      size_t qdim = i.second * j.second;
      map<Qnum, size_t>::iterator it = out.find(qn);
      if (it != out.basis().end()){
        assert(out.has_basis(it->first));
        (it->second) += qdim;
      } else {
        out.basis_.emplace(qn, qdim);
      }
    }
  }
  return out;
}

ostream& operator << (ostream& os, const Basis& o) {
  switch(o.type()){
    case Basis::IN:{
      os << "IN: " << "\n";
      break;
    }
    case Basis::OUT:{
      os << "OUT: " << "\n";
      break;
    }
    default:
      throw logic_error("Basis::ostream failed");
  }
  for (auto&& i : o.basis()) {
    os << i.first << " ,Dim = " << i.second << "\n";
  }
  return os;
}

void Basis::broadcast(int source) {
  if (mpi__->size() != 1) {
    const int iproc = mpi__->world_rank();
#ifndef NDEBUG 
    if (iproc != source)
      assert(basis_.size() == 0);
    else
      assert(basis_.size() != 0);
#endif

  // Type //
    {
      int t = static_cast<int>(type_);
      mpi__->broadcast(&t, 1u, source);
      if (iproc != source)
        type_ = static_cast<Type>(t);
    }

    // Map //
    {
      size_t msize = size();
      mpi__->broadcast(&msize, 1u, source);
    
      std::unique_ptr<int[]> qnums(new int[msize]);   
      std::unique_ptr<size_t[]> qdims(new size_t[msize]);
    
      if (iproc == source) {
        size_t icount = 0;
        for(const auto& i : basis_) {
          qnums[icount] = i.first.u1();
          qdims[icount] = i.second;
          ++icount; 
        }
      }

      mpi__->broadcast(qnums.get(), msize, source);
      mpi__->broadcast(qdims.get(), msize, source);
    
      if (iproc != source) {
        for (size_t icount = 0; icount != msize; ++ icount) {
          basis_.emplace(Qnum(qnums[icount], true), qdims[icount]);
        }
      }
    }
  }
  return;
}

void Basis::init(const vector<Qnum>& o) {
  assert(o.size() != 0);
  if (!basis_.empty())
     basis_.clear();
  for (auto&& i : o) {
    map<Qnum, size_t>::iterator it = find(i);
    if (it !=  basis().end()){
      assert(has_basis(it->first));
      ++(it->second);
    } else {
      basis_.emplace(i, 1);
    }
  }
}
