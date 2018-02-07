#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <set>
#include "../shared_tensor/tag_generator.h"
#include "basis.h"
#include "qtensor.h"

using namespace std;


Qtensor::Qtensor(const Qnum charge_in) : charge_(charge_in) {
}

Qtensor::Qtensor(const Qnum charge_in, const vector<Basis>& index_in, const bool build) : charge_(charge_in), index_(index_in) {
  if (build)
    init();
}

Qtensor::Qtensor(const vector<Basis>& index_in, const bool build) : Qtensor(Qnum::zero(), index_in, build) {
}

Qtensor::Qtensor(const Qtensor& obj) : charge_(obj.charge_), index_(obj.index_) {
  obj.load();
  for_each(obj.data_.begin(), obj.data_.end(), [&](const map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair){
    data_.emplace(pair.first, pair.second->copy());
  });
  obj.dump();
}

Qtensor::Qtensor(Qtensor&& obj) noexcept : charge_(obj.charge_), index_(move(obj.index_)), data_(move(obj.data_)), disk_(obj.disk_), modified_(obj.modified_) {
  file_.swap(obj.file_);
}

Qtensor& Qtensor::operator = (const Qtensor& obj) {
  if (this != &obj) {
    charge_ = obj.charge_;
    index_ = obj.index_;
    data_.clear();
    obj.load();
    for_each(obj.data_.begin(), obj.data_.end(), [&](const map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair){
      data_.emplace(pair.first, pair.second->copy());
    });
    obj.dump();
  }
  return *this;
}

Qtensor& Qtensor::operator = (Qtensor&& obj) noexcept {
  if (this != &obj) {
    charge_ = obj.charge_;
    index_ = move(obj.index_);
    data_ = move(obj.data_);
    disk_ = obj.disk_;
    modified_ = obj.modified_;
    file_.swap(obj.file_);
  }
  return *this;
}

shared_ptr<Tensor> Qtensor::get_block() const {
  assert(block_size() == 1);
  load();
  auto out = data_.begin()->second->copy();
  dump();
  return out;
}

std::shared_ptr<Tensor> Qtensor::get_block(const std::vector<Qnum>& o) const {
  auto it = find_block(o);
  if (it != data_.end()) {
    load();
    auto out = it->second->copy();
    dump();
    return out;
  }
  return nullptr;
}

bool Qtensor::put_block(const Tensor& tensor) {
  if (block_size() == 1) {
    load();
    assert(tensor.extent() == data_.begin()->second->extent());
    data_.begin()->second = tensor.copy();
    modified_ = true;
    dump();
    return true;
  }
  return false;
}

bool Qtensor::put_block(const vector<Qnum>& indexQ, const Tensor& tensor) {
  auto it = find_block(indexQ);
  if (it != data().end()) {
    load();
    assert(tensor.extent() == it->second->extent());
    it->second = tensor.copy();
    modified_ = true;
    dump();
    return true;
  }
  return false;
}

bool Qtensor::insert_block(const vector<Qnum>& indexQ, const Tensor& tensor) {
#ifndef NDEBUG 
  assert(!has_block(indexQ));
  assert(indexQ.size() == rank());
  Qnum sum_q = charge();
  for (size_t i = 0; i != rank(); ++i) {
    assert(index(i).has_basis(indexQ[i]));
    if (index(i).type() == Basis::IN)
      sum_q += indexQ[i];
    else if (index(i).type() == Basis::OUT)
      sum_q -= indexQ[i];
  } 
  assert(sum_q == Qnum::zero());
#endif
  load();
  auto out = data_.emplace(indexQ, tensor.copy()).second;
  modified_ = true;
  dump();
  return out;
}

void Qtensor::randomize(const double lbound, const double ubound) {
  load();
  for_each(data_.begin(), data_.end(), 
    [&lbound, &ubound](map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair) { 
    pair.second->randomize(lbound, ubound); 
  });
  modified_ = true;
  dump();
  return;
}

void Qtensor::fill(const double value) {
  load();
  for_each(data_.begin(), data_.end(), 
    [&value](map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair) { 
    pair.second->fill(value); 
  });
  modified_ = true;
  dump();
  return;
}

void Qtensor::reverse() {
  charge_ = -charge_;
  for_each(index_.begin(), index_.end(), [](Basis& b){
    b = b.reverse();
  });
}

double Qtensor::max() const {
  load();
  auto e = max_element(data_.begin(), data_.end(), 
    [](const map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair1, const map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair2){ 
    return pair1.second->max() < pair2.second->max(); 
  });
  const double out = e->second->max();
  dump();
  return out;
}

ostream& operator << (ostream& os, const Qtensor& obj) {
  os << "Charge: " << obj.charge() << "\n";
  for (auto&& i : obj.index()) {
  //  os << i;
  //}
    switch(i.type()){
      case Basis::IN:{
        os << "IN  " << " | ";
        break;
      }
      case Basis::OUT:{
        os << "OUT " << " | ";
        break;
      }
      default:
        throw logic_error("Qtensor::ostream failed");
    }
  } os << "\n";

  if (obj.block_size() != 0) {
    obj.load();
    for (auto&& pair : obj.data()) {
      for (auto&& p1 : pair.first) {
        os << p1 << " | ";
      }
      os << "\n" << *pair.second << "\n";
    }
    obj.dump();

  } else {
    os << "Qtensor Block size is 0" << "\n";
  }
  return os;
}

void Qtensor::broadcast(int source) {
  if (mpi__->size() != 1) { 
    const int iproc = mpi__->rank();

#ifndef NDEBUG
    if (iproc != source) {
      assert(index_.size() == 0);
      assert(data_.size() == 0);
    }  else {
      assert(index_.size() != 0);
      assert(data_.size() != 0);
    }
#endif

    // charge_ //
    charge_.broadcast(source);

    // index_ //
    {
      size_t nrank = (iproc == source) ? rank() : 0u;
      mpi__->broadcast(&nrank, 1u, source);
      if (iproc != source)
        index_.resize(nrank);
      for (size_t i = 0; i != nrank; ++i)
        index_[i].broadcast(source);
    }

    // disk_ //
    {  
      int dvalue = (iproc == source) ? static_cast<int>(disk_) : 0;
      mpi__->broadcast(&dvalue, 1u, source);
      if (iproc != source)
        disk_ = static_cast<bool>(dvalue);
    }

    // data_ //
    { 
       const size_t nrank = rank();
       size_t bsize = block_size();
       mpi__->broadcast(&bsize, 1u, source);
       std::unique_ptr<int[]> qnums(new int[bsize*nrank]);
     
       if (iproc == source) {
         size_t icount = 0;
         for(const auto& i : data_) {
           for (size_t j = 0; j != nrank; ++j) {
             qnums[icount++] = i.first[j].u1();
           }
         }
         assert(icount == bsize*nrank);
       }
       mpi__->broadcast(qnums.get(), bsize*nrank, source);

       if (iproc != source) {
         size_t icount = 0;
         for (size_t i = 0; i != bsize; ++i) {
           vector<Qnum> indexQ(nrank);
           vector<size_t> indexD(nrank);
           for (size_t j = 0; j != nrank; ++j) {
             indexQ[j] = Qnum(qnums[icount++], true);
             indexD[j] = index_[j].qdim(indexQ[j]);
           }
           data_.emplace(indexQ, make_shared<Tensor>(indexD));
         }
       }

       broadcast_pack(source);
    }

  }
  return;
}

void Qtensor::broadcast_pack(int source) {
  if (mpi__->size() != 1) {
    load();
    const size_t total_element = (rank() != 0) ? accumulate(data_.begin(), data_.end(), 0, [](size_t value, const map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair) { return value + pair.second->size(); }) : 0;

    auto vec = make_shared<Vec>(total_element);
    { //pack
      size_t position = 0;        
      for (const auto& i : data_) {
        const size_t isize = i.second->size();
        copy_n(i.second->pdata(), isize, vec->element_ptr(position));
        position += isize;
      }   
      assert(position == total_element);    
    }

    mpi__->broadcast(vec->pdata(), vec->size(), source);
    
    { //unpack
      size_t position = 0;
      for (auto&& i : data_) {
        const size_t isize = i.second->size();
        copy_n(vec->element_ptr(position), isize, i.second->pdata());
        position += isize;
      }
      assert(position == total_element);
    }

    if (disk_) {
      if (file_.length() == 0) {
        assert(modified_ == true);
        file_ = tag__->get();
      }
      dump();
    }
  }
  return;
}

void Qtensor::allreduce() {
  if (mpi__->size() != 1) {
    assert(disk_ == false);
    Qtensor tmp(charge(), index());
    tmp += *this;
    tmp.allreduce_pack();
    *this = move(tmp);
  }
  return; 
}

void Qtensor::allreduce_pack() {
  if (mpi__->size() != 1) {
    assert(disk_ == false);
    const size_t total_element = (rank() != 0) ? accumulate(data_.begin(), data_.end(), 0, [](size_t value, const map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair) { return value + pair.second->size(); }) : 0;
    
    auto vec = make_shared<Vec>(total_element);
    { //pack
      size_t position = 0;    
      for (const auto& i : data_) {
        const size_t isize = i.second->size();
        copy_n(i.second->pdata(), isize, vec->element_ptr(position));
        position += isize;
      }
      assert(position == total_element);    
    }

    mpi__->allreduce(vec->pdata(), vec->size());
    
    { //unpack
      size_t position = 0;
      for (auto&& i : data_) {
        const size_t isize = i.second->size();
        copy_n(vec->element_ptr(position), isize, i.second->pdata());
        position += isize;
      }
      assert(position == total_element);
    }
  }
  return; 
}

void Qtensor::load() const {
  if (disk_ && !modified_) {  // I assume if modified_ is false, it has been writtened via flush or dump
    assert(file_.length()!=0);
    read();
  }
  return;
}

void Qtensor::dump() const {
  if (disk_) {
    assert(file_.length() != 0);
    if (modified_) {
      write();
      modified_ = false;
    }

    // free
    for (auto&& pair : data_)
      pair.second = nullptr;
  }
  return;
}

void Qtensor::flush() { //can only be called by from outside
  assert(disk_ == false);
  assert(modified_ == true);
  disk_ = true;
  if (file_.length() == 0)
    file_ = tag__->get();
  dump();
  return;
}

// Private //
void Qtensor::init() { //init all symmetry allowed blocks
  const size_t nrank = rank();

  //a vector of iterator
  vector<map<Qnum, size_t>::const_iterator> it(nrank);
  for (size_t i = 0; i != nrank; ++i)
    it[i] = index()[i].basis().begin();

  //Loop over all possible combination
  while (it[0] != index()[0].basis().end()) {
    //Main work//
    Qnum sum_q = charge();

    for (size_t i = 0; i != nrank; ++i) {
      if (index()[i].type() == Basis::IN) {
        sum_q += it[i]->first;
      } else if (index()[i].type() == Basis::OUT) {
        sum_q -= it[i]->first;
      }   
    }

    if (sum_q == Qnum::zero()) {
      vector<Qnum> indexQ; 
      vector<size_t> indexD; 
      for (auto&& i : it) {
        indexQ.push_back(i->first);
        indexD.push_back(i->second);
      }
      data_.emplace(indexQ, make_shared<Tensor>(indexD));
    }
 
    //Loop//
    ++it[nrank-1];
    for (size_t i = nrank - 1; (it[i] == index()[i].basis().end()) &&  i != 0; --i) {
      it[i] = index()[i].basis().begin();
      ++it[i-1];
    }
  }
  assert(block_size() != 0); 
}

void Qtensor::write() const {
  ofstream output;
  output.open(file_, ios::binary);
  if (!output.good())
    throw runtime_error("Qtensor::write fails " + file_);
#ifndef NDEBUG
  for (auto&& pair : data_)
    assert(pair.second != nullptr);
#endif
  for (const auto& pair : data_) {
    const Tensor& ref = *pair.second;
    //Rank
    const size_t rank = ref.rank();
    binary_write(output, &rank, 1); 

    //Extent
    const std::vector<size_t>& extent = ref.extent();
    binary_write(output, &extent[0], rank);

    //Element
    const double* pdata = ref.pdata();
    binary_write(output, pdata, ref.size());
  }
  output.close();
  return;
}

void Qtensor::read() const {
   ifstream input;
   input.open(file_, ios::binary);
   if (!input.good())
     throw runtime_error("Qtensor::read fails " + file_);
#ifndef NDEBUG
  for (auto&& pair : data_)
    assert(pair.second == nullptr);
#endif
  for (auto&& pair : data_) {
    size_t rank;
    binary_read(input, &rank, 1); 
    
    std::vector<size_t> extent(rank);
    binary_read(input, &extent[0], rank);
    
    pair.second = make_shared<Tensor>(extent);
    double* pdata = pair.second->pdata();
    binary_read(input, pdata, pair.second->size()); 
  }
  input.close();

  return;
}

bool Qtensor::sgn_permute(const vector<size_t>& fidx){
  // Ref: https://stackoverflow.com/questions/20702782/efficiently-determine-the-parity-of-a-permutation
  // Odd  -> True
  // Even -> False
  const size_t ns = fidx.size();
  if (ns > 1) {
    vector<size_t> idx(ns);
    iota(idx.begin(), idx.end(), 0); 
    sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return fidx[a] < fidx[b]; });   

    int count = 0;
    long visited = 0;
    for (size_t i = 0; i != ns; ++i) {
      long mask = (1L<<i);
      if ((visited & mask) != 0)
        continue;
      visited |= mask;
      for (size_t j = idx[i]; (visited & (1L<<j)) == 0; j = idx[j]) {
        visited |= (1L<<j);
        ++count;
      }   
    }   
    return (count%2 == 1); 
  }
  return false;
}
