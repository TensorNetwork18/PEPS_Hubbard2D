#include <fstream>
#include <array>
#include "../binary_io.h"
#include "../mpi/mpi_interface.h"
#include "peps.h"
#include "etensor.h"

using namespace std;


PEPS::PEPS(const size_t nr, const size_t nc, const vector<Qtensor>& o) : nrows_(nr), ncols_(nc) {
  data_.reserve(o.size());
  for_each(o.begin(), o.end(), [&](const Qtensor& e){
    assert(e.block_size() != 0);
    assert(e.rank() == 5);
    assert(e.index(0).type() == Basis::OUT); //l
    assert(e.index(1).type() == Basis::OUT); //u
    assert(e.index(2).type() == Basis::OUT); //s
    assert(e.index(3).type() == Basis::IN);  //d
    assert(e.index(4).type() == Basis::IN);  //r
    data_.push_back(e); 
  }); 
}

PEPS::PEPS(const PEPS& o) : nrows_(o.nrows_), ncols_(o.ncols_), data_(o.data_) {
}

PEPS::PEPS(PEPS&& o) : nrows_(o.nrows_), ncols_(o.ncols_), data_(move(o.data_)) {
}


PEPS& PEPS::operator = (const PEPS& o) {
  if (this != &o) {
    data_ = o.data_;
  }
  return *this;
}

PEPS& PEPS::operator = (PEPS&& o) {
  if (this != &o) {
    data_ = move(o.data_);
  }
  return *this;
}

shared_ptr<PEPS> PEPS::conjugate_copy() const{
  auto out = make_shared<PEPS>(*this);
  for (size_t i = 0; i != nrows(); ++i) {
    for (size_t j = 0; j != ncols(); ++j) {
      out->site(i,j) = out->conjugate_site(i,j);
    }
  }
  return out;
}

Qtensor PEPS::conjugate_site(const size_t i, const size_t j) const {
  assert(i < nrows_ && j < ncols_);
  auto out = site(i,j).permute({4,3,2,1,0},false);
  for (auto&& k : out.index()) {
    k = k.reverse();
  }
  return out;
}

void PEPS::randomize(const double lbound, const double ubound) {
  if (mpi__->rank() == 0) {
    for (auto&& i : data_)
      i.randomize(lbound, ubound);
  }

  for (auto&& i : data_)
    i.broadcast_pack(0);
}

void PEPS::fill(const double a) {
  for (auto&& i : data_)
    i.fill(a);
}

void PEPS::scale(const double s) {
  for (auto&& i : data_)
    i.scale(s);
}

Qtensor PEPS::rdmL(const size_t i, const size_t j) const {
  const Qtensor& conjugate = conjugate_site(i,j);
  Qtensor out;
  out.contract(1.0, conjugate, {3,2,1,0}, site(i,j), {1,2,3,4}, false);
  return out;
}

Qtensor PEPS::rdmR(const size_t i, const size_t j) const {
  const Qtensor& conjugate = conjugate_site(i,j);
  Qtensor out;
  out.contract(1.0, conjugate, {4,3,2,1}, site(i,j), {0,1,2,3}, false);
  return out;
}

Qtensor PEPS::rdmU(const size_t i, const size_t j) const {
  const Qtensor& conjugate = conjugate_site(i,j);
  Qtensor out;
  out.contract(1.0, conjugate, {4,2,1,0}, site(i,j), {0,2,3,4}, false);
  return out;
}

Qtensor PEPS::rdmD(const size_t i, const size_t j) const {
  const Qtensor& conjugate = conjugate_site(i,j);
  Qtensor out;
  out.contract(1.0, conjugate, {4,3,2,0}, site(i,j), {0,1,2,4}, false);
  return out;
}

Qtensor PEPS::rdmS(const size_t i, const size_t j) const {
  const Qtensor& conjugate = conjugate_site(i,j);
  Qtensor out;
  out.contract(1.0, conjugate, {4,3,1,0}, site(i,j), {0,1,3,4}, false);
  return out;
}

ostream& operator << (ostream& os, const PEPS& o) {
  if (o.length() != 0) {
    assert(o.ncols() != 0 && o.nrows() != 0);
    size_t k = 0;
    for (auto&& i : o.data()) {
      os << "*** PEPS SITE(" << k/o.ncols() << "," << k%o.ncols() << ") ***" << "\n" << i;
      ++k;
    }
  } else{
    os << "PEPS length is 0" << "\n";
  }
  return os;
}

void PEPS::save() const {
  if (mpi__->rank() == 0) {
    for (size_t i = 0; i != nrows(); ++i) {
      for (size_t j = 0; j != ncols(); ++j) {
        ofstream output;
        output.open("wfn." + to_string(i) + "_" + to_string(j), ios::binary);
        //The Number of Blocks
        const size_t nblock = site(i,j).block_size();
        binary_write(output, &nblock, 1);
        for (const auto& pair : site(i,j).data()) {
          //Key
          const auto& indexQ = pair.first;
          array<int,5> qnums;
          for (size_t k = 0; k != 5; ++k)
            qnums[k] = indexQ[k].u1();
          binary_write(output, &qnums[0], 5u);

          //Value
          const Tensor& tensor = *pair.second;
          const size_t nsize = tensor.size();
          binary_write(output, &nsize, 1u);

          const size_t d0 = tensor.extent(0);
          const size_t d1 = tensor.extent(1);
          const size_t d2 = tensor.extent(2);
          const size_t d3 = tensor.extent(3);
          const size_t d4 = tensor.extent(4);
          const double* pdata = tensor.pdata();
          size_t icount = 0;
          for (size_t r4 = 0; r4 != d4; ++r4) {
            for (size_t r3 = 0; r3 != d3; ++r3) {
              for (size_t r2 = 0; r2 != d2; ++r2) {
                for (size_t r1 = 0; r1 != d1; ++r1) { 
                  for (size_t r0 = 0; r0 != d0; ++r0, ++icount) {
                    binary_write(output, &r0, 1u);
                    binary_write(output, &r1, 1u);
                    binary_write(output, &r2, 1u);
                    binary_write(output, &r3, 1u);
                    binary_write(output, &r4, 1u);
                    binary_write(output, pdata+icount, 1u);
                  }
                }
              }
            }
          }
        }
        output.close(); 
      }
    }
  }
  return;
}

void PEPS::load(const string& filename) {
  if (mpi__->rank() == 0) {
    for (size_t i = 0; i != nrows(); ++i) {
      for (size_t j = 0; j != ncols(); ++j) {
        ifstream input;
        input.open(filename + "." + to_string(i) + "_" + to_string(j), ios::binary);
        if (!input.good())
          throw runtime_error("PEPS::load could not open the load file " + filename + "." + to_string(i) + "_" + to_string(j));

        //Get number of Blocks
        size_t nblock;
        binary_read(input, &nblock, 1u);

        for (size_t iblock = 0; iblock != nblock; ++iblock) {        
          //Get Key
          vector<Qnum> indexQ(5);
          for (size_t k = 0; k != 5; ++k) {
            int q;
            binary_read(input, &q, 1u);
            Qnum::Parity p = (q%2 == 0) ? Qnum::EVEN : Qnum::ODD; 
            indexQ[k].assign(q, p);
          }

          //Get Value
          shared_ptr<Tensor> ptr = site(i,j).get_block(indexQ);
          if (ptr) {
            size_t nsize;
            binary_read(input, &nsize, 1u);

            for (size_t isize = 0; isize != nsize; ++ isize) {
              array<size_t, 5> qindex;
              binary_read(input, &qindex[0], 5u);
              double e;
              binary_read(input, &e, 1u);
              if (qindex[0] < ptr->extent(0) && qindex[1] < ptr->extent(1) && qindex[2] < ptr->extent(2) && qindex[3] < ptr->extent(3) && qindex[4] < ptr->extent(4))
                ptr->element(qindex[0], qindex[1], qindex[2], qindex[3], qindex[4]) = e;
            }
            site(i,j).put_block(indexQ, *ptr);
          }
        }
        input.close();
      }
    }
  }
  for (auto&& i : data_)
    i.broadcast_pack(0);
  return;
}
