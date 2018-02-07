#include <numeric>
#include "../mpi/mpi_interface.h"
#include "mps2.h"
#include "compress2.h"

using namespace std;


MPS2::MPS2(Type t) : enable_shared_from_this<MPS2>(), type_(t) {
}

MPS2::MPS2(Type t, const vector<shared_ptr<const Etensor>>& o, bool n) : enable_shared_from_this<MPS2>(), type_(t), normalized_(n) {
  assert(o.size() != 0);
  data_.reserve(o.size());
  size_t i1, i2;
  switch (type()) {
    case U:{
      i1 = 2; i2 = 3;
      break;
    }
    case D:{
      i1 = 4; i2 = 5;
      break;
    }
    default:
      throw logic_error("MPS2::MPS2 failed");
  }
  //Construct Boundary MPS
  Basis bO(Basis::OUT, 1u);
  Basis bI(Basis::IN, 1u);
  for_each(o.begin(), o.end(), [&](const shared_ptr<const Etensor>& e){
    assert(e->data().rank() != 0);
    assert(e->data().block_size() != 0);
    auto&& edata = e->data();
    Qtensor tmp({bO, edata.index(i1).reverse(), edata.index(i2).reverse(), bI});
    tmp.fill(1.0);
    data_.push_back(move(tmp));
  });
}

MPS2::MPS2(const MPS2& o) : enable_shared_from_this<MPS2>(), type_(o.type_), normalized_(o.normalized_), center_(o.center_), data_(o.data_) {
}

MPS2::MPS2(MPS2&& o) : enable_shared_from_this<MPS2>(), type_(o.type_), normalized_(o.normalized_), center_(o.center_), data_(move(o.data_)) {
  o.normalized_ = false;
}

MPS2& MPS2::operator = (const MPS2& o) {
  if (this != &o) {
    normalized_ = o.normalized_;
    center_ = o.center_;
    data_ = o.data_;
  }
  return *this;
}

MPS2& MPS2::operator = (MPS2&& o) {
  if (this != &o) {
    normalized_ = o.normalized_;
    center_ = o.center_;
    data_ = move(o.data_);
    o.normalized_ = false;
  }
  return *this;
}

void MPS2::flip_type() {
  switch (type()) {
    case U:{
      type_ = D;
      break;
    }
    case D:{
      type_ = U;
      break;
    }
    default:
      throw logic_error("MPS2::operator = failed");
  } 
}


double MPS2::left_normalize(const bool sgn) { 
  double coeff = 1.0;
  if (!normalized()) {
    const size_t na = data_.front().rank()-3;
    vector<size_t> a0(na);
    iota(a0.begin(), a0.end(), 0u);
    vector<size_t> a1(na);
    iota(a1.begin(), a1.end(), 1u);

    for (auto it = data_.begin(); it != data_.end(); ++it) {
      Qtensor tmp;
      {
        //svd//
        vector<Qtensor> ans;
        if (sgn) {
          ans = it->svd(3, Basis::OUT, "U");
        } else {
          ans = it->svd(3, Basis::IN, "U");
        }
        *it = move(ans[0]);
        it->flush();
        //rescale sigular values//
        const double max_s = ans[1].max();
        if (max_s > 1.0E3) {
          const double k = round(log2(max_s));
          const double scal = pow(2,k);
          const double scal_ = pow(2,-k);
          ans[1].scale(scal_);
          coeff *= scal;
        }
        tmp.contract(1.0, ans[1], {1}, ans[2], {0}, sgn);
      }
      //multiply to next site//
      if (it+1 != data_.end()) {
        const Qtensor hold = move(*(it+1));
        (it+1)->contract(1.0, tmp, a1, hold, a0, sgn);
        (it+1)->flush();
      } else {
        assert(tmp.get_block()->size() == 1);
        coeff *= tmp.get_block()->element(0,0);
      }
    }
    normalized_ = true;

  } else {

    assert(center_ < length());

    for(auto it = next(data_.begin(), center_); it != data_.end(); ++it) {
      Qtensor tmp;
      {
        //svd//
        vector<Qtensor> ans;
        if (sgn) {
          ans = it->svd(3, Basis::OUT, "U");
        } else {
          ans = it->svd(3, Basis::IN, "U");
        }
        *it = move(ans[0]);
        it->flush();
      
        //rescale sigular values//
        const double max_s = ans[1].max();
        if (max_s > 1.0E3) {
          const double k = round(log2(max_s));
          const double scal = pow(2,k);
          const double scal_ = pow(2,-k);
          ans[1].scale(scal_);
          coeff *= scal;
        }
        tmp.contract(1.0, ans[1], {1}, ans[2], {0}, sgn);   
      }
      //multiply to next site//
      if (it+1!=data_.end()) {
        const Qtensor hold = move(*(it+1));
        (it+1)->contract(1.0, tmp, {1}, hold, {0}, sgn);
        (it+1)->flush();
      } else {
        assert(tmp.get_block()->size() == 1);
        coeff *= tmp.get_block()->element(0,0);
      }
    }
  }
  center_ = length()-1;
  return coeff;
}

/*
double MPS2::right_normalize(const bool sgn) {
  double coeff = 1.0;
  if (!normalized()) {
    const size_t na = data_.front().rank()-3;
    vector<size_t> a0(na);
    iota(a0.begin(), a0.end(), 0u);
    for (auto it = data_.rbegin(); it != data_.rend(); ++it) {
      const size_t idx = (it+1 != data_.rend()) ? na : 1;
      vector<Qtensor> ans;
      //svd//
      if (sgn) {
        ans = it->svd(idx, Basis::OUT, "V");
      } else {
        ans = it->svd(idx, Basis::IN, "V");
      }
      *it = move(ans[2]);

      //rescale sigular values//
      const double max_s = ans[1].max();
      if (max_s > 1.0E5) {
        const double k = round(log2(max_s));
        const double scal = pow(2,k);
        const double scal_ = pow(2,-k);
        ans[1].scale(scal_);
        coeff *= scal;
      }   
      
      //multiply to next site//
      if (it+1 != data_.rend()) {
        Qtensor tmp;
        tmp.contract(1.0, ans[0], {na}, ans[1], {0}, sgn);
        size_t nr = (it+1)->rank();
        vector<size_t> a1(na);
        iota(a1.begin(), a1.end(), nr-na);
        (it+1)->contract(1.0, *(it+1), a1, tmp, a0, sgn);
      } else {
        Qtensor tmp;
        tmp.contract(1.0, ans[0], {1}, ans[1], {0}, sgn);
        assert(tmp.get_block()->size() == 1);
        coeff *= tmp.get_block()->element(0,0);
      }
    }
    normalized_ = true;
  } else {
    assert(center_ < length());
    const size_t shift = length()-center_-1;
    for (auto it = next(data_.rbegin(), shift); it != data_.rend(); ++it) {
      vector<Qtensor> ans;
      //svd//
      if (sgn) {
         ans = it->svd(1, Basis::OUT, "V");
      } else {
         ans = it->svd(1, Basis::IN, "V");
      }
      *it = move(ans[2]);

      //rescale sigular values//
      const double max_s = ans[1].max();
      if (max_s > 1.0E5) {
        const double k = round(log2(max_s));
        const double scal = pow(2,k);
        const double scal_ = pow(2,-k);
        ans[1].scale(scal_);
        coeff *= scal;
      }   
     
      //multiply to next site//
      Qtensor tmp;
      tmp.contract(1.0, ans[0], {1}, ans[1], {0}, sgn);
      if (it+1 != data_.rend()) {
        (it+1)->contract(1.0, *(it+1), {3}, tmp, {0}, sgn);
      } else {
        assert(tmp.get_block()->size() == 1);
        coeff *= tmp.get_block()->element(0,0);
      }
    }
  }
  center_ = 0;
  return coeff;
}
*/

void MPS2::mixed_normalize(const size_t i, const bool sgn) {
  assert(normalized());
  assert(i < length());
  assert(center_ < length());

  if (center_ < i) {
    for (auto it = next(data_.begin(), center_); it != next(data_.begin(), i) && it+1 != data_.end() ; ++it) {
      Qtensor tmp;
      {
        vector<Qtensor> ans;
        if (sgn) {
          ans = it->svd(3, Basis::OUT, "U");
        } else {
          ans = it->svd(3, Basis::IN, "U");
        }
        *it = move(ans[0]);
        it->flush();
        tmp.contract(1.0, ans[1], {1}, ans[2], {0}, sgn);
      }
      const Qtensor hold = move(*(it+1));
      (it+1)->contract(1.0, tmp, {1}, hold, {0}, sgn);
      (it+1)->flush();
    }
  } else if (center_ > i) {
    const size_t shift = length() - center_ - 1;
    const size_t shift2 = length() - i - 1;
    for (auto it = next(data_.rbegin(), shift); it != next(data_.rbegin(), shift2) && it+1!=data_.rend(); ++it) {
      Qtensor tmp;
      {
        vector<Qtensor> ans;
        if (sgn) {
          ans = it->svd(1, Basis::OUT, "V");
        } else {
          ans = it->svd(1, Basis::IN, "V");
        }
        *it = move(ans[2]);
        it->flush();
        tmp.contract(1.0, ans[0], {1}, ans[1], {0}, sgn);
      }
      const Qtensor hold = move(*(it+1));
      (it+1)->contract(1.0, hold, {3}, tmp, {0}, sgn);
      (it+1)->flush();
    }
  }
  center_ = i;
}

shared_ptr<MPS2> MPS2::product(const vector<shared_ptr<const Etensor>>& o) const {
  const size_t nl = length();
  assert(o.size() == nl);
  auto out = make_shared<MPS2>();
  out->type_ = type();
  vector<size_t> a;
  vector<size_t> fswaps;
  switch (type()) {
    case U:{
      for (size_t i = 0; i != nl; ++i) {
        Qtensor tmp;
        tmp.contract(1.0, site(i), {1u}, o[i]->conj(), {2u});
        tmp.contract(1.0, tmp, {1u,6u}, o[i]->wfn(), {1u,2u});
        out->data_.emplace_back(tmp.permute({0u,4u,5u,3u,6u,1u,2u,7u}));
        out->data_.back().flush();
      }
      for (size_t i = 0; i != nl; ++i) {
        if (this->site(i).parity()) {
          fswaps.push_back(i);
        }
        if (o[i]->data().parity()) {
          fswaps.push_back(i+nl);
        }
      }
      break;
    }
    case D:{
      for (size_t i = 0; i != nl; ++i) {
        Qtensor tmp;
        tmp.contract(1.0, o[i]->conj(), {1u}, site(i), {1u});
        tmp.contract(1.0, tmp, {3u,5u}, o[i]->wfn(), {2u,3u});
        out->data_.emplace_back(tmp.permute({2u,5u,3u,1u,6u,0u,7u,4u}));
        out->data_.back().flush();
      }
      for (size_t i = 0; i != nl; ++i) {
        if (o[i]->data().parity()) {
          fswaps.push_back(i);
        }
        if (this->site(i).parity()) {
          fswaps.push_back(i+nl);
        }
      }
      break;
    }
    default:
      throw logic_error("MPS2::product failed");
  }

  //contract boundary indices
  const bool change_sgn = sgn_permute(fswaps);
  //auto&& o1 = out->data_.front();
  //Qtensor q1({o1.index(0).reverse(), o1.index(1).reverse()}, true);
  //q1.fill(1.0);
  //o1.contract(1.0, q1, {0,1}, o1, {0,1});
  //o1.flush();
  {
    const Qtensor hold = move(out->data_[0]);
    Qtensor q1({hold.index(0).reverse(), hold.index(1).reverse()}, true);
    q1.fill(1.0);
    out->data_[0].contract(1.0, q1, {0,1}, hold, {0,1});
    out->data_[0].flush();
  }
  //auto&& ol = out->data_.back();
  //Qtensor ql({ol.index(5).reverse(), ol.index(6).reverse()}, true);
  //if (change_sgn) {
  //  ql.fill(-1.0);
  //} else {
  //  ql.fill( 1.0);
  //}
  //ol.contract(1.0, ql, {0,1}, ol, {5,6});
  //ol.flush();
  {
    const Qtensor hold = move(out->data_[out->data_.size()-1]);
    Qtensor ql({hold.index(5).reverse(), hold.index(6).reverse()}, true);
    if (change_sgn)
      ql.fill(-1.0);
    else
      ql.fill( 1.0);
    out->data_[out->data_.size()-1].contract(1.0, ql, {0,1}, hold, {5,6});
    out->data_[out->data_.size()-1].flush();
  }
  return out;
}

shared_ptr<MPS2> MPS2::compress(const map<Qnum, size_t>& qlist, const double thresh, const size_t max_iter, const bool mute) {
  const size_t new_dim = accumulate(qlist.begin(), qlist.end(), 0, [](size_t value, const std::map<Qnum, size_t>::value_type& e){ return value + e.second; });
  const size_t max_dim = site(length()/2).index()[0].total_dim(); //bond dimension of the center site
  
  if (max_dim > new_dim) {
    Compress2 obj(share(), qlist, mute);
    obj.iteration(thresh, max_iter);
    return obj.output();
  }
  return share();
}

ostream& operator << (ostream& os, const MPS2& o) {
  if (o.length() != 0) {
    size_t k = 0;
    for (auto&& i : o.data()) {
      os << "*** MPS2 SITE " << k++ << " : " << o.length() << " ***" << "\n" << i;
    }
  } else {
    os << "MPS2 length is 0" << "\n";
  }
  return os; 
}

void MPS2::broadcast(int source) {
  if (mpi__->size() != 1) { 
    const int iproc = mpi__->world_rank();

#ifndef NDEBUG
    if (iproc != source) {
      assert(data_.size() == 0);
    }  else {
      assert(data_.size() != 0);
    }
#endif

    // Type //
    {
      int t = static_cast<int>(type_);
      mpi__->broadcast(&t, 1u, source);
      if (iproc != source)
        type_ = static_cast<Type>(t);
    }

    // Normalized //  
    if (iproc != source) {
      normalized_ = true;
    } else {
      assert(normalized_);
    }

    // Length //
    { 
      size_t nl = (iproc == source) ? length() : 0u;
      mpi__->broadcast(&nl, 1u, source);
      if (iproc != source)
        data_.resize(nl);
    }
  
    // Center // Must be left_normalized
    if (iproc != source) {
      center_ = length()-1;
    } else {
      assert(center_ == length()-1);
    }
  
    // Data //
    for (size_t i = 0; i != length(); ++i)
       data_[i].broadcast(source);

    for (size_t i = 0; i != length(); ++i)
      assert(data_[i].disk() == true);
  }
  return;
}

bool MPS2::sgn_permute(const vector<size_t>& fidx){
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
