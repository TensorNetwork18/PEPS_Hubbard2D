#include "tmatrix2.h"

using namespace std;


TMatrix2::TMatrix2(Type tp, const shared_ptr<const MPS2>& u, const shared_ptr<const MPS2>& d, const vector<shared_ptr<const Etensor>>& o, const bool f) : type_(tp), uptr_(u), dptr_(d), pmpo_(o), sgn_(f) {
  size_t nl = uptr_->length();
  assert(nl != 0u);
  assert(dptr_->length() == nl);
  assert(pmpo_.size() == nl || pmpo_.size() == 0u);
  assert(uptr_->type() == MPS2::U);
  assert(dptr_->type() == MPS2::D);
  data_.resize(nl);

  switch (type()) { //build boundary sites
    case L: {
      vector<size_t> fswaps; 
      if (!pmpo_.empty()) {
        vector<Basis> b = {uptr_->site(0u).index(0u).reverse(), pmpo_[0u]->data().index(0u).reverse(), pmpo_[0u]->data().index(1u).reverse(), dptr_->site(0u).index(0u).reverse()};
        base_ =  make_shared<Qtensor>(b, true);
        if (sgn_) {
          for (size_t i = 0; i != nl; ++i) {
            if (uptr_->site(i).parity()) fswaps.push_back(i);
            if (pmpo_[i]->data().parity()) fswaps.push_back(i+nl);
            if (dptr_->site(i).parity()) fswaps.push_back(i+2u*nl);
          }
        }
      } else {
        vector<Basis> b = {uptr_->site(0u).index(0).reverse(), dptr_->site(0u).index(0).reverse()};
        base_ =  make_shared<Qtensor>(b, true);
        if (sgn_) {
          for (size_t i = 0; i != nl; ++i) {
            if (uptr_->site(i).parity()) fswaps.push_back(i);
            if (dptr_->site(i).parity()) fswaps.push_back(i+nl);
          }
        }
      }
      (sgn_ && sgn_permute(fswaps)) ? base_->fill(-1.0) : base_->fill(1.0);
      break; 
    }
    case R: {
      const size_t s = nl-1u;
      if (!pmpo_.empty()) {
        vector<Basis> b = {uptr_->site(s).index(3u).reverse(), pmpo_[s]->data().index(6u).reverse(), pmpo_[s]->data().index(7u).reverse(), dptr_->site(s).index(3u).reverse()};
        base_ = make_shared<Qtensor>(b, true);
      } else {
        vector<Basis> b = {uptr_->site(s).index(3).reverse(), dptr_->site(s).index(3).reverse()};
        base_ = make_shared<Qtensor>(b, true);
      }
      base_->fill(1.0);
      break;
    }
    default:
      throw logic_error("Tmatrix2::Tmatrix2 failed");   
  }
}

TMatrix2::TMatrix2(Type tp, const shared_ptr<const MPS2>& u, const vector<shared_ptr<const Etensor>>& o, const shared_ptr<const MPS2>& d, const bool f) : TMatrix2(tp, u, d, o, f){
}

Qtensor& TMatrix2::evaluate(const size_t i) {
  const size_t nl = length();
  if (i < nl) {
    if (site(i).rank() != 0) {
      assert(site(i).block_size() != 0);
      return site(i);
    } else {
      switch (type()) {
        case L: {
          return evaluate_from_left(i);
          break;
        }
        case R: {
          return evaluate_from_right(i);
          break;
        }
        default:
          throw logic_error("Tmatrix2::evaluate failed"); 
      }
    }
  } else {
    if (!data_.empty()) {
      data_.clear();
      data_.resize(nl);
    }
  }
  return base();
}

Qtensor& TMatrix2::evaluate_from_left(const size_t i) {
  const auto& last = evaluate(i-1);
  Qtensor tmp;
  if(!pmpo_.empty()){
    assert(last.rank() == 4);
    tmp.contract(1.0, last, {0}, uptr_->site(i), {0}, sgn());
    tmp.contract(1.0, tmp, {0u,3u}, pmpo_[i]->conj(), {3u,2u}, sgn());
    tmp.contract(1.0, tmp, {0u,2u,6u}, pmpo_[i]->wfn(), {0u,1u,2u}, sgn());
    tmp.contract(1.0, tmp, {0u,3u,4u}, dptr_->site(i), {0u,1u,2u},sgn());
  } else{
    assert(last.rank() == 2);
    tmp.contract(1.0, last, {0}, uptr_->site(i), {0}, sgn());
    tmp.contract(1.0, tmp, {0,1,2}, dptr_->site(i), {0,1,2}, sgn());
  }
  site(i) = move(tmp);
  return site(i);
}

Qtensor& TMatrix2::evaluate_from_right(const size_t i) {
  const auto& last = evaluate(i+1);
  Qtensor tmp;
  if (!pmpo_.empty()) {
    assert(last.rank() == 4);
    tmp.contract(1.0, dptr_->site(i), {3}, last, {3}, sgn());
    tmp.contract(1.0, pmpo_[i]->conj(), {0u,1u}, tmp, {4u,1u}, sgn());
    tmp.contract(1.0, pmpo_[i]->wfn(), {2u,3u,4u}, tmp, {2u,4u,6u}, sgn());
    tmp.contract(1.0, uptr_->site(i), {1u,2u,3u}, tmp, {2u,1u,5u}, sgn());
    tmp = tmp.permute({0u,2u,1u,3u}, sgn());
  } else {
    assert(last.rank() == 2);
    tmp.contract(1.0, dptr_->site(i), {3}, last, {1}, sgn());
    tmp.contract(1.0, uptr_->site(i), {1,2,3}, tmp, {1,2,3}, sgn());
  }
  site(i) = move(tmp);
  return site(i);
}

bool TMatrix2::sgn_permute(const vector<size_t>& fidx){
  // Odd  -> True
  // Even -> False
  const size_t ns = fidx.size();
  vector<size_t> idx(ns);
  iota(idx.begin(), idx.end(), 0); 
  sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return fidx[a] < fidx[b]; }); 
  
  vector<bool> visited(ns, false);
  bool sgn = false;
  for(size_t i=0; i<ns; ++i){
    if(!visited[i]){
      size_t next = i;
      int count = 0;
      while(!visited[next]){
        visited[next] = true;
        next = idx[next];
        ++count;
      }   
      if(count%2==0)
        sgn = !sgn;
    }   
  }
  return sgn;
}
