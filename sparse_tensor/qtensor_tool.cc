#include <iostream>
#include <functional>
#include <utility>
#include <vector>
#include <set>
#include <numeric>
#include "basis.h"
#include "qtensor.h"

using namespace std;


bool Qtensor::operator == (const Qtensor& obj) const {
  bool out = (charge_ == obj.charge_ && index_ == obj.index_ && block_size() == obj.block_size());
  if (out) {
    load();
    obj.load();
    auto it1 = data_.begin();
    auto it2 = obj.data_.begin();
    for(; it1!=data_.end() && it2!=obj.data_.end(); ++it1, ++it2) {
      out &= (it1->first == it2->first && it1->second->is_equal(*it2->second));
      if (!out) 
        break;
    }
    dump();
    obj.dump();
  }
  return out;
}

Qtensor& Qtensor::operator += (const Qtensor& obj) {
  assert(charge_ == obj.charge_);
  assert(index_ == obj.index_);

  load();
  obj.load();
  for (auto&& i : obj.data_) {
    auto it = data_.find(i.first);
    if (it != data_.end()) {
      it->second->ax_plus_y(1.0, *i.second);
    } else {
      data_.emplace(i.first, i.second->copy());
    }
  }
  modified_ = true;
  dump();
  obj.dump();

  return *this;
}

Qtensor& Qtensor::operator -= (const Qtensor& obj) { 
  assert(charge_ == obj.charge_);
  assert(index_ == obj.index_);

  load();
  obj.load();
  for (auto&& i : obj.data_) {
    shared_ptr<Tensor> minus_copy = i.second->copy();
    minus_copy->scale(-1.0);
    
    auto it = data_.find(i.first);
    if (it != data_.end()) {
      it->second->ax_plus_y(1.0, *minus_copy);
    } else {
      data_.emplace(i.first, move(minus_copy));
    }
  }
  modified_ = true;
  dump();
  obj.dump();
  return *this;
}

vector<Qtensor> Qtensor::geigs(vector<double>& ene, const Qtensor& rhs, const size_t nstate) const {
  assert(charge_==Qnum::zero());
  const size_t rank2 = rank()/2;

  #ifndef NDEBUG
  assert(rank()%2 == 0);
  assert(index_ == rhs.index_);
  auto j = rhs.index_.begin();
  for (auto i = index_.begin(); i != next(index_.begin(),rank2); ++i, ++j) {
    assert(i->type() != next(i, rank2)->type());
    assert(i->basis() == next(i, rank2)->basis());
  }
  #endif

  //Sort by Qnum & Build block matrix accordingly
  vector<vector<Qnum>> qkeys;
  transform(data_.begin(), data_.end(), back_inserter(qkeys), 
    [](const map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair) { return pair.first; 
  });

  set<vector<Qnum>> qc_set;
  while (qkeys.size() != 0) {
    auto it = qkeys.begin();
    Qnum qn = Qnum::zero();
    for( size_t i = 0; i != rank2; ++i) {
      if (index()[i].type()==Basis::IN) {
        qn += it->operator[](i);
      } else if (index()[i].type()==Basis::OUT) {
        qn -= it->operator[](i);
      }
    }
    //Search for blocks with Qnum qn equals zero
    if (qn == Qnum::zero()) {
      vector<Qnum> c(it->begin(), it->begin()+rank2);
      qc_set.insert(c);
    }
    qkeys.erase(it);
  }

  //Block EIG_SYM
  vector<shared_ptr<const Tensor>> block_ham;
  vector<shared_ptr<const Tensor>> block_rhs;
  for (auto&& i : qc_set) {
    vector<Qnum> qid(i);
    for (auto&& j : qc_set) {
      vector<Qnum> qid2(i);
      qid2.insert(qid2.end(), j.begin(), j.end());
      block_ham.push_back(get_block(qid2));
      block_rhs.push_back(rhs.get_block(qid2));
    }
  }

  const size_t qcols = qc_set.size();
  //Calculate col offsets
  vector<size_t> coffsets(qcols+1);
  vector<vector<size_t>> cextents(qcols);
  size_t ns = 0;
  for (size_t i = 0; i != qcols; ++i) {
    coffsets[i] = ns;
    size_t nc = 0;
    for (auto k = block_ham.begin()+i; distance(block_ham.begin(),k) < (int)block_ham.size(); k += qcols) {
      if (*k) {
        cextents[i].assign((*k)->extent().begin()+rank2, (*k)->extent().end());
        nc = accumulate((*k)->extent().begin()+rank2, (*k)->extent().end(), 1, multiplies<size_t>());
        break;
      }
    } 
    assert(nc != 0);
    ns += nc; 
  }
  coffsets[qcols] = ns;
    
  //Construct Matrix & Diagonalize
  Matrix ham(coffsets.back(), coffsets.back());
  Matrix norm(coffsets.back(), coffsets.back());
  auto it_ham = block_ham.begin();
  auto it_rhs = block_rhs.begin();
  //Row-wise loop//
  for (size_t i = 0; i != qcols; ++i) {
    for (size_t j = 0; j != qcols; ++j) {
      if (*it_ham)
        ham.copy_block(coffsets[i], coffsets[j], coffsets[i+1] - coffsets[i], coffsets[j+1] - coffsets[j], (*it_ham)->pdata());
      ++it_ham;
      if (*it_rhs)
        norm.copy_block(coffsets[i], coffsets[j], coffsets[i+1] - coffsets[i], coffsets[j+1] - coffsets[j], (*it_rhs)->pdata());
      ++it_rhs;
    }
    
  }
  ham.symmetrize();
  norm.symmetrize();
  Matrix wfn;
  auto eigenvalues = ham.geig_sym(norm, wfn, 1.0E-13);
  if (eigenvalues.size() < nstate)
    std::runtime_error("Qtensor::geigs failed");
  
  //Return
  ene.resize(nstate);
  copy_n(eigenvalues.begin(), nstate, ene.begin());
  //Compute out
  vector<Qtensor> out;
  for (size_t istate = 0; istate != nstate; ++istate) {
    shared_ptr<Matrix> psi = wfn.cut_col(istate, istate+1);
    Qtensor iout;
    iout.index_.assign(index().begin(), index().begin()+rank2);
    auto itq = qc_set.begin();
    for (size_t i = 0; i != qcols; ++i, ++itq) {
      const auto& m = psi->cut_row(coffsets[i], coffsets[i+1]);
      auto ext = cextents[i];
      auto tptr = make_shared<Tensor>(ext);
      tptr->set_elements(m->size(), m->pdata());
      iout.data_.emplace(*itq, move(tptr));
    }
    assert(iout.data_.size() != 0);
    out.push_back(move(iout));
  }
  return out;
}

vector<Qtensor> Qtensor::split(const int order) const {
  #ifndef NDEBUG
  assert(charge() == Qnum::zero());
  assert(rank() == 2u);
  assert(index(0u).type() != index(1u).type());
  assert(index(0u).basis() == index(1u).basis());
  #endif

  Basis new_basis(index(0u).type());
  Qtensor out0({index(1u), new_basis}, false);
  Qtensor out1({new_basis.reverse(), index(0u)}, false);


  load(); 
  for (const auto& e : data_) {
    Matrix mat(e.second->extent(0), e.second->extent(1));
    mat.set_elements(e.second->size(), e.second->pdata());
    vector<shared_ptr<Matrix>> ans = mat.split(order);

    const auto& a0 = ans[0];
    const auto& a1 = ans[1];
    const Qnum qn = e.first[0];
    assert(e.first[1] == qn);

    const size_t qdim = a0->ncols();
    assert(a1->nrows() == qdim);
    out0.index_.back().insert(qn, qdim);
    out1.index_.front().insert(qn, qdim);
    
    auto t0 = make_shared<Tensor>(a0->nrows(), a0->ncols());
    t0->set_elements(a0->size(), a0->pdata());
    out0.data_.emplace(e.first, move(t0));    

    auto t1 = make_shared<Tensor>(a1->nrows(), a1->ncols());
    t1->set_elements(a1->size(), a1->pdata());
    out1.data_.emplace(e.first, move(t1));
  }
  dump();
  return {out0, out1};
}

vector<Qtensor> Qtensor::sqrt2(vector<double>& valueS) const {
#ifndef NDEBUG
  assert(charge() == Qnum::zero());
  assert(rank() == 2u);
  assert(index(0u).type() != index(1u).type());
#endif

  Basis new_basis(index(1u).type());
  Qtensor out0({index(0u), new_basis}, false);
  Qtensor out1({new_basis.reverse(), index(1u)}, false);
  
  valueS.clear();
  load();
  for (const auto& e : data_) {
    Matrix mat(e.second->extent(0), e.second->extent(1));
    mat.set_elements(e.second->size(), e.second->pdata());
    vector<double> tmp;
    vector<shared_ptr<Matrix>> ans = mat.sqrt2(tmp);
    valueS.insert(valueS.end(), tmp.begin(), tmp.end());

    const auto& a0 = ans[0];
    const auto& a1 = ans[1];
    const Qnum qn = e.first[0];
    assert(e.first[1] == qn);

    const size_t qdim = a0->ncols();
    assert(a1->nrows() == qdim);
    out0.index_.back().insert(qn, qdim);
    out1.index_.front().insert(qn, qdim);

    auto t0 = make_shared<Tensor>(a0->nrows(), a0->ncols());
    t0->set_elements(a0->size(), a0->pdata());
    out0.data_.emplace(e.first, move(t0));

    auto t1 = make_shared<Tensor>(a1->nrows(), a1->ncols());
    t1->set_elements(a1->size(), a1->pdata());
    out1.data_.emplace(e.first, move(t1));
  }
  dump();
  return {out0, out1};
}

void Qtensor::scale(const double a) {
  load();
  for_each(data_.begin(), data_.end(), [&](map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair){
    pair.second->scale(a);   
  });
  dump();
}


/* Not use at this moment
Qtensor Qtensor::eigs(double& ene, const Qtensor& guess) const {
  //TODO Seperate Quantum Number
  assert(charge() == Qnum::zero());
  size_t idx = rank()/2;
  #ifndef NDEBUG
  assert(rank()%2 == 0);
  for (auto i = index_.begin(); i != next(index_.begin(),idx); ++i) {
    assert(i->type() != next(i,idx)->type());
    assert(i->basis() == next(i,idx)->basis());
  }
  assert(idx == guess.rank());
  for (size_t i = 0; i != idx; ++i) {
    assert(index_[i] == guess.index(i));
  }
  #endif

  Qtensor out;
  out.index_.assign(index().begin(), index().begin()+idx);
  //Sort by Qnum & Build block matrix accordingly
  vector<vector<Qnum>> qkeys;
  transform(data_.begin(), data_.end(), back_inserter(qkeys), 
    [](const map<vector<Qnum>, shared_ptr<Tensor>>::value_type& pair) { return pair.first; 
  });
  set<vector<Qnum>> qc_set;
  while (qkeys.size() != 0) {
    auto it = qkeys.begin();
    Qnum qn = Qnum::zero();
    for (size_t i = 0; i != idx; ++i) {
      if (index()[i].type() == Basis::IN) {
        qn += it->operator[](i);
      } else if (index()[i].type() == Basis::OUT) {
        qn -= it->operator[](i);
      }   
    }
    //Search for blocks with Qnum qn equals zero
    if (qn == Qnum::zero()) {
      vector<Qnum> c(it->begin(), it->begin()+idx);
      qc_set.insert(c);
    }
    qkeys.erase(it);
  }

  //Block EIG_SYM
  vector<shared_ptr<const Tensor>> block_ham;
  vector<shared_ptr<const Tensor>> block_psi;
  for (auto&& i : qc_set) {
    vector<Qnum> qid(i);
    block_psi.push_back(guess.get_block(qid));
    for (auto&& j : qc_set) {
      vector<Qnum> qid2(i);
      qid2.insert(qid2.end(), j.begin(), j.end());
      block_ham.push_back(get_block(qid2));
    }
  }

  const size_t qcols = qc_set.size();
  //Calculate col offsets
  vector<size_t> coffsets(qcols+1);
  vector<vector<size_t>> cextents(qcols);
  size_t ns = 0;
  for (size_t i = 0; i != qcols; ++i) {
    coffsets[i] = ns;
    size_t nc = 0;
    for (auto k = block_ham.begin()+i; distance(block_ham.begin(),k) < (int)block_ham.size(); k += qcols) {
      if (*k) {
        cextents[i].assign((*k)->extent().begin()+idx, (*k)->extent().end());
        nc = accumulate((*k)->extent().begin()+idx, (*k)->extent().end(), 1, multiplies<size_t>());
        break;
      }
    } 
    assert(nc != 0);
    ns += nc; 
  }
  coffsets[qcols] = ns;
    
  //Construct Matrix & Diagonalize
  Matrix ham(coffsets.back(), coffsets.back());
  Matrix psi(coffsets.back(),1);
  auto it_ham = block_ham.begin();
  auto it_psi = block_psi.begin();
  //Row-wise loop//
  for (size_t i = 0; i != qcols; ++i) {
    if (*it_psi)
      psi.copy_block(coffsets[i], 0, coffsets[i+1] - coffsets[i], 1, (*it_psi)->pdata());
    ++it_psi;
    for (size_t j = 0; j != qcols; ++j){
      if (*it_ham)
        ham.copy_block(coffsets[i], coffsets[j], coffsets[i+1] - coffsets[i], coffsets[j+1] - coffsets[j], (*it_ham)->pdata());
      ++it_ham;
    }
  }
  auto eigenvalues = ham.eig_davidson(psi,10);
  ene = eigenvalues[0];

  //compute out
  auto itq = qc_set.begin();
  for (size_t i = 0; i != qcols; ++i, ++itq) {
    const auto& m = psi.cut_row(coffsets[i], coffsets[i+1]);
    auto ext = cextents[i];
    auto tptr = make_shared<Tensor>(ext);
    tptr->set_elements(m->size(), m->pdata());
    out.insert_block(*itq, tptr);
  }
  assert(out.data_.size() != 0);
  return out;
}
*/
