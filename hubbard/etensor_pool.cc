#include "etensor_pool.h"

using namespace std;

Hubbard2::EtensorPool::EtensorPool (const shared_ptr<const PEPS>& psi) : psi_(psi) {
  init();
  const size_t nr = psi_->nrows();
  const size_t nc = psi_->ncols();
  for (size_t i = 0; i != nr; ++i) {
    for (size_t j = 0; j != nc; ++j) {
      alphaC_.emplace_back(make_shared<Etensor>(psi_->site(i,j), op_[0]));
      alphaA_.emplace_back(make_shared<Etensor>(psi_->site(i,j), op_[1]));
      betaC_.emplace_back(make_shared<Etensor>(psi_->site(i,j), op_[2]));
      betaA_.emplace_back(make_shared<Etensor>(psi_->site(i,j), op_[3]));
      doublon_.emplace_back(make_shared<Etensor>(psi_->site(i,j), op_[4]));
      identity_.emplace_back(make_shared<Etensor>(psi_->site(i,j), op_[5]));
    }
  }
}

void Hubbard2::EtensorPool::update(const size_t i, const size_t j) const {
  assert(i < psi_->nrows() && j < psi_->ncols());
  alphaC(i,j)->update();
  alphaA(i,j)->update();
  betaC(i,j)->update();
  betaA(i,j)->update();
  doublon(i,j)->update();
  identity(i,j)->update();
  return;
}

vector<shared_ptr<const Etensor>> Hubbard2::EtensorPool::get_row_alphaC(const size_t irow, const size_t jcol) const {
  const size_t nc = psi_->ncols();
  assert(jcol < nc);
  vector<shared_ptr<const Etensor>> out(nc);
  for (size_t j = 0; j != nc; ++j) {
    out[j] = (j==jcol) ? alphaC(irow, j) : identity(irow,j);
  }
  return out;
}

vector<shared_ptr<const Etensor>> Hubbard2::EtensorPool::get_row_alphaA(const size_t irow, const size_t jcol) const {
  const size_t nc = psi_->ncols();
  assert(jcol < nc);
  vector<shared_ptr<const Etensor>> out(nc);
  for (size_t j = 0; j != nc; ++j) {
    out[j] = (j==jcol) ? alphaA(irow, j) : identity(irow,j);
  }
  return out;
}

vector<shared_ptr<const Etensor>> Hubbard2::EtensorPool::get_row_betaC(const size_t irow, const size_t jcol) const {
  const size_t nc = psi_->ncols();
  assert(jcol < nc);
  vector<shared_ptr<const Etensor>> out(nc);
  for (size_t j = 0; j != nc; ++j) {
    out[j] = (j==jcol) ? betaC(irow, j) : identity(irow,j);
  }
  return out;
}

vector<shared_ptr<const Etensor>> Hubbard2::EtensorPool::get_row_betaA(const size_t irow, const size_t jcol) const {
  const size_t nc = psi_->ncols();
  assert(jcol < nc);
  vector<shared_ptr<const Etensor>> out(nc);
  for (size_t j = 0; j != nc; ++j) {
    out[j] = (j==jcol) ? betaA(irow, j) : identity(irow,j);
  }
  return out;
}

vector<shared_ptr<const Etensor>> Hubbard2::EtensorPool::get_row_doublon(const size_t irow, const size_t jcol) const {
  const size_t nc = psi_->ncols();
  assert(jcol < nc);
  vector<shared_ptr<const Etensor>> out(nc);
  for (size_t j = 0; j != nc; ++j) {
    out[j] = (j==jcol) ? doublon(irow, j) : identity(irow,j);
  }
  return out;
}

vector<shared_ptr<const Etensor>> Hubbard2::EtensorPool::get_row_identity(const size_t irow) const {
  const size_t nc = psi_->ncols();
  vector<shared_ptr<const Etensor>> out(nc);
  for (size_t j = 0; j != nc; ++j) {
    out[j] = identity(irow,j);
  }
  return out;
}

vector<shared_ptr<const Etensor>> Hubbard2::EtensorPool::get_row_alphaC_alphaA(const size_t irow, const size_t jcol1, const size_t jcol2) const {
  assert(jcol2 < psi_->ncols() && jcol1 < jcol2);
  vector<shared_ptr<const Etensor>> out(get_row_identity(irow));
  out[jcol1] = alphaC(irow, jcol1);
  out[jcol2] = alphaA(irow, jcol2);
  return out;
}

vector<shared_ptr<const Etensor>> Hubbard2::EtensorPool::get_row_alphaA_alphaC(const size_t irow, const size_t jcol1, const size_t jcol2) const {
  assert(jcol2 < psi_->ncols() && jcol1 < jcol2);
  vector<shared_ptr<const Etensor>> out(get_row_identity(irow));
  out[jcol1] = alphaA(irow, jcol1);
  out[jcol2] = alphaC(irow, jcol2);
  return out;
}

vector<shared_ptr<const Etensor>> Hubbard2::EtensorPool::get_row_betaC_betaA(const size_t irow, const size_t jcol1, const size_t jcol2) const {
  assert(jcol2 < psi_->ncols() && jcol1 < jcol2);
  vector<shared_ptr<const Etensor>> out(get_row_identity(irow));
  out[jcol1] = betaC(irow, jcol1);
  out[jcol2] = betaA(irow, jcol2);
  return out;
}

vector<shared_ptr<const Etensor>> Hubbard2::EtensorPool::get_row_betaA_betaC(const size_t irow, const size_t jcol1, const size_t jcol2) const {
  assert(jcol2 < psi_->ncols() && jcol1 < jcol2);
  vector<shared_ptr<const Etensor>> out(get_row_identity(irow));
  out[jcol1] = betaA(irow, jcol1);
  out[jcol2] = betaC(irow, jcol2);
  return out;
}

void Hubbard2::EtensorPool::init() {
  Qnum q0(0, true);
  Qnum q1(1, true);
  Qnum q2(2, true);
  Basis s(Basis::OUT, {q0,q1,q1,q2});

  Qtensor alphaC(q1, {s, s.reverse()}, false);
  auto ptr = make_shared<Tensor>(2u,1u);
  ptr->element(0,0) =  1.0;
  alphaC.insert_block({q1,q0}, *ptr);
  ptr = make_shared<Tensor>(1u,2u);
  ptr->element(0,1) =  1.0;
  alphaC.insert_block({q2,q1}, *ptr);

  Qtensor alphaA(-q1, {s, s.reverse()}, false);
  ptr = make_shared<Tensor>(1u,2u);
  ptr->element(0,0) =  1.0;
  alphaA.insert_block({q0,q1}, *ptr);
  ptr = make_shared<Tensor>(2u,1u);
  ptr->element(1,0) =  1.0;
  alphaA.insert_block({q1,q2}, *ptr);

  Qtensor betaC(q1, {s, s.reverse()}, false);
  ptr = make_shared<Tensor>(2u,1u);
  ptr->element(1,0) =  1.0;
  betaC.insert_block({q1,q0}, *ptr);
  ptr = make_shared<Tensor>(1u,2u);
  ptr->element(0,0) = -1.0;
  betaC.insert_block({q2,q1}, *ptr);

  Qtensor betaA(-q1, {s, s.reverse()}, false);
  ptr = make_shared<Tensor>(1u,2u);
  ptr->element(0,1) =  1.0;
  betaA.insert_block({q0,q1}, *ptr);
  ptr = make_shared<Tensor>(2u,1u);
  ptr->element(0,0) = -1.0;
  betaA.insert_block({q1,q2}, *ptr);

  Qtensor doublon({s, s.reverse()}, false);
  ptr = make_shared<Tensor>(1u,1u);
  ptr->element(0,0) = 1.0;
  doublon.insert_block({q2,q2}, *ptr);

  Qtensor identity({s, s.reverse()}, false);
  ptr = make_shared<Tensor>(1u,1u);
  ptr->element(0,0) = 1.0;
  identity.insert_block({q0,q0}, *ptr);
  identity.insert_block({q2,q2}, *ptr);
  ptr = make_shared<Tensor>(2u,2u);
  ptr->element(0,0) = 1.0;
  ptr->element(1,1) = 1.0;
  identity.insert_block({q1,q1}, *ptr);

  op_[0] = move(alphaC);
  op_[1] = move(alphaA);
  op_[2] = move(betaC);
  op_[3] = move(betaA);
  op_[4] = move(doublon);
  op_[5] = move(identity);
  return;
}
