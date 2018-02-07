#include <iostream>
#include <cmath>
#include <algorithm>
#include "tensor.h"
#include "mkl_interface.h"

using namespace std;


bool Tensor::is_equal (const Tensor& o, const double thresh) const {
  bool out = (extent() == o.extent());
  if (out){
    Tensor tmp = *this - o;
    out &= all_of(tmp.pdata(), tmp.pdata()+tmp.size(), [&](double i){return abs(i) < thresh;});
  }
  return out;
}

vector<shared_ptr<Tensor>> Tensor::svd(const size_t indx) const {
  assert(indx < rank() && indx != 0);
  size_t m = 1; for_each( extent().begin(), extent().begin()+indx, [&](size_t i) {m *= i;} );
  size_t n = 1; for_each( extent().begin()+indx, extent().end(), [&](size_t i) {n *= i;} );
  vector<size_t> ext_u(extent().begin(), extent().begin()+indx);
  vector<size_t> ext_vt(extent().begin()+indx, extent().end());
  const size_t c = (m<n) ? m : n;
  ext_u.push_back(c);
  ext_vt.insert(ext_vt.begin(), c);
  auto s = make_shared<Tensor>(c);
  auto u = make_shared<Tensor>(ext_u);
  auto vt = make_shared<Tensor>(ext_vt);
  dgesvd(m, n, *this, *s, *u, *vt);
  return {u, s, vt};
}

shared_ptr<Tensor> Tensor::permute(const vector<size_t>& indx) const {
  const size_t n = rank();
  assert(indx.size() == n);
  vector<size_t> tmp(n);

  for (size_t i = 0; i != n; ++i) tmp[i] = extent(indx[i]); 
  auto out = make_shared<Tensor>(tmp);

  const size_t s = size();
  for (size_t i = 0; i != s; ++i) {
    size_t dist = i;
    for (size_t j = 0; j != n; ++j) {
      tmp[j] = dist%extent(j);
      dist = dist/extent(j);
    }
    dist = tmp[indx[n-1]];
    for (size_t j = n-1; j-->0;) {
      dist *= out->extent(j);
      dist += tmp[indx[j]];
    }
    *(out->pdata() + dist) = *(pdata() + i);
  }
  return out;
}

void Tensor::contract(const double alpha, const Tensor& A, const vector<size_t>& aA, const Tensor& B, const vector<size_t>& aB) {
  const size_t na = aA.size();
  #ifndef NDEBUG
  assert(na == aB.size());
  for (size_t i = 0; i != na; ++i)
    assert(A.extent(aA[i]) == B.extent(aB[i]));
  #endif
  auto refA = make_shared<Tensor>();
  auto refB = make_shared<Tensor>();
  if (na != 0) {
    vector<size_t> indxA(A.rank());
    iota(indxA.begin(), indxA.end(), 0);
    vector<size_t> indxB(B.rank());
    iota(indxB.begin(), indxB.end(), 0);
    vector<size_t> removeA(aA);
    vector<size_t> removeB(aB);
    sort(removeA.begin(), removeA.end());
    sort(removeB.begin(), removeB.end());
    for (size_t i = na; i-- > 0;) {
      indxA.erase(indxA.begin()+removeA[i]);
      indxB.erase(indxB.begin()+removeB[i]);
    }
    indxA.insert(indxA.end(), aA.begin(), aA.end());
    indxB.insert(indxB.begin(), aB.begin(), aB.end());
    refA = A.permute(indxA);
    refB = B.permute(indxB);
  } else {
    refA = A.copy();
    refB = B.copy();
  }
  vector<size_t> tmp;
  const size_t nr = A.rank()+B.rank()-2*na;
  if (nr == 0) {
    tmp.resize(1);
    tmp[0] = 1;
  } else {
    tmp.reserve(nr);
    tmp.insert(tmp.end(), refA->extent().begin(), refA->extent().end()-na);
    tmp.insert(tmp.end(), refB->extent().begin()+na, refB->extent().end());
  }
  Tensor out(tmp);
  size_t m = 1; for_each( refA->extent().begin(), refA->extent().end()-na, [&](size_t i) {m *= i;} );
  size_t n = 1; for_each( refB->extent().begin()+na, refB->extent().end(), [&](size_t i) {n *= i;} );
  size_t k = 1; for_each( refB->extent().begin(), refB->extent().begin()+na, [&](size_t i) {k *= i;} );
  const double beta = 0.0;
  dgemm("N", "N", m, n, k, alpha, *refA, *refB, beta, out);
  *this = move(out);
}


void Tensor::scale(const double a) {
  dscal(a, *this);
}

double Tensor::max() const {
  double out = *max_element(pdata(), pdata()+size(), [](const double& e1, const double& e2){ 
    return e1 < e2; }); 
  return out;
}

double Tensor::max_abs() const {
  double out = *max_element(pdata(), pdata()+size(), [](const double& e1, const double& e2){ 
    return abs(e1) < abs(e2); }); 
  return abs(out);
}

double Tensor::sum_diag() const {
  double out = 0.0;
  assert(rank() == 2u);
  assert(extent(0u) == extent(1u));
  for (size_t ii = 0; ii != extent(0u); ++ii) {
    out+=element(ii,ii);
  }
  return out;
}

ostream& operator << (ostream& os, const Tensor& o) {
  if (o.pdata()) {
    const size_t ns = o.size();
    const size_t nr = o.rank();
    for (size_t i = 0; i != ns; ++i) {
      size_t dist = i;
      for (size_t j=0; j != nr; ++j) {
        os << setw(3) << fixed << dist%o.extent(j) << " ";
        dist = dist/o.extent(j);       
      }
      os << setw(15) << fixed << setprecision(8) << *(o.pdata()+i) << "\n";
    }
  } else {
    os << "It's a Null Tensor" << "\n";
  }
  return os;
}
