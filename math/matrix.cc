#include <iostream>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <cmath>
#include <map>
#include "matrix.h"
#include "mkl_interface.h"
#include "davidson.h"

using namespace std;


Matrix Matrix::operator * (const Matrix& o) const {
  const size_t m = this->nrows();
  const size_t k = this->ncols();
  const size_t n = o.ncols();
  assert(o.nrows() == k);
  Matrix out(m,n);
  const double alpha = 1.0;
  const double beta = 0.0;
  dgemm("N", "N", m, n, k, alpha, *this, o, beta, out);
  return out;
}

Matrix Matrix::operator % (const Matrix& o) const {
  const size_t m = this->ncols();
  const size_t k = this->nrows();
  const size_t n = o.ncols();
  assert(o.nrows() == k);
  Matrix out(m,n);
  const double alpha = 1.0;
  const double beta = 0.0;
  dgemm("T", "N", m, n, k, alpha, *this, o, beta, out);
  return out;
}

bool Matrix::operator == (const Matrix& o)const {
  return (nrows() == o.nrows()) && (ncols() == o.ncols()) \
    && (equal(pdata(), pdata()+size(), o.pdata()));
}

//
void Matrix::copy_block(const size_t row_start, const size_t col_start, const size_t row_size, const size_t col_size, const Matrix& o) {
  assert(row_start+row_size <= nrows());
  assert(col_start+col_size <= ncols());
  for (size_t i = col_start, j = 0; i != col_start + col_size; ++i, ++j)
    copy_n(o.element_ptr(0,j), row_size, element_ptr(row_start, i)); 
}

void Matrix::copy_block(const size_t row_start, const size_t col_start, const size_t row_size, const size_t col_size, const double* o) {
  //Caution: check size of o[] 
  assert(row_start+row_size <= nrows());
  assert(col_start+col_size <= ncols());
  for (size_t i = col_start, j = 0; i != col_start + col_size; ++i, ++j)
    copy_n(o+j*row_size, row_size, element_ptr(row_start, i));
}

shared_ptr<Matrix> Matrix::resize(const size_t row_size, const size_t col_size) {
  assert(row_size >= nrows() && col_size >= ncols());
  auto out = make_shared<Matrix>(row_size, col_size);
  for (size_t i = 0; i != ncols(); ++i)
    copy_n(element_ptr(0, i), nrows(), out->element_ptr(0,i));
  return out;
}

shared_ptr<Matrix> Matrix::get_submatrix(const size_t row_start, const size_t col_start, const size_t row_size, const size_t col_size) const {
  assert(row_start + row_size <= nrows());
  assert(col_start + col_size <= ncols());
  auto out = make_shared<Matrix>(row_size, col_size);
  for (size_t i = col_start, j = 0; i < col_start + col_size; i++, j++)
     copy_n(element_ptr(row_start, i), row_size, out->element_ptr(0, j));
  return out;
}

shared_ptr<Matrix> Matrix::cut_row(const size_t row_start, const size_t row_end) const {
  return get_submatrix(row_start, 0, row_end-row_start, ncols());
}

shared_ptr<Matrix> Matrix::cut_col(const size_t col_start, const size_t col_end) const {
  return get_submatrix(0, col_start, nrows(), col_end-col_start);
}

void Matrix::symmetrize() {
  const size_t n = nrows();
  assert(ncols() == n);
  for (size_t i = 0; i != n; ++i)
    for (size_t j = i+1; j != n; ++j)
      element(i, j) = element(j, i) = 0.5*(element(i, j)+element(j, i));
}

void Matrix::add_diag(const double a, const size_t nstart, const size_t nend) {
  assert(nstart <= nend);
  assert(nend <= nrows() && nend <= ncols());
  for (size_t ii = nstart; ii != nend; ++ii) element(ii,ii) += a;
}

void Matrix::unit() { 
  fill(0);
  const size_t n = min(nrows(), ncols());
  add_diag(1.0, 0, n);
}

bool Matrix::is_symmetric(const double thresh) const {
  shared_ptr<Matrix> tmp = transpose();
  *tmp -= *this;
  return tmp->rms() < thresh;
}

bool Matrix::is_antisymmetric(const double thresh) const {
  shared_ptr<Matrix> tmp = transpose();
  *tmp += *this;
  return tmp->rms() < thresh;
}

bool Matrix::is_identity(const double thresh) const {
  shared_ptr<Matrix> tmp = clone();
  tmp->unit();
  *tmp -= *this;
  return tmp->rms() < thresh;
}

double Matrix::dot_product(const Matrix& o) const {
  return ddot(*this, o);
}

double Matrix::norm() const { 
  return sqrt(dot_product(*this));
} 

double Matrix::variance() const { 
  return dot_product(*this) / (nrows() * ncols()); 
}

double Matrix::rms() const { 
  return sqrt(variance()); 
}

double Matrix::sum_diag() const {
  double out = 0.0;
  assert(nrows() == ncols());
  for (size_t ii = 0; ii != nrows(); ++ii) {
    out+=element(ii,ii);
  }
  return out;
}

void Matrix::scale(const double a) {
   dscal(a, *this); 
}

shared_ptr<Matrix> Matrix::transpose() const {
  return transpose_mkl(nrows(), ncols(), *this);
}

void Matrix::eig_sym(Vec& eig) {
  const size_t n = nrows();
  assert(ncols() == n);
  Vec tmp(n);
  dsyevd(n ,*this, tmp);
  eig = move(tmp);
}

shared_ptr<Matrix> Matrix::tildex(const double thresh) const {
  shared_ptr<Matrix> out = copy();
  const size_t n = ncols();
  Vec eig(n);
  out->eig_sym(eig);
  size_t m = 0;
  for (size_t i = 0; i != n; ++i) {
    if (eig(i) > thresh) {
      const double e = 1.0/sqrt(eig(i));
      transform(out->element_ptr(0,i), out->element_ptr(0,i+1), out->element_ptr(0,m++), [&e](double a){ return a*e; });
    }
  }
  out = out->cut_col(0,m);
  return out;
}

vector<shared_ptr<Matrix>> Matrix::split(const int a, const double thresh) const {
  //return VD^(1/a) D^(1/a)Vt
  assert(is_symmetric());
  const size_t n = nrows();
  assert(ncols() == n); 
  Vec eig(n);
  shared_ptr<Matrix> out1 = copy();
  dsyevd(n ,*out1, eig);
  shared_ptr<Matrix> out2 = out1->copy();
  size_t m = 0;
  for (size_t i = 0; i != n; ++i) {
    if (eig(i) > thresh) {
      //const double f = sqrt(sqrt(eig(i)));
      const double f = pow(eig(i), 1.0/a);
      transform(out1->element_ptr(0,i), out1->element_ptr(0,i+1), out1->element_ptr(0,m), [&f](double a){ return a*f; }); 
      const double e = 1.0/f;
      transform(out2->element_ptr(0,i), out2->element_ptr(0,i+1), out2->element_ptr(0,m++), [&e](double a){ return a*e; }); 
    }
  }
  out1 = out1->cut_col(0,m);
  out2 = out2->cut_col(0,m);
  out2 = out2->transpose();
  return {out1, out2};
}

vector<shared_ptr<Matrix>> Matrix::sqrt2(std::vector<double>& sval) const {
  //return us^0.5 s^0.5vt
  const size_t ns = min(nrows(), ncols());
  auto ans = copy()->svd(ns); //Disable SVD truncation
  shared_ptr<Matrix> out1 = ans[0];
  shared_ptr<Matrix> out2 = ans[2]->transpose();
  sval.resize(ns); 
  copy_n(ans[1]->pdata(), ns, sval.begin());
  for (size_t i = 0; i != ans[1]->nrows(); ++i) {
    const double f = sqrt(ans[1]->element(i,0));
    for_each(out1->element_ptr(0,i), out1->element_ptr(0,i+1), [&f](double& val) { val*=f; });
    for_each(out2->element_ptr(0,i), out2->element_ptr(0,i+1), [&f](double& val) { val*=f; });
  }
  out2 = out2->transpose();
  return {out1, out2};
}

vector<double> Matrix::eig_davidson(Matrix& psi, const size_t max_subspace, const size_t max_iter, const double thresh, const bool mute) const {
  //psi is initial guess on input, eigenvectors on exit
  const size_t ndim = psi.nrows();
  assert(nrows() == ncols());
  assert(ndim == nrows());
  const size_t nstate = psi.ncols();
  Davidson<Matrix,Matrix> davidson(nstate, max_subspace);

  vector<bool> conv(nstate, false);
  vector<double> out(nstate, 0.0);

  unique_ptr<double[]> denom(new double[ndim]);
  for (size_t i = 0; i != ndim; ++i) {
    denom[i] = this->element(i,i);
  }

  for (size_t iter = 0; iter != max_iter; ++iter) {
    auto sigma = make_shared<Matrix>(*this * psi);

    vector<shared_ptr<const Matrix>> sigman;
    vector<shared_ptr<const Matrix>> ccn;
    for (size_t i = 0; i != nstate; ++i) {
      if (!conv[i]) {
        sigman.push_back(sigma->cut_col(i,i+1u));
        ccn.push_back(psi.cut_col(i,i+1u));
      }   
      else{
        sigman.push_back(shared_ptr<const Matrix>());
        ccn.push_back(shared_ptr<const Matrix>());
      }   
    }   
    const vector<double> energies = davidson.compute(ccn, sigman);
    // residual
    vector<shared_ptr<Matrix>> errvec = davidson.residual();
    vector<double> errors;
    for (size_t i = 0; i != nstate; ++i) {
      errors.push_back(errvec[i]->rms());
      conv[i] = errors[i] < thresh;
    }   
    if (any_of(conv.begin(), conv.end(), [] (const bool t) { return (!t); })) {
      for (size_t ist = 0; ist != nstate; ++ist) {
        if (conv.at(ist)) continue;
        auto tmp_cc = make_shared<Matrix>(ndim, 1); 
        double* target_array = tmp_cc->pdata();
        double* source_array = errvec[ist]->pdata();
        const double en = energies[ist];
        for (size_t i = 0; i != ndim; ++i) {
          target_array[i] = source_array[i] / min(en - denom[i], -0.1);
        }
        double nrm = tmp_cc->norm();
        double scal = (nrm > 1.0e-15 ? 1.0/nrm : 0.0);
        tmp_cc->scale(scal);
        copy_n(tmp_cc->pdata(), ndim, psi.element_ptr(0, ist));
      }
    }

    if (!mute) {
      if (nstate!=1 && iter) cout << endl;
      for (size_t i = 0; i != nstate; ++i) {
        cout << setw(7) << iter << setw(3) << i << setw(2) << (conv[i] ? "*" : " ")
                                << setw(17) << fixed << setprecision(8) << energies[i] << "   "
                                << setw(10) << scientific << setprecision(2) << errors[i] << "   " << endl;
      }
    }
    copy_n(energies.begin(), energies.size(), out.begin());
    if (all_of(conv.begin(), conv.end(), [] (const bool t) { return (t); })) break;
  }

  vector<shared_ptr<Matrix>> eigenstates = davidson.civec();
  for (size_t i = 0; i != nstate; ++i) {
    copy_n(eigenstates[i]->pdata(), ndim, psi.element_ptr(0, i));
  }
  return out;
}

vector<double> Matrix::geig_sym(const Matrix& norm, Matrix& psi, const double ratio) const {
  const size_t nr = nrows();
  assert(ncols() == nr);
  assert(norm.nrows() == nr);
  assert(norm.ncols() == nr);
  vector<double> out;
  
  //Step1// //half-inverse
  Vec s(nr);
  Matrix rhs(norm);
  dsyevd(nr, rhs, s);

#ifndef NDEBUG
  double pmin = 0.0;
  for (size_t i = 0; i != nr; ++i) {
    if (s(i) >= 0 ) {
      pmin = s(i);
      break; 
    }
  }
  cout << "Condition Number: " << scientific << s(nr-1)/pmin << endl;
#endif

  const double thresh1 = s.element(nr-1) * ratio;
  size_t n2 = 0;
  for (size_t i = 0; i != nr; ++i) {
    double* istart = rhs.element_ptr(0,i);
    double* iend = rhs.element_ptr(0,i+1);
    if (s(i) > thresh1) {
      const double e = 1.0/sqrt(s(i));
      for_each(istart, iend, [&e](double& val) { val*=e; });
    } else{
      ++n2;
    }
  }
  const Matrix m = rhs%(*this)*rhs;

  if (n2 != 0) {
    const size_t n1 = nr - n2;
    //Step2//

#ifndef NDEBUG 
    cout << "Normal H: " << nr << "->" << n1 << endl;
#endif 
    //Normal H//
    Matrix heff = move(*m.get_submatrix(n2,n2,n1,n1));
    Vec eig(n1);
    dsyevd(n1, heff, eig);

    Matrix wfn(nr, n1);
    wfn.copy_block(n2, 0, n1, n1, heff);

    psi = rhs*wfn;
    out.resize(n1);
    copy_n(eig.pdata(), n1, out.begin());
  } else {
#ifndef NDEBUG 
    cout << "Normal H & S: " << nr << endl;
#endif
    //Normal H & S//
    Vec eig(nr);
    dsyevd(nr, m, eig);
    psi = rhs*m;
    out.resize(nr);
    copy_n(eig.pdata(), nr, out.begin());   
  }
  return out;
}

vector<shared_ptr<Matrix>> Matrix::svd(const size_t dim, const double thresh) const {
  const size_t m = nrows();
  const size_t n = ncols();
  const size_t c = (m<n) ? m:n;
  auto s = make_shared<Matrix>(c,1);
  auto u = make_shared<Matrix>(m,c);
  auto vt = make_shared<Matrix>(c,n);
  dgesvd(m, n, *this, *s, *u, *vt);
  
  for (size_t j = 0; j != min(m,n); ++j) {
    if (s->element(j,0) < thresh) {
      //Zero out
      s->element(j,0) = 0.0;
      for (size_t i = 0; i != m; ++i)
        u->element(i,j) = 0.0;
      for (size_t k = 0; k != n; ++k)
        vt->element(j,k) = 0.0;
      //Truncate
      if (j > dim) {
        s = s->cut_row(0,j);
        u = u->cut_col(0,j);
        vt = vt->cut_row(0,j);
        break;
      }
    }
  }
  return {u,s,vt};
}

/*vector<shared_ptr<Matrix>> Matrix::svdd() const {
  const size_t m = nrows();
  const size_t n = ncols();
  const size_t c = (m<n) ? m:n;
  auto s = make_shared<Matrix>(c,1);
  auto u = make_shared<Matrix>(m,c);
  auto vt = make_shared<Matrix>(c,n);
  dgesdd(m, n, *this, *s, *u, *vt);
  return {u,s,vt};
}*/

vector<shared_ptr<Matrix>> Matrix::jsvd(const size_t dim, const double thresh) const {
  const size_t m = nrows();
  const size_t n = ncols();
  shared_ptr<Matrix> s, u, vt;
  if (m > n) {
    s = make_shared<Matrix>(n,1);
    u = make_shared<Matrix>(m,n);
    auto v = make_shared<Matrix>(n,n);
    dgejsv(m, n, *this, *s, *u, *v);
    vt = v->transpose();
  } else {
    shared_ptr<Matrix> tmp = transpose();
    s = make_shared<Matrix>(m,1);
    u = make_shared<Matrix>(m,m);
    auto v = make_shared<Matrix>(n,m);
    dgejsv(n, m, *tmp, *s, *v, *u);
    vt = v->transpose();
  }
  // Truncate small sigular values
  if (min(m,n) > dim) {
    for (size_t i = dim; i != min(m,n); ++i) {
      if (s->element(i,0) < thresh) {
        s = s->cut_row(0,i);
        u = u->cut_col(0,i);
        vt = vt->cut_row(0,i);
        break;
      }
    }
  }
  return {u,s,vt};
}

vector<shared_ptr<Matrix>> Matrix::svdj(const size_t dim, const double thresh) const {
  const size_t m = nrows();
  const size_t n = ncols();
  shared_ptr<Matrix> s, u, vt; 
  if (m > n) {
    s = make_shared<Matrix>(n,1);
    u = copy();
    auto v = make_shared<Matrix>(n,n);
    dgesvj(m, n, *u, *s, *v);
    vt = v->transpose();
  } else {
    s = make_shared<Matrix>(m,1);
    u = make_shared<Matrix>(m,m);
    auto v = transpose();
    dgesvj(n, m, *v, *s, *u);
    vt = v->transpose();
  }
  // Truncate small sigular values
  if (min(m,n) > dim) {
    for (size_t i = dim; i != min(m,n); ++i) {
      if (s->element(i,0) < thresh) {
        s = s->cut_row(0,i);
        u = u->cut_col(0,i);
        vt = vt->cut_row(0,i);
        break;
      }   
    }   
  }
  return {u,s,vt};
}


void Matrix::qr() {
  const size_t m = nrows();
  const size_t n = ncols();
  dorgqr(m, n, *this);
}


shared_ptr<Matrix> Matrix::qr_full() {
  const size_t m = nrows();
  const size_t n = ncols();
  auto out = make_shared<Matrix>(m,m);
  out->unit();
  dormqr(m, n, *this, *out);
  return out;
}

shared_ptr<Matrix> Matrix::solve(const Matrix& o) const {
  assert(nrows() == o.nrows());

  Matrix tmp(o);
  auto out = this->copy();
  dgesv(nrows(), ncols(), tmp, *out);

  return out;
}

shared_ptr<Matrix> Matrix::solve(const Matrix& o, const int nmin) const {
  assert(nrows() == o.nrows());

  Matrix tmp(o);
  auto out = this->copy();
  unique_ptr<int[]> ipiv(new int[nmin]);
  int info;
  int n = out->nrows();
  int m = out->ncols();
  dgesv_(&nmin, &m, tmp.pdata(), &n, ipiv.get(), out->pdata(), &n, &info);
  if (info != 0) {
    std::cout << "info: " << info <<std::endl; 
    throw std::runtime_error("dgesv_ failed");
  }

  return out;
}


ostream& operator << (ostream& os, const Matrix& o) {
  if (o.pdata()) {
    for (size_t i = 0; i < o.nrows(); ++i) {
      for (size_t j = 0; j < o.ncols(); ++j) {
        os << setw(15) << fixed << setprecision(8) << o.element(i,j) << " ";
        }
      os << "\n";
    }
  } else{
    os << "It's a Null Matrix" << "\n";  
  }
  return os;
}
