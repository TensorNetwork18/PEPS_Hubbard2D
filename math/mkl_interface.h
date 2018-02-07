#ifndef __MKL_INTERFACE_H
#define __MKL_INTERFACE_H

#include <memory>
#include <stdexcept>
#include <algorithm>
#include "mkl_extern.h"
#include "tensor_base.h"

using namespace std;

template<class T>
void dscal(const double& a, T& out) {
  const int n = static_cast<int>(out.size());
  const int incx = 1;
  dscal_(&n, &a, out.pdata(), &incx);
}

template<class T>
shared_ptr<T> transpose_mkl(const int rows, const int cols, const T& inp) {
  auto out = make_shared<T>(cols, rows);
  const double alpha = 1.0;
  mkl_domatcopy_("C","T", &rows, &cols, &alpha, inp.pdata(), &rows, out->pdata(), &cols);
  return out;
}

template<class U>
void dlapmt(const bool f, const int m, const int n, const int* k, U& out) {
  assert(static_cast<int>(out.size()) == m*n);
  dlapmt_(&f, &m, &n, out.pdata(), &m, k);
}

template<class T, class U>
void daxpy(const double& a, const T& inp, U& out) {
  const int n = inp.size();
  assert(static_cast<int>(out.size()) == n);
  const int incx = 1;
  const int incy = 1;
  daxpy_(&n, &a, inp.pdata(), &incx, out.pdata(), &incy);    
}


template<class T1, class T2>
double ddot(const T1& inp1, const T2& inp2) {
  const int n = inp1.size();
  assert(static_cast<int>(inp2.size()) == n);
  const int incx = 1;
  const int incy = 1;
  return ddot_(&n, inp1.pdata(), &incx, inp2.pdata(), &incy);
}

template<class T1, class T2, class U>
void dgemm(const char* transa, const char* transb, const int& m, const int& n, const int& k, const double& alpha, const T1& inp1, const T2& inp2, const double& beta, U& out) {
  assert(static_cast<int>(inp1.size()) == m*k);
  assert(static_cast<int>(inp2.size()) == k*n);
  if (transa[0] == 'N' && transb[0] == 'N') {
    dgemm_(transa, transb, &m, &n, &k, &alpha, inp1.pdata(), &m, inp2.pdata(), &k, &beta, out.pdata(), &m);
  } else if(transa[0] == 'T' && transb[0] == 'N') {
    dgemm_(transa, transb, &m, &n, &k, &alpha, inp1.pdata(), &k, inp2.pdata(), &k, &beta, out.pdata(), &m);
  } else {
    throw std::logic_error("this case is not implemented in dgemm");
  }
}

template<class T1, class T2>
void dsyevd(const int n, T1& out1, T2& out2) {
  assert(static_cast<int>(out1.size()) == n*n);
  int lwork = -1;
  int liwork = -1;
  int iwkopt, info;
  double wkopt;
  std::unique_ptr<double[]> w(new double[n]);
  dsyevd_("V", "U", &n, out1.pdata(), &n, w.get(), &wkopt, &lwork, &iwkopt, &liwork, &info);
  lwork = static_cast<int>(wkopt);
  std::unique_ptr<double[]> work(new double[lwork]);
  liwork = iwkopt;
  std::unique_ptr<int[]> iwork(new int[liwork]);
  dsyevd_("V", "U", &n, out1.pdata(), &n, w.get(), work.get(), &lwork, iwork.get(), &liwork, &info);
  if (info != 0) throw std::runtime_error("dsyevd_ failed");
  const int m = (n < static_cast<int>(out2.size())) ? n : out2.size();
  std::copy_n(w.get(), m, out2.pdata());
}

template<class T, class U1, class U2>
void dgesvd(const int m, const int n, const T& inp, U1& out1, U2& out2, U2& out3) {
  assert(static_cast<int>(inp.size()) == m*n);
  int lwork = -1;
  int info;
  double wkopt;
  if (m > n) {
    dgesvd_("S", "A", &m, &n, inp.pdata(), &m, out1.pdata(), out2.pdata(), &m, out3.pdata(), &n, &wkopt, &lwork, &info );
    lwork = static_cast<int>(wkopt);
    std::unique_ptr<double[]> work(new double[lwork]);
    dgesvd_("S", "A", &m, &n, inp.pdata(), &m, out1.pdata(), out2.pdata(), &m, out3.pdata(), &n, work.get(), &lwork, &info );
  } else {
    dgesvd_("A", "S", &m, &n, inp.pdata(), &m, out1.pdata(), out2.pdata(), &m, out3.pdata(), &m, &wkopt, &lwork, &info );
    lwork = static_cast<int>(wkopt);
    std::unique_ptr<double[]> work(new double[lwork]);
    dgesvd_("A", "S", &m, &n, inp.pdata(), &m, out1.pdata(), out2.pdata(), &m, out3.pdata(), &m, work.get(), &lwork, &info );
  }    
  if (info != 0) {
    std::cout << "info: " << info <<std::endl;
    throw std::runtime_error("dgesvd_ failed");
  }
}

template<class T, class U1, class U2>
void dgesdd(const int m, const int n, const T& inp, U1& out1, U2& out2, U2& out3) {
  assert(static_cast<int>(inp.size()) == m*n);
  int lwork = -1;
  int info;
  double wkopt;
  if (m > n) {
    std::unique_ptr<int[]> iwork(new int[8*n]);
    dgesdd_("S", &m, &n, inp.pdata(), &m, out1.pdata(), out2.pdata(), &m, out3.pdata(), &n, &wkopt, &lwork, iwork.get(), &info );
    lwork = static_cast<int>(wkopt);
    std::unique_ptr<double[]> work(new double[lwork]);
    dgesdd_("S", &m, &n, inp.pdata(), &m, out1.pdata(), out2.pdata(), &m, out3.pdata(), &n, work.get(), &lwork, iwork.get(), &info );
  } else {
    std::unique_ptr<int[]> iwork(new int[8*m]);
    dgesdd_("S", &m, &n, inp.pdata(), &m, out1.pdata(), out2.pdata(), &m, out3.pdata(), &m, &wkopt, &lwork, iwork.get(), &info );
    lwork = static_cast<int>(wkopt);
    std::unique_ptr<double[]> work(new double[lwork]);
    dgesdd_("S", &m, &n, inp.pdata(), &m, out1.pdata(), out2.pdata(), &m, out3.pdata(), &m, work.get(), &lwork, iwork.get(), &info );
  }    
  if (info != 0) {
    std::cout << "info: " << info <<std::endl;
    throw std::runtime_error("dgesdd_ failed");
  }
}

template<class T, class U1, class U2>
void dgejsv(const int m, const int n, const T& inp, U1& out1, U2& out2, U2& out3) {
  assert(static_cast<int>(inp.size()) == m*n);
  assert(m >= n);
  int info;
  int lwork = std::max(2*m+n, 6*n+2*n*n);
  std::unique_ptr<double[]> work(new double[std::max(7, lwork)]);
  std::unique_ptr<int[]> iwork(new int[std::max(3, m+3*n)]);
  dgejsv_("F", "U", "V", "R", "N", "P", &m, &n, inp.pdata(), &m, out1.pdata(), out2.pdata(), &m, out3.pdata(), &n, work.get(), &lwork, iwork.get(), &info);
  if (info != 0) {
    std::cout << "info: " << info <<std::endl;
    throw std::runtime_error("dgejsv_ failed");
  }
}

template<class T, class U1, class U2> 
void dgesvj(const int m, const int n, T& inp, U1& out1, U2& out2) {
  assert(static_cast<int>(inp.size()) == m*n);
  assert(m >= n); 
  int info;
  int lwork = std::max(6, m+n);
  std::unique_ptr<double[]> work(new double[lwork]);
  dgesvj_("G", "U", "V", &m, &n, inp.pdata(), &m, out1.pdata(), &n, out2.pdata(), &n, work.get(), &lwork, &info);
  if (info != 0) {
    std::cout << "info: " << info <<std::endl;
    throw std::runtime_error("dgesvj_ failed");
  }
}

/*template<class U1, class U2>
void dgeqrf(const int m, const int n, U1& out1, U2& out2) {
  assert((int)out1.size() == m*n);
  int lwork = -1;
  int info;
  double wkopt;
  dgeqrf_(&m, &n, out1.pdata(), &m, out2.pdata(), &wkopt, &lwork, &info);
  lwork = (int)wkopt;
  std::unique_ptr<double[]> work(new double[lwork]);
  dgeqrf_(&m, &n, out1.pdata(), &m, out2.pdata(), work.get(), &lwork, &info); 
  if (info != 0) throw std::runtime_error("dgeqrf_ failed");
}*/

template<class U1>
void dorgqr(const int m, const int n, U1& out1) {
  assert(static_cast<int>(out1.size()) == m*n);
  assert(m >= n);
  int lwork = -1; 
  int info;
  double wkopt;
  std::unique_ptr<double[]> tau(new double[n]);
  dgeqrf_(&m, &n, out1.pdata(), &m, tau.get(), &wkopt, &lwork, &info);
  lwork = static_cast<int>(wkopt);
  std::unique_ptr<double[]> work(new double[lwork]);
  dgeqrf_(&m, &n, out1.pdata(), &m, tau.get(), work.get(), &lwork, &info); 
  if (info != 0) throw std::runtime_error("dgerqf_ failed");
  
  lwork = -1;
  dorgqr_(&m, &n, &n, out1.pdata(), &m, tau.get(), &wkopt, &lwork, &info);
  lwork = static_cast<int>(wkopt);
  std::unique_ptr<double[]> work2(new double[lwork]);
  dorgqr_(&m, &n, &n, out1.pdata(), &m, tau.get(), work2.get(), &lwork, &info);
  if (info != 0) throw std::runtime_error("dorgqr_ failed");
}

template<class U1, class U2>
void dormqr(const int m, const int n, U1& out1, U2& out2) {
  assert(static_cast<int>(out1.size()) == m*n);
  assert(m >= n);
  int lwork = -1; 
  int info;
  double wkopt;
  std::unique_ptr<double[]> tau(new double[n]);
  dgeqrf_(&m, &n, out1.pdata(), &m, tau.get(), &wkopt, &lwork, &info);
  lwork = static_cast<int>(wkopt);
  std::unique_ptr<double[]> work(new double[lwork]);
  dgeqrf_(&m, &n, out1.pdata(), &m, tau.get(), work.get(), &lwork, &info); 
  if (info != 0) throw std::runtime_error("dgerqf_ failed");

  lwork = -1; 
  dormqr_("L", "N", &m, &m, &n, out1.pdata(), &m, tau.get(), out2.pdata(), &m, &wkopt, &lwork, &info);
  lwork = static_cast<int>(wkopt);
  std::unique_ptr<double[]> work2(new double[lwork]);
  dormqr_("L", "N", &m, &m, &n, out1.pdata(), &m, tau.get(), out2.pdata(), &m, work2.get(), &lwork, &info);
  if (info != 0) throw std::runtime_error("dormqr_ failed");
}

template<class U1, class U2>
void dsysv(const int n, const int nrhs, U1& out1, U2& out2) {
  assert(out1.size() == n*n);
  assert(out2.size() == n*nrhs);
  int lwork = -1;
  int info;
  double wkopt;
  std::unique_ptr<int[]> ipiv(new int[n]);
  dsysv_("L", &n, &nrhs, out1.pdata(), &n, ipiv.get(), out2.pdata(), &n, &wkopt, &lwork, &info );
  lwork = static_cast<int>(wkopt);
  std::unique_ptr<double[]> work(new double[lwork]);
  dsysv_("L", &n, &nrhs, out1.pdata(), &n, ipiv.get(), out2.pdata(), &n, work.get(), &lwork, &info );
  if (info != 0) {
    std::cout << "info: " << info <<std::endl; 
    throw std::runtime_error("dsysv_ failed");
  }
}

template<class T, class U>
void dgesv(const int n, const int m, const T& lhs, U& out) {
  assert(static_cast<int>(lhs.size()) == n*n);
  assert(static_cast<int>(out.size()) == n*m);
  int info;
  std::unique_ptr<int[]> ipiv(new int[n]);
  dgesv_(&n, &m, lhs.pdata(), &n, ipiv.get(), out.pdata(), &n, &info);
  if (info != 0) {
    std::cout << "info: " << info <<std::endl; 
    throw std::runtime_error("dgesv_ failed");
  }
}

template<class U1, class U2>
void dgebal(const int n, U1& out1, U2& out2) {
  assert(static_cast<int>(out1.size()) == n*n);
  assert(static_cast<int>(out2.size()) == n);
  int ilo, ihi, info;
  dgebal_("S", &n, out1.pdata(), &n, &ilo, &ihi, out2.pdata(), &info);
  if (info != 0) throw std::runtime_error("dgebal_ failed");
}

template<class T, class U>
void dgebak(const char* side, const int row, const int col, const T& inp, U& out) {
  assert(static_cast<int>(inp.size()) == row);
  assert(static_cast<int>(out.size()) == row*col);
  int ilo = 1;
  int info;
  dgebak_("S", side, &row, &ilo, &row, inp.pdata(), &col, out.pdata(), &row, &info);
  if (info != 0) throw std::runtime_error("dgebak_ failed");
}

template<class U1, class U2, class U3>
void dsygv(const int itype, const int n, U1& out1, U2& out2, U3& out3) { // out1->H, out2->N, out3->E
  assert(static_cast<int>(out1.size()) == n*n);
  assert(static_cast<int>(out2.size()) == n*n);
  assert(static_cast<int>(out3.size()) == n);
  int lwork = -1;
  int info;
  double wkopt;
  dsygv_(&itype, "V", "U", &n, out1.pdata(), &n, out2.pdata(), &n, out3.pdata(), &wkopt, &lwork, &info); 
  lwork = static_cast<int>(wkopt);
  std::unique_ptr<double[]> work(new double[lwork]);
  dsygv_(&itype, "V", "U", &n, out1.pdata(), &n, out2.pdata(), &n, out3.pdata(), work.get(), &lwork, &info);
  if (info != 0) {
    std::cout << "info: " << info <<std::endl; 
    throw std::runtime_error("dsygv_ failed");
  }
}

template<class U>
void dtrtri(const char* uplo, const int n, U& out) {
  assert(out.size() == n*n);
  int info;
  dtrtri_(uplo, "N", &n, out.pdata(), &n, &info);
  if (info != 0) throw std::runtime_error("dtrtri_ failed");
}

template<class U>
double dpocon(const char* uplo, const int n, const U& out) { //out: symmetric positive matrix & return Cholesky factorization
  assert(static_cast<int>(out.size()) == n*n);
  
  //compute 1-norm
  double norm1 = 0.0;
  for (int i = 0; i != n; ++i) {
    double sum = 0.0;
    std::for_each(out.pdata()+i*n, out.pdata()+(i+1)*n, [&sum](const double& val) { sum += abs(val); });
    norm1 = std::max(sum, norm1);
  }
  
  //Cholesky factorization
  int info;
  dpotrf_(uplo, &n, out.pdata(), &n, &info);
  if (info != 0) {
    std::cout << "info: " << info <<std::endl;
    throw std::runtime_error("dpotrf_ failed");
  }
  
  //Compute condition number
  double rcond;
  std::unique_ptr<double[]> work(new double[3*n]);
  std::unique_ptr<int[]> iwork(new int[n]);
  dpocon_(uplo, &n, out.pdata(), &n, &norm1, &rcond, work.get(), iwork.get(), &info);
  if (info != 0) {
    std::cout << "info: " << info <<std::endl;
    throw std::runtime_error("dpocon_ failed");
  }
  return rcond;
}
#endif
