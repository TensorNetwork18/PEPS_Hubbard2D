#ifndef _MKL_EXTERN_H
#define _MKL_EXTERN_H

extern "C" {
  void dscal_(const int* n, const double* a, double* x, const int* incx);

  void mkl_domatcopy_ (const char* ordering, const char* trans, const int* rows, const int* cols, 
    const double* alpha, const double* a, const int* lda, double* B, const int* ldb);

  void dlapmt_(const bool*, const int*, const int*, double*, const int*, const int*);

  void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);

  double ddot_(const int* n, const double* x, const int* incx, const double* y, const int* incy);

  void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
    const double* alpha, const double* a, const int* lda, const double* b, const int* ldb,
    const double* beta, double* c, const int* ldc);

  void dsyevd_(const char* jobz, const char* uplo, const int* n, double* a,
    const int* lda, double* w, double* work, const int* lwork, 
    int* iwork, const int* liwork, int* info );

  void dgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n,
    double* a, const int* lda, double* s, double* u, const int* ldu,
    double* vt, const int* ldvt, double* work, const int* lwork, int* info);

  void dgesdd_(const char* jobz, const int* m, const int* n, 
    double* a, const int* lda, double* s, double* u, const int* ldu, 
    double* vt, const int* ldvt, double* work, const int* lwork, int* iwork, int* info);

  void dgejsv_(const char* joba, const char* jobu, const char* jobv, const char* jobr, const char* jobt, const char* jobp,
    const int* m, const int* n, double* a, const int* lda, double* sva, double* u, const int* ldu, double* v, const int* ldv,
    double* work, const int* lwork, int* iwork, int* info);

  void dgesvj_(const char* joba, const char* jobu, const char* jobv, 
    const int* m, const int* n, double* a, const int* lda, double* sva, const int* mv, double* v, const int* ldv,
    double* work, const int* lwork, int* info);

  void dgeqrf_(const int* m, const int* n, double* a, const int* lda, double* tau, double* work, const int* lwork, int* info);
  
  void dorgqr_(const int* m, const int* n, const int* k, double* a, const int* lda, double* tau, double* work, const int* lwork, int* info);

  void dormqr_(const char* side, const char* trans, const int* m, const int* n, const int* k, double* a, const int* lda, double* tau,
    double* c, const int* ldc, double* work, const int* lwork, int* info);

  void dsysv_(const char* uplo, const int* n, const int* nrhs, double* a, const int* lda,
    int* ipiv, double* b, const int* ldb, double* work, int* lwork, int* info );

  void dsygv_(const int* itype, const char* jobz, const char* uplo, const int* n, 
    double* a, const int* lda, double *b, const int* ldb, double *w, double* work, int* lwork, int* info);

  void dgebal_(const char* job, const int* n, double* a, const int* lda, int* ilo, int* ihi, double* scale, int* info);

  void dgebak_(const char* job, const char* side, const int* n, const int* ilo, const int* ihi, double* scale, const int* m, double* v, const int* ldv, int* info);

  void dtrtri_(const char* uplo, const char* diag, const int* n, double* a, const int* lda, int* info);

  void dpotrf_(const char* uplo, const int* n, double* a, const int* lda, int* info);
  
  void dpocon_(const char* uplo, const int* n, const double* a, const int* lda, const double* anorm, double* rcond, double* work, int* iwork, int* info);
  
  void dgesv_(const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* b, const int* ldb, int* info);
  //void dpoequ_(const int* n, double* a, const int* lda, double* s, double* scond, double* amax, int* info);

  /*void dtrevc_(const char* side, const char* howmny, const bool* select, const int* n, const int* t, const int* ldt, 
    double* vl, const int* ldvl, double* vr, const int* ldvr, const int* mm, int* m, double* work, int* info);

  void dsygvx_(const int* itype, const char* jobz, const char* range, const char* uplo, 
    const int* n, double* a,  const int* lda, double* b, const int* ldb, const double* vl, const double* vu, const int* il, const int* iu, const double* abstol, 
    int* m, double* w, double* z, const int* ldz, double* work, const int* lwork, int* iwork, int* ifail, int* info);*/
}

#endif /* _MKL_EXTERN_H  */
