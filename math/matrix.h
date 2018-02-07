#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "tensor_base.h"
#include "vec.h"

class Vec;
class Matrix: public Tensor_base<double> {
  public:
    virtual ~Matrix()=default;
    
    Matrix() : Tensor_base<double>() { }
    Matrix(const size_t m, const size_t n) : Tensor_base<double>(m,n) { }

    Matrix(const Matrix& o) : Tensor_base<double>(o) { }
    Matrix(Matrix&& o) : Tensor_base<double>(std::move(o)) { }
    Matrix& operator = (const Matrix& o) { Tensor_base<double>::operator=(o); return *this; }
    Matrix& operator = (Matrix&& o) { Tensor_base<double>::operator=(o); return *this; }

    size_t size() const { return nrows()*ncols(); }
    size_t nrows() const { return this->extent(0); }
    size_t ncols() const { return this->extent(1); }

    std::shared_ptr<Matrix> clone() const { return std::make_shared<Matrix>(nrows(), ncols()); }
    std::shared_ptr<Matrix> copy() const { return std::make_shared<Matrix>(*this); }

    void ax_plus_y (const double a, const Matrix& o) { assert(nrows() == o.nrows() && ncols() == o.ncols()); daxpy(a, o, *this); }
    Matrix& operator += (const Matrix& o) { this->ax_plus_y(1.0, o); return *this; }
    Matrix operator + (const Matrix& o) const { Matrix out(*this); out += o; return out; }
    Matrix& operator -= (const Matrix& o) { this->ax_plus_y(-1.0, o); return *this; }
    Matrix operator - (const Matrix& o) const { Matrix out(*this); out -= o; return out; }
    Matrix& operator *= (const Matrix& o) { *this = *this * o; return *this; }
    Matrix operator * (const Matrix&) const;
    Matrix operator % (const Matrix&) const;
    bool operator == (const Matrix&)const;

    void copy_block(const size_t, const size_t, const size_t, const size_t, const Matrix&);
    void copy_block(const size_t, const size_t, const size_t, const size_t, const double* a);
    std::shared_ptr<Matrix> resize(const size_t, const size_t);
    std::shared_ptr<Matrix> get_submatrix(const size_t, const size_t, const size_t, const size_t) const;
    std::shared_ptr<Matrix> cut_row(const size_t, const size_t) const;
    std::shared_ptr<Matrix> cut_col(const size_t, const size_t) const;

    void symmetrize();
    void add_diag(const double, const size_t, const size_t);
    void unit();
    bool is_symmetric(const double thresh = 1.0e-11) const;
    bool is_antisymmetric(const double thresh = 1.0e-11) const;
    bool is_identity(const double thresh = 1.0e-11) const; 

    double dot_product(const Matrix& o) const;
    double norm() const;
    double variance() const;
    double rms() const;
    double sum_diag() const;

    void scale(const double);
    std::shared_ptr<Matrix> transpose() const;
    std::shared_ptr<Matrix> tildex(const double thresh = 1.0e-13) const; // psudo half_inverse
    std::vector<std::shared_ptr<Matrix>> split(const int a = 2, const double thresh = 1.0e-13) const;
    std::vector<std::shared_ptr<Matrix>> sqrt2(std::vector<double>&) const;
    void eig_sym(Vec&);
    std::vector<std::shared_ptr<Matrix>> svd(const size_t dim = 1, const double thresh = 1.0e-13) const;
    std::vector<std::shared_ptr<Matrix>> jsvd(const size_t dim = 1, const double thresh = 1.0e-13) const;
    std::vector<std::shared_ptr<Matrix>> svdj(const size_t dim = 1, const double thresh = 1.0e-13) const;
    void qr();
    std::shared_ptr<Matrix> qr_full();
    std::shared_ptr<Matrix> solve(const Matrix&) const;
    std::shared_ptr<Matrix> solve(const Matrix&, const int) const;
    std::vector<double> eig_davidson(Matrix& psi, const size_t max_subspace = 2, const size_t max_iter = 100, const double thresh = 1.0e-11, const bool mute = true) const;
    std::vector<double> geig_sym(const Matrix& norm, Matrix& psi, const double ratio = 1.0e-11) const;
    friend std::ostream& operator << (std::ostream& os, const Matrix&);
};

#endif
