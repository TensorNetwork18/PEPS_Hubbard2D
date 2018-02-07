#ifndef DAVIDSON_H_
#define DAVIDSON_H_

#include <vector>
#include <memory>
#include <cassert>
#include "vec.h"
#include "matrix.h"


template <typename T, typename U>
class Davidson {
  protected:
    struct BasisPair {
      public:
        std::shared_ptr<const T> cc;
        std::shared_ptr<const U> sigma;
        BasisPair() { }
        BasisPair(std::shared_ptr<const T> a, std::shared_ptr<const U> b) : cc(a), sigma(b) { }
    };

    size_t nstate_;
    size_t max_;
    size_t size_;

    std::vector<std::shared_ptr<BasisPair>> basis_;

    // Hamiltonian
    std::shared_ptr<Matrix> mat_;
    // eivenvalues
    Vec vec_;
    // an eigenvector
    std::shared_ptr<Matrix> eig_;
    // overlap matrix
    std::shared_ptr<Matrix> overlap_;
  public:
    // Davidson with periodic collapse of the subspace
    Davidson()=default;
    Davidson(size_t n, size_t max) : nstate_(n), max_((max+1)*n), size_(0), vec_(max_) {
      if (max < 2) throw std::runtime_error("Davidson diagonalization requires at least two trial vectors per root.");
    }

    double compute(std::shared_ptr<const T> cc, std::shared_ptr<const U> cs) {
      assert(nstate_ == 1);
      return compute(std::vector<std::shared_ptr<const T>>{cc},
                     std::vector<std::shared_ptr<const U>>{cs}).front();
    }

    std::vector<double> compute(std::vector<std::shared_ptr<const T>> cc, std::vector<std::shared_ptr<const U>> cs) {
      // reset the convergence flags
      std::vector<bool> converged(nstate_, false);

      // new pairs
      std::vector<std::shared_ptr<BasisPair>> newbasis;
      assert(cc.size() == nstate_ && cs.size() == nstate_);
      for (size_t ic = 0; ic < nstate_; ++ic) {
        assert(!cc[ic] == !cs[ic]);
        if (cc[ic] && cs[ic]) {
          newbasis.push_back(std::make_shared<BasisPair>(cc[ic], cs[ic]));
        } else{
          converged[ic] = true;
        }
      }

      // adding new matrix elements
      {
        const size_t n = newbasis.size();
        mat_ = mat_ ? mat_->resize(size_+n, size_+n) : std::make_shared<Matrix>(n, n);
        overlap_ = overlap_ ? overlap_->resize(size_+n, size_+n) : std::make_shared<Matrix>(n, n);
      }
      basis_.insert(basis_.end(), newbasis.begin(), newbasis.end());
      for (auto& ib : newbasis) {
        ++size_;
        size_t i = 0;
        for (auto& b : basis_) {
          if (i+1 > size_)
            break;
          mat_->element(i, size_-1) = b->cc->dot_product(*(ib->sigma));
          mat_->element(size_-1, i) = mat_->element(i, size_-1);

          overlap_->element(i, size_-1) = b->cc->dot_product(*(ib->cc));
          overlap_->element(size_-1, i) = overlap_->element(i, size_-1);
          ++i;
        }
      }
      
      // canonical orthogonalization
      std::shared_ptr<const Matrix> ovlp_scr = overlap_->tildex();
      if (ovlp_scr->ncols() < nstate_)
        throw std::runtime_error("Too much linear dependency in guess vectors provided to DavidsonDiag; cannot obtain the requested number of states.");
      
      // diagonalize matrix to get
      eig_ = std::make_shared<Matrix>(*ovlp_scr % *mat_ * *ovlp_scr);
      eig_->eig_sym(vec_);
      eig_ = std::make_shared<Matrix>(*ovlp_scr * *eig_->cut_col(0,nstate_));
      
      // first basis vector is always the current best guess
      std::vector<std::shared_ptr<T>> cv = civec();
      std::vector<std::shared_ptr<U>> sv = sigmavec();
      for (size_t i = 0; i != nstate_; ++i) {
        basis_[i] = std::make_shared<BasisPair>(cv[i], sv[i]);
      }

      // due to this, we need to transform mat_ and overlap_
      auto trans = eig_->resize(eig_->nrows(), eig_->nrows());
      for (size_t i = nstate_; i != eig_->nrows(); ++i)
        trans->element(i,i) = 1.0;
      mat_ = std::make_shared<Matrix>(*trans % *mat_ * *trans);
      overlap_ = std::make_shared<Matrix>(*trans % *overlap_ * *trans);

      eig_->fill(0.0);
      for (size_t i = 0; i != nstate_; ++i) {
        eig_->element(i, i) = 1.0;
      }
      
      // possibly reduce the dimension
      assert(size_ == basis_.size());
      if (size_ > max_ - nstate_) {
        std::map<size_t, size_t> remove;
        const size_t soff = size_ - newbasis.size();
        for (size_t i = 0; i != nstate_; ++i) {
          if (converged[i]) continue;
          // a vector with largest weight will be removed.
          size_t n = 0;
          double abs = 1.0e10;
          for (size_t j = nstate_; j < soff; ++j) {
            if (std::abs(trans->element(j, i)) < abs) {
              if (remove.find(j) != remove.end()) continue;
              abs = std::abs(trans->element(j, i));
              n = j;
            }
          }
          remove.emplace(n, soff+remove.size());
        }
        assert(newbasis.size() == remove.size());
        //std::cout << "    ** throwing out " << remove.size() << " trial vectors **" << std::endl;
        for (auto m : remove) {
          basis_[m.first] = basis_[m.second];
          mat_->copy_block(0, m.first, size_, 1, *mat_->get_submatrix(0, m.second, size_, 1));
          mat_->copy_block(m.first, 0, 1, size_, *mat_->get_submatrix(m.second, 0, 1, size_));
          overlap_->copy_block(0, m.first, size_, 1, *overlap_->get_submatrix(0, m.second, size_, 1));
          overlap_->copy_block(m.first, 0, 1, size_, *overlap_->get_submatrix(m.second, 0, 1, size_));

          trans->copy_block(m.first, 0, 1, nstate_, *trans->get_submatrix(m.second, 0, 1, nstate_));
        }
        basis_ = std::vector<std::shared_ptr<BasisPair>>(basis_.begin(), basis_.end()-remove.size());
        size_ = basis_.size();
        mat_ = mat_->get_submatrix(0, 0, size_, size_);
        overlap_ = overlap_->get_submatrix(0, 0, size_, size_);
      }
      vector<double> out(nstate_);
      std::copy_n(vec_.pdata(), nstate_, out.begin());
      return out;
    }

    // perhaps can be cleaner.
    std::vector<std::shared_ptr<U>> residual() {
      std::vector<std::shared_ptr<U>> out;
      for (size_t i = 0; i != nstate_; ++i) {
        auto tmp = basis_.front()->sigma->clone();
        size_t k = 0;
        for (auto& iv : basis_) {
          if (std::abs(eig_->element(k++,i)) > 1.0e-16)
            tmp->ax_plus_y(-vec_(i)*eig_->element(k-1,i), *(iv->cc));
        }
        k = 0;
        for (auto& iv : basis_) {
          if(std::abs(eig_->element(k++,i)) > 1.0e-16)
            tmp->ax_plus_y(eig_->element(k-1,i), *(iv->sigma));
        }
        out.push_back(tmp);
      }
      return out;
    }

    // returns ci vector
    std::vector<std::shared_ptr<T>> civec() {
      std::vector<std::shared_ptr<T>> out;
      for (size_t i = 0; i != nstate_; ++i) {
        auto tmp = basis_.front()->cc->clone();
        size_t k = 0;
        for (auto& iv : basis_) {
          tmp->ax_plus_y(eig_->element(k++,i), *(iv->cc));
        }
        out.push_back(tmp);
      }
      return out;
    }

    // return sigma vector
    std::vector<std::shared_ptr<U>> sigmavec() {
      std::vector<std::shared_ptr<U>> out;
      for (size_t i=0; i!=nstate_; ++i) {
        auto tmp = basis_.front()->sigma->clone();
        size_t k = 0;
        for (auto& iv : basis_) {
          tmp->ax_plus_y(eig_->element(k++,i), *(iv->sigma));
        }
        out.push_back(tmp);
      }
      return out;
    }
};
#endif
