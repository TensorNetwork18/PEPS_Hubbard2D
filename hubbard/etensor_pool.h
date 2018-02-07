#ifndef _ETENSOR_POOL_H_
#define _ETENSOR_POOL_H_

#include <cassert>
#include <memory>
#include <array>
#include "../peps/peps.h"
#include "../peps/etensor.h"

namespace Hubbard2 {
class EtensorPool {
  private:
    std::shared_ptr<const PEPS> psi_;
    std::array<Qtensor, 6> op_; //alphaC, alphaA, betaC, betaA, doublon, identity
    std::vector<std::shared_ptr<Etensor>> alphaC_;
    std::vector<std::shared_ptr<Etensor>> alphaA_;
    std::vector<std::shared_ptr<Etensor>> betaC_;
    std::vector<std::shared_ptr<Etensor>> betaA_;
    std::vector<std::shared_ptr<Etensor>> doublon_;
    std::vector<std::shared_ptr<Etensor>> identity_;

  public:
    virtual ~EtensorPool()=default;
    EtensorPool()=default;
    EtensorPool(const std::shared_ptr<const PEPS>&);

    const std::shared_ptr<const PEPS> psi() const { return psi_; }
    void update(const size_t, const size_t) const;
    const std::shared_ptr<Etensor>& alphaC(const size_t i, const size_t j) const { assert(i < psi_->nrows() && j < psi_->ncols()); return alphaC_[i*psi_->ncols()+j]; }
    const std::shared_ptr<Etensor>& alphaA(const size_t i, const size_t j) const { assert(i < psi_->nrows() && j < psi_->ncols()); return alphaA_[i*psi_->ncols()+j]; }
    const std::shared_ptr<Etensor>& betaC(const size_t i, const size_t j) const { assert(i < psi_->nrows() && j < psi_->ncols()); return betaC_[i*psi_->ncols()+j]; }
    const std::shared_ptr<Etensor>& betaA(const size_t i, const size_t j) const { assert(i < psi_->nrows() && j < psi_->ncols()); return betaA_[i*psi_->ncols()+j]; }
    const std::shared_ptr<Etensor>& doublon(const size_t i, const size_t j) const { assert(i < psi_->nrows() && j < psi_->ncols()); return doublon_[i*psi_->ncols()+j]; }
    const std::shared_ptr<Etensor>& identity(const size_t i, const size_t j) const { assert(i < psi_->nrows() && j < psi_->ncols()); return identity_[i*psi_->ncols()+j]; }

    std::vector<std::shared_ptr<const Etensor>> get_row_alphaC(const size_t, const size_t) const;
    std::vector<std::shared_ptr<const Etensor>> get_row_alphaA(const size_t, const size_t) const;
    std::vector<std::shared_ptr<const Etensor>> get_row_betaC(const size_t, const size_t) const;
    std::vector<std::shared_ptr<const Etensor>> get_row_betaA(const size_t, const size_t) const;
    std::vector<std::shared_ptr<const Etensor>> get_row_doublon(const size_t, const size_t) const;
    std::vector<std::shared_ptr<const Etensor>> get_row_identity(const size_t) const;
    std::vector<std::shared_ptr<const Etensor>> get_row_alphaC_alphaA(const size_t, const size_t, const size_t) const;
    std::vector<std::shared_ptr<const Etensor>> get_row_alphaA_alphaC(const size_t, const size_t, const size_t) const;
    std::vector<std::shared_ptr<const Etensor>> get_row_betaC_betaA(const size_t, const size_t, const size_t) const;
    std::vector<std::shared_ptr<const Etensor>> get_row_betaA_betaC(const size_t, const size_t, const size_t) const;
    
  private:
    void init();
};
}
#endif
