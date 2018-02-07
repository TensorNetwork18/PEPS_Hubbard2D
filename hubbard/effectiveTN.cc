#include "effectiveTN.h"

using namespace std;

Hubbard2::EffectiveTN::EffectiveTN (const double u_term, const shared_ptr<Block>& uBlock, const shared_ptr<Block>& dBlock,  const EtensorPool& etensor_pool, const int irow) {
  assert(uBlock->type() == MPS2::U);
  assert(dBlock->type() == MPS2::D);
  const CanonicalMPS& identityU = uBlock->identity();
  const CanonicalMPS& identityD = dBlock->identity();
  const vector<shared_ptr<const Etensor>>& identityM = etensor_pool.get_row_identity(irow);

  for (auto&& element : uBlock->blocks()) {
    switch(element.type) {
      //-1.0 * AlphaC * AlphaA * I
      case CanonicalMPS::AlphaC : {
        const size_t jcol = element.position;
        ham_global_.emplace_back(pairTM(-1.0, element, etensor_pool.get_row_alphaA(irow, jcol), identityD));
        break;
      }
      // 1.0 * AlphaA * AlphaC * I 
      case CanonicalMPS::AlphaA : {
        const size_t jcol = element.position;
        ham_global_.emplace_back(pairTM( 1.0, element, etensor_pool.get_row_alphaC(irow, jcol), identityD));
        break;
      }
      //-1.0 * BetaC * BetaA * I
      case CanonicalMPS::BetaC : {
        const size_t jcol = element.position;
        ham_global_.emplace_back(pairTM(-1.0, element, etensor_pool.get_row_betaA(irow, jcol), identityD));
        break;
      }
      // 1.0 * BetaA * BetaC * I    
      case CanonicalMPS::BetaA : {
        const size_t jcol = element.position;
        ham_global_.emplace_back(pairTM( 1.0, element, etensor_pool.get_row_betaC(irow, jcol), identityD));
        break;
      }
      // 1.0 * Ham * I * I
      case CanonicalMPS::Ham : {
        if (element.exist())
          ham_local_.emplace_back(pairTM( 1.0, element, identityM, identityD));
        break;
      }
      default : {
        continue;
      }
    }
  }

  for (auto&& element : dBlock->blocks()) {
    switch(element.type) {
      //-1.0 * I * AlphaC * AlphaA
      case CanonicalMPS::AlphaA : {
        const size_t jcol = element.position;
        ham_global_.emplace_back(pairTM(-1.0, identityU, etensor_pool.get_row_alphaC(irow, jcol), element));
        break;
      }
      // 1.0 * I * AlphaA * AlphaC 
      case CanonicalMPS::AlphaC : {
        const size_t jcol = element.position;
        ham_global_.emplace_back(pairTM( 1.0, identityU, etensor_pool.get_row_alphaA(irow, jcol), element));
        break;
      }
      //-1.0 * I * BetaC * BetaA
      case CanonicalMPS::BetaA : {
        const size_t jcol = element.position;
        ham_global_.emplace_back(pairTM(-1.0, identityU, etensor_pool.get_row_betaC(irow, jcol), element));
        break;
      }
      // 1.0 * I * BetaA * BetaC   
      case CanonicalMPS::BetaC : {
        const size_t jcol = element.position;
        ham_global_.emplace_back(pairTM( 1.0, identityU, etensor_pool.get_row_betaA(irow, jcol), element));
        break;
      }
      // 1.0 * I * I * Ham
      case CanonicalMPS::Ham : {
        if (element.exist())
          ham_local_.emplace_back(pairTM( 1.0,  identityU, identityM, element));
        break;
      }
      default : {
        continue;
      }
    }
  }

  const size_t nc = etensor_pool.psi()->ncols();
  for (size_t jcol = 0; jcol != nc-1; ++jcol) {
    //-1.0 * I * AlphaC_AlphaA * I
    ham_global_.emplace_back(pairTM(-1.0, identityU, etensor_pool.get_row_alphaC_alphaA(irow, jcol, jcol+1), identityD));
    // 1.0 * I * AlphaA_AlphaC * I
    ham_global_.emplace_back(pairTM( 1.0, identityU, etensor_pool.get_row_alphaA_alphaC(irow, jcol, jcol+1), identityD));
    //-1.0 * I * BetaC_BetaA * I
    ham_global_.emplace_back(pairTM(-1.0, identityU, etensor_pool.get_row_betaC_betaA(irow, jcol, jcol+1), identityD));
    // 1.0 * I * BetaA_BetaC * I
    ham_global_.emplace_back(pairTM( 1.0, identityU, etensor_pool.get_row_betaA_betaC(irow, jcol, jcol+1), identityD));
  }
  for (size_t jcol = 0; jcol != nc; ++jcol) {
    //   u * I * Doublon * I
    ham_global_.emplace_back(pairTM(u_term, identityU, etensor_pool.get_row_doublon(irow, jcol), identityD));
  }

  // 1.0 * I * I * I
  norm_.emplace_back(pairTM( 1.0, identityU, identityM, identityD));
}

pair<Qtensor, Qtensor> Hubbard2::EffectiveTN::evaluate(size_t j) {
  assert(norm_.size() == 1);
  Qtensor out_N(norm_[0].evaluate(j));

  Qtensor out_H(out_N.index(), true);
  const bool do_parallel = (ham_global_.size() != 1 && mpi__->size() != 1) ? true : false;
  
  if (do_parallel) { 
    const int njobs = ham_global_.size();
    const int nproc = mpi__->world_size();
    const int iproc = mpi__->world_rank();

    // distribute //
    const int ntask = (iproc < njobs%nproc) ? njobs / nproc + 1 : njobs / nproc;
    const int ltask = ham_local_.size();
    for (int itask = 0; itask != ntask; ++itask) {
      const int job_id = itask * nproc + iproc;
      const auto& job = ham_global_[job_id];
      out_H += job.evaluate(j);
    }

    for (int itask = 0; itask != ltask; ++itask) {
      const auto& job = ham_local_[itask];
      out_H += job.evaluate(j);
    }
    // allreduce
    out_H.allreduce_pack();

  } else {

    for (const auto& job : ham_global_) {
      out_H += job.evaluate(j);
    }
    
    for (const auto& job : ham_local_) {
      out_H += job.evaluate(j);
    }
  }

  return make_pair(move(out_H), move(out_N));
}

Qtensor Hubbard2::EffectiveTN::evaluate2(size_t j) {
  assert(norm_.size() == 1);
  return norm_[0].evaluate2(j);
}
