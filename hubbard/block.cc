#include "../mpi/mpi_interface.h"
#include "block.h"

Hubbard2::Block::Block(const double u_term, const std::shared_ptr<Block>& last_block_ptr, const EtensorPool& etensor_pool, const size_t irow, const std::map<Qnum, size_t>& compress_map) {

  if (last_block_ptr) {
    assert(irow < etensor_pool.psi()->nrows());
    type_ = last_block_ptr->type();

    vector< tuple<CanonicalMPS::Type, double, CanonicalMPS, vector<shared_ptr<const Etensor>>, int> > job_global;
    vector< tuple<CanonicalMPS::Type, double, CanonicalMPS, vector<shared_ptr<const Etensor>>, int> > job_local;
    for (const auto& element : last_block_ptr->blocks()) {
      switch(element.type) {
        case CanonicalMPS::AlphaC : {
          //AlphaC * AlphaA -> Ham
          const double fac = (type() == MPS2::U) ? -1.0 : 1.0;
          const size_t jcol = element.position;
          job_global.emplace_back(CanonicalMPS::Ham, fac, element, etensor_pool.get_row_alphaA(irow,jcol), -1);
          break;
        }

        case CanonicalMPS::AlphaA : {
          //AlphaA * AlphaC -> Ham
          const double fac = (type() == MPS2::U) ? 1.0 : -1.0;
          const size_t jcol = element.position;
          job_global.emplace_back(CanonicalMPS::Ham, fac, element, etensor_pool.get_row_alphaC(irow, jcol), -1);
          break;
        }

        case CanonicalMPS::BetaC : {
          //BetaC * BetaA -> Ham
          const double fac = (type() == MPS2::U) ? -1.0 : 1.0;
          const size_t jcol = element.position;
          job_global.emplace_back(CanonicalMPS::Ham, fac, element, etensor_pool.get_row_betaA(irow, jcol), -1);
          break;
        }

        case CanonicalMPS::BetaA : {
          //BetaA * BetaC -> Ham
          const double fac = (type() == MPS2::U) ? 1.0 : -1.0;
          const size_t jcol = element.position;
          job_global.emplace_back(CanonicalMPS::Ham, fac, element, etensor_pool.get_row_betaC(irow, jcol), -1);
          break;
        }

        case CanonicalMPS::Identity : {
          //Identity * Identity -> Identity
          job_global.emplace_back(CanonicalMPS::Identity, 1.0, element, etensor_pool.get_row_identity(irow), -1);
          const int nc = etensor_pool.psi()->ncols();
          for (int jcol = 0; jcol != nc; ++jcol) {
            //Identity * AlphaC -> AlphaC
            job_global.emplace_back(CanonicalMPS::AlphaC,   1.0, element, etensor_pool.get_row_alphaC(irow, jcol), jcol);
            //Identity * AlphaA -> AlphaA
            job_global.emplace_back(CanonicalMPS::AlphaA,   1.0, element, etensor_pool.get_row_alphaA(irow, jcol), jcol);
            //Identity * BetaC -> BetaC
            job_global.emplace_back(CanonicalMPS::BetaC,    1.0, element, etensor_pool.get_row_betaC(irow, jcol), jcol);
            //Identity * BetaA -> BetaA
            job_global.emplace_back(CanonicalMPS::BetaA,    1.0, element, etensor_pool.get_row_betaA(irow, jcol), jcol);
            //u * Identity * Doublon -> Ham
            job_global.emplace_back(CanonicalMPS::Ham, u_term, element, etensor_pool.get_row_doublon(irow, jcol), -1);
          }
          for (int jcol = 0; jcol != nc-1; ++jcol) {
            //-1.0 * Identity * AlphaCAlphaA -> Ham
            job_global.emplace_back(CanonicalMPS::Ham, -1.0, element, etensor_pool.get_row_alphaC_alphaA(irow, jcol, jcol+1), -1);
            // 1.0 * Identity * AlphaAAlphaC -> Ham
            job_global.emplace_back(CanonicalMPS::Ham,  1.0, element, etensor_pool.get_row_alphaA_alphaC(irow, jcol, jcol+1), -1);
            //-1.0 * Identity * BetaCBetaA -> Ham
            job_global.emplace_back(CanonicalMPS::Ham, -1.0, element, etensor_pool.get_row_betaC_betaA(irow, jcol, jcol+1), -1);
            // 1.0 * Identity * BetaABetaC -> Ham
            job_global.emplace_back(CanonicalMPS::Ham,  1.0, element, etensor_pool.get_row_betaA_betaC(irow, jcol, jcol+1), -1);
          }
          break;
        }

        case CanonicalMPS::Ham : {
          //Ham * Identity -> Ham
          if (element.exist())
            job_local.emplace_back(CanonicalMPS::Ham, 1.0, element, etensor_pool.get_row_identity(irow), -1);
          break;
        }

        default:
          throw logic_error("Fail upon Hubbard2::Block::Block");
      }
    }

    //Do the work
    const bool compress = (compress_map.size() != 0);
    const bool do_parallel = (job_global.size() != 1 && mpi__->size() != 1) ? true : false;

    if (do_parallel) {
      const int njobs = job_global.size();
      blocks_.resize(njobs + job_local.size());

      const int nproc = mpi__->world_size();
      const int iproc = mpi__->world_rank();
      
      // distribute //
      const int ntask = (iproc < njobs%nproc) ? njobs / nproc + 1 : njobs / nproc;
      const int ltask = job_local.size();
      for (int itask = 0; itask != ntask; ++itask) {
        const int job_id = itask * nproc + iproc;
        const auto& job = job_global[job_id];
        auto new_mps_ptr = get<2>(job).mps_ptr -> product(get<3>(job));
        auto new_coeff = get<1>(job) * get<2>(job).coeff * new_mps_ptr->left_normalize();
        if (compress)
          new_mps_ptr = new_mps_ptr->compress(compress_map); 
        blocks_[job_id] = CanonicalMPS(get<0>(job), new_coeff, new_mps_ptr, get<4>(job));
      }

      for (int itask = 0; itask != ltask; ++itask) {
        const int job_id = njobs + itask;
        const auto& job = job_local[itask];
        auto new_mps_ptr = get<2>(job).mps_ptr -> product(get<3>(job));
        auto new_coeff = get<1>(job) * get<2>(job).coeff * new_mps_ptr->left_normalize();
        if (compress)
          new_mps_ptr = new_mps_ptr->compress(compress_map);
        blocks_[job_id] = CanonicalMPS(get<0>(job), new_coeff, new_mps_ptr, get<4>(job));
      }
     
      // broadcast //
      for (int job_id = 0; job_id != njobs; ++job_id) {
        blocks_[job_id].broadcast(job_id%nproc);
      }

    } else {

      for (const auto& job : job_global) {
        auto new_mps_ptr = get<2>(job).mps_ptr -> product(get<3>(job));
        auto new_coeff = get<1>(job) * get<2>(job).coeff * new_mps_ptr->left_normalize();
        if (compress)
          new_mps_ptr = new_mps_ptr->compress(compress_map);
        blocks_.emplace_back(CanonicalMPS(get<0>(job), new_coeff, new_mps_ptr, get<4>(job)));
      }

      for (const auto& job : job_local) {
        auto new_mps_ptr = get<2>(job).mps_ptr -> product(get<3>(job));
        auto new_coeff = get<1>(job) * get<2>(job).coeff * new_mps_ptr->left_normalize();
        if (compress)
          new_mps_ptr = new_mps_ptr->compress(compress_map);
        blocks_.emplace_back(CanonicalMPS(get<0>(job), new_coeff, new_mps_ptr, get<4>(job)));
      }
    }
    
  } else {
    throw logic_error("Hubbard2::Block::Block failed");
  }
}

Hubbard2::Block::Block(MPS2::Type type, const EtensorPool& etensor_pool, const size_t irow) : type_(type) {
  assert(irow == 0u || irow == etensor_pool.psi()->nrows()-1u);
  blocks_.emplace_back(CanonicalMPS(CanonicalMPS::Identity, 1.0, make_shared<MPS2>(type_, etensor_pool.get_row_identity(irow)), -1));
}
