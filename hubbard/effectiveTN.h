#ifndef _EFF_TN_
#define _EFF_TN_

#include <utility>
#include "../peps/mps2.h"
#include "../peps/tmatrix2.h"
#include "etensor_pool.h"
#include "block.h"

namespace Hubbard2 {
class EffectiveTN {
  private:
    struct pairTM {
      double coeff;
      std::shared_ptr<TMatrix2> left;
      std::shared_ptr<TMatrix2> right;
      virtual ~pairTM()=default;
      pairTM()=default;
      pairTM(const double f, const CanonicalMPS& u, const std::vector<std::shared_ptr<const Etensor>>& m, const CanonicalMPS& d) {
        assert(u.mps_ptr->type() == MPS2::U);
        assert(d.mps_ptr->type() == MPS2::D);
        coeff = f * u.coeff * d.coeff;
        left = make_shared<TMatrix2>(TMatrix2::L, u.mps_ptr, d.mps_ptr, m);
        right = make_shared<TMatrix2>(TMatrix2::R, u.mps_ptr, d.mps_ptr, m);
      }

      Qtensor evaluate(const size_t j) const {
        const Qtensor& utensor = left->uptr()->site(j);
        const Qtensor& dtensor = left->dptr()->site(j);
        const Qtensor& ctensor = left->pmpo()[j]->op();
        const Qtensor& ltensor = left->evaluate(j-1);
        const Qtensor& rtensor = right->evaluate(j+1);
        //Contract order: [(l * u) * c] * (d * r)
        Qtensor tmp_l;
        tmp_l.contract(1.0, ltensor, {0}, utensor, {0});
        tmp_l.contract(1.0, tmp_l, {}, ctensor, {});
 
        Qtensor tmp_r;
        tmp_r.contract(1.0, dtensor, {3}, rtensor, {3});

        Qtensor out;
        out.contract(coeff, tmp_l, {2,5}, tmp_r, {0,3});
        assert(out.charge() == Qnum::zero());
        return out; 
      }

      Qtensor evaluate2(const size_t j) const {
        const Qtensor& utensor = left->uptr()->site(j);
        const Qtensor& dtensor = left->dptr()->site(j);
        const Qtensor& ltensor = left->evaluate(j-1);
        const Qtensor& rtensor = right->evaluate(j+1);
        //Contract order: (l * u) * (d * r)
        Qtensor tmp_l;
        tmp_l.contract(1.0, ltensor, {0}, utensor, {0});
 
        Qtensor tmp_r;
        tmp_r.contract(1.0, dtensor, {3}, rtensor, {3});

        Qtensor out;
        out.contract(coeff, tmp_l, {2,5}, tmp_r, {0,3});
        assert(out.charge() == Qnum::zero());
        return out; 
      }

    };
    vector<pairTM> ham_global_;
    vector<pairTM> ham_local_;
    vector<pairTM> norm_;

  public:
    virtual ~EffectiveTN()=default;
    EffectiveTN()=default;
    EffectiveTN(const double, const shared_ptr<Block>&, const shared_ptr<Block>&,  const EtensorPool&, const int);
    
    std::pair<Qtensor, Qtensor> evaluate(size_t j);
    Qtensor evaluate2(size_t j);
};
}
#endif
