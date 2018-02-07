#ifndef _MINIMIZE_E_H_
#define _MINIMIZE_E_H_

#include <utility>
#include <memory>
#include <vector>
#include "../peps/peps.h"
#include "etensor_pool.h"
#include "block.h"
#include "effectiveTN.h"

namespace Hubbard2 {
class MinimizeE {
  private:
    double u_term_;
    std::shared_ptr<PEPS> psi_;
    EtensorPool etensor_pool_;
    
    std::map<Qnum, size_t> compress_map_;

    std::vector<std::shared_ptr<Block>> pBlocksU_;
    std::vector<std::shared_ptr<Block>> pBlocksD_;
    std::shared_ptr<EffectiveTN> pEffTN_;
    
    std::vector<double> sweep_ene_;
    double last_average_ene_;
    bool conv_;

    double noise_ = 0.0;
  public:
    virtual ~MinimizeE()=default;
    MinimizeE()=default;

    MinimizeE(const double o1, const std::shared_ptr<PEPS>& o2) : MinimizeE(o1, o2, {}) { }
    MinimizeE(const double, const std::shared_ptr<PEPS>&, const std::map<Qnum, size_t>&);
    
    MinimizeE(const MinimizeE&)=delete;
    void operator = (const MinimizeE&)=delete;
    void compute(); 
    std::pair<bool, double> iteration(const double perturb = 0.0, const double thresh = 1.0E-5, const size_t max_iter = 1000);
    void up_to_down_iter(const size_t);
    void down_to_up_iter(const size_t);
    void left_to_right_iter(const size_t);
    void right_to_left_iter(const size_t);
    void update_one_site(const size_t, const size_t, const bool update = true);
    double check_convergence();
    std::shared_ptr<PEPS> output() const;
};
}
#endif
