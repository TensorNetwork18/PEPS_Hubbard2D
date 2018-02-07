#ifndef _GAUGE_H_
#define _GAUGE_H_

#include <vector>
#include <memory>
#include "peps.h"

class Gauge {
  private:
    std::shared_ptr<PEPS> psi_;
    std::vector<double> svalue_;
    std::vector<double> last_svalue_;
    std::vector<Qtensor> gaugeL_;
    std::vector<Qtensor> gaugeR_;
    std::vector<Qtensor> gaugeU_;
    std::vector<Qtensor> gaugeD_;
    bool conv_;
    double diff_;
  public:
    virtual ~Gauge()=default;
    Gauge()=default;
    Gauge(const std::shared_ptr<PEPS>&);
    Gauge(const Gauge&)=delete;
    void operator = (const Gauge&)=delete;

    bool iteration(const double thresh = 1.0E-5, const size_t max_iter = 1000);
    double compute_LR(size_t, size_t);
    double compute_UD(size_t, size_t);
    void update();
    double check_convergence();
    std::shared_ptr<PEPS> output() const;
};
#endif
