#include "../timer.h"
#include "gauge.h"

using namespace std;

Gauge::Gauge(const std::shared_ptr<PEPS>& psi) : psi_(psi) {
  
}

bool Gauge::iteration(const double thresh, const size_t max_iter) {
  Timer time_gauge;
  const size_t nr = psi_->nrows();
  const size_t nc = psi_->ncols();
  conv_ = false;
  for (size_t iter = 0; iter != max_iter; ++iter) {
    if (!svalue_.empty())
      svalue_.clear();
    if (!gaugeL_.empty())
      gaugeL_.clear();
    if (!gaugeR_.empty())
      gaugeR_.clear();
    if (!gaugeU_.empty())
      gaugeU_.clear();
    if (!gaugeD_.empty())
      gaugeD_.clear();
    for (size_t i = 0; i != nr; ++i) {
      for (size_t j = 0; j != nc-1u; ++j) {
        svalue_.push_back(compute_LR(i,j));
      }
    }
    for (size_t i = 0; i != nr-1u; ++i) {
      for (size_t j = 0; j != nc; ++j) {
        svalue_.push_back(compute_UD(i,j));
      }
    }
    update();
    conv_ = (check_convergence() < thresh);
    if (conv_)
      break;
  }
  cout << "  Time for Gauge-Fixing:    " << time_gauge.tick() << endl; 
  return conv_;
}

double Gauge::compute_LR (const size_t i, const size_t j) {
  // (i,j) X^-1  X     Y     Y^-1 (i,j+1)
  // (i,j) riR[1] riR[0] riL[0] riL[1] (i,j+1)
  const size_t nr = psi_->nrows();
  const size_t nc = psi_->ncols();
  const bool c1 = (i == 0 || i == nr);
  const bool c2 = (j == 0 || j == nc);
  const bool c3 = (j+1 == 0 || j+1 == nc);
  int eta_r = 4, eta_l = 4;
  if (c1 && c2) {
    eta_r = 2;
  } else if (c1 || c2) {
    eta_r = 3;
  }
  const auto& riR = psi_->rdmR(i,j).split(eta_r);
  if (c1 && c3) {
    eta_l = 2;
  } else if (c1 || c3) {
    eta_l = 3;
  }
  const auto& riL = psi_->rdmL(i,j+1).split(eta_l);
  Qtensor tmp;
  tmp.contract(1.0, riR[0], {0u}, riL[0], {0u});
  vector<double> sv;
  const auto& tmp2 = tmp.sqrt2(sv);
  Qtensor gR;
  gR.contract(1.0, riR[1], {0}, tmp2[0], {0});
  gaugeR_.push_back(move(gR));
  Qtensor gL;
  gL.contract(1.0, tmp2[1], {1}, riL[1], {0});
  gaugeL_.push_back(move(gL));
  double sum = 0.0;
  assert(sv.size() != 0); 
  for (auto&& e : sv)
    sum += e;
//cout << "i, j " << i << " " << j << " :: " << sum/static_cast<double>(sv.size()) << endl;
  return sum/static_cast<double>(sv.size());
}

double Gauge::compute_UD(const size_t i, const size_t j) {
  const size_t nr = psi_->nrows();
  const size_t nc = psi_->ncols();
  const bool c1 = (i == 0 || i == nr);
  const bool c2 = (j == 0 || j == nc);
  const bool c3 = (i+1 == 0 || i+1 == nr);
  int eta_d = 4, eta_u = 4;
  if (c1 && c2) {
    eta_d = 2;
  } else if (c1 || c2) {
    eta_d = 3;
  }
  const auto& riD = psi_->rdmD(i,j).split(eta_d);
  if (c2 && c3) {
    eta_u = 2;
  } else if (c2 || c3) {
    eta_u = 3;
  }
  const auto& riU = psi_->rdmU(i+1,j).split(eta_u);
  Qtensor tmp;
  tmp.contract(1.0, riD[0], {0u}, riU[0], {0u});
  vector<double> sv;
  const auto& tmp2 = tmp.sqrt2(sv);
  Qtensor gD;
  gD.contract(1.0, riD[1], {0}, tmp2[0], {0});
  gaugeD_.push_back(move(gD));
  Qtensor gU;
  gU.contract(1.0, tmp2[1], {1}, riU[1], {0});
  gaugeU_.push_back(move(gU));
  double sum = 0.0;
  assert(sv.size() != 0);
  for (auto&& e : sv)
    sum += e;
//cout << "i, j " << i << " " << j << " :: " << sum/static_cast<double>(sv.size()) << endl;
  return sum/static_cast<double>(sv.size());
}

void Gauge::update() {
  const size_t nr = psi_->nrows();
  const size_t nc = psi_->ncols();
  for (size_t i = nr; i-- > 0;) {
    for (size_t j = nc-1u; j-- > 0;) {
      Qtensor tmp;
      tmp.contract(1.0, psi_->site(i,j), {4}, gaugeR_.back(), {0});
      psi_->site(i,j) = move(tmp);
      gaugeR_.pop_back();
      tmp.contract(1.0, gaugeL_.back(), {1}, psi_->site(i,j+1), {0});
      psi_->site(i,j+1) = move(tmp);  
      gaugeL_.pop_back();
    }
  }
  for (size_t i = nr-1u; i-- > 0;) {
    for (size_t j = nc; j-- > 0;) {
      Qtensor tmp;
      tmp.contract(1.0, psi_->site(i,j), {3}, gaugeD_.back(), {0});
      psi_->site(i,j) = tmp.permute({0,1,2,4,3});
      gaugeD_.pop_back();
      tmp.contract(1.0, gaugeU_.back(), {1}, psi_->site(i+1,j), {1});
      psi_->site(i+1,j) = tmp.permute({1,0,2,3,4});
      gaugeU_.pop_back();
    }
  }
}

double Gauge::check_convergence() {
  const size_t ns = svalue_.size();
  assert(ns != 0);
  if (last_svalue_.size() == 0)
    last_svalue_.resize(ns, 0.0);
  double diff = 0.0;
  for (size_t i = 0; i != ns; ++i) {
    diff += abs(last_svalue_[i] - svalue_[i])/svalue_[i];  
  }
  last_svalue_ = move(svalue_);
  svalue_.clear();
  diff_ = diff/static_cast<double>(ns);
//cout << diff_ << endl;
  return diff_;
}

shared_ptr<PEPS> Gauge::output() const {
  if (conv_) {
    cout << "  *** PEPS fix-point reached ***" << endl;
    cout << "  Diff = " << setw(15) << fixed << setprecision(8) << diff_ << endl << endl;
    return psi_;
  } else{
    cout << "  *** Failure Convergence PEPS Gauge-fixing ***" << endl;
    cout << "  Diff = " << setw(15) << fixed << setprecision(8) << diff_ << endl << endl;
    return psi_;
  }
}
