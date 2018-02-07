#include <algorithm>
#include "compress2.h"

using namespace std;


Compress2::Compress2(const shared_ptr<const MPS2>& o1, const vector<Qnum>& qnums, const vector<size_t>& qdims, const bool m) : ref_(o1), mute_(m) {
  const size_t nq = qnums.size();
  assert(qdims.size()==nq);
  if (!qlist_.empty()) qlist_.clear();
  for (size_t i=0; i!=nq; ++i) {
    qlist_.emplace(qnums[i], qdims[i]);
  }
  init();
}

Compress2::Compress2(const shared_ptr<const MPS2>& o1, const map<Qnum, size_t>& o2, const bool m) : ref_(o1), qlist_(o2), mute_(m) {
  init();
}

bool Compress2::iteration(const double thresh, const size_t max_iter) {
  const size_t nl = ref_->length();
  if (!sweep_err_.empty()) {
    sweep_err_.clear();
  }
  last_average_err_ = 0.0;
  conv_ = false;
  for (size_t iter = 0; iter != max_iter; ++iter) {
    right_to_left_iter(nl-1);
    left_to_right_iter(0);
    conv_ = (check_convergence() < thresh);
    if (conv_) break;
  }
  for (auto&& o : data_->data()) {
    o.reverse();
  }
  data_->flip_type();
  return conv_;
}

void Compress2::left_to_right_iter(const size_t isite) {
  const size_t nl = ref_->length();
  for (size_t i = isite; i != nl-1; ++i) {
    sweep_err_.emplace_back(update_one_site(i));
    data_->mixed_normalize(i+1, false);
  }
}

void Compress2::right_to_left_iter(const size_t isite) {
  for (size_t i = isite+1; --i > 0;) {
    sweep_err_.emplace_back(update_one_site(i));
    data_->mixed_normalize(i-1, false);
  }
}

double Compress2::update_one_site(const size_t isite){
  Qtensor err;

  switch (ref_->type()) {
    case MPS2::U: {
      Qtensor out;
      out.contract(1.0, left_->evaluate(isite-1), {0}, ref_->site(isite), {0}, false);
      out.contract(1.0, out, {3}, right_->evaluate(isite+1), {0}, false);
      err.contract(1.0, out, {0,1,2,3}, data_->site(isite), {0,1,2,3}, false);
      out.reverse();
      data_->site(isite) = move(out);
      break;
    }
    case MPS2::D: {
      Qtensor out;
      out.contract(1.0, left_->evaluate(isite-1), {1}, ref_->site(isite), {0}, false);
      out.contract(1.0, out, {3}, right_->evaluate(isite+1), {1}, false);
      err.contract(1.0, out, {0,1,2,3}, data_->site(isite), {0,1,2,3}, false);
      out.reverse();
      data_->site(isite) = move(out);
      break;
    }
    default:
      throw logic_error("Compress2::update_one_site failed");
  }

  assert(err.get_block()->size() == 1);
  const double error = 1.0 - err.get_block()->element(0,0);
  if (!mute_) 
    cout << "site: " << isite << " :: " << setw(15) << fixed << setprecision(8) << error << endl;
  return error;
}

double Compress2::check_convergence() {
  const double average_err = accumulate(sweep_err_.begin(), sweep_err_.end(), 0.0)/(double)sweep_err_.size();
  const double out = abs(last_average_err_ - average_err);
  last_average_err_ = average_err;
  sweep_err_.clear();
  return out;
}

shared_ptr<MPS2> Compress2::output() const {
  if (conv_) {
    if (!mute_) {
      cout << "  *** MPS2 Compress2 Converged ***" << endl;
      cout << "  Average_Err = " << setw(15) << fixed << setprecision(8) << last_average_err_ << endl;
    }
  } else {
    cout << "  *** Failure Convergence MPS2 Compress2  ***" << endl;
    cout << "  Average_Err = " << setw(15) << fixed << setprecision(8) << last_average_err_ << endl;
  }
  return data_;
}

void Compress2::init(){
  const size_t nl = ref_->length();
  assert(ref_->normalized());
  switch (ref_->type()) {
    case MPS2::U: {
      data_ = make_shared<MPS2>(MPS2::D);
      truncate_copy();
      const double s = data_->left_normalize(false);
      if (s < 0)
        data_->site(nl-1).scale(-1.0); 
      left_.reset(new TMatrix2(TMatrix2::L, ref_, data_, {}, false));
      right_.reset(new TMatrix2(TMatrix2::R, ref_, data_, {}, false));
      break;
    }
    case MPS2::D: {
      data_ = make_shared<MPS2>(MPS2::U);
      truncate_copy();
      const double s = data_->left_normalize(false);
      if (s < 0)
        data_->site(nl-1).scale(-1.0);
      left_.reset(new TMatrix2(TMatrix2::L, data_, ref_, {}, false));
      right_.reset(new TMatrix2(TMatrix2::R, data_, ref_, {}, false));
      break;
    }
  }
}

void Compress2::truncate_copy() {
  for (const auto& o : ref_->data()) {
    Qtensor out(-o.charge(), o.index(), false);
    for (auto&& i : out.index()) {
      i = i.reverse();
    }

    auto&& b0 = out.index()[0].basis();
    for (auto it = b0.begin(); it!=b0.end();) {
      auto it2 = qlist_.find(it->first);
      if (it2 != qlist_.end()) {
        it->second = min(it->second, it2->second);
        ++ it;
      } else {
        it = b0.erase(it);
      }
    }

    auto&& b3 = out.index()[3].basis();
    for (auto it = b3.begin(); it!=b3.end();) {
      auto it2 = qlist_.find(it->first);
      if (it2 != qlist_.end()) {
        it->second = min(it->second, it2->second);
        ++it; 
      } else {
        it = b3.erase(it);
      }   
    }

    o.load();
    for (auto it = o.data().begin(); it != o.data().end(); ++it) {
      auto it0 = qlist_.find(it->first[0u]);
      auto it3 = qlist_.find(it->first[3u]);
      if (it0 != qlist_.end() && it3 != qlist_.end()) {
        const Tensor& tensor = *it->second;
        const size_t d0 = min(tensor.extent(0u), it0->second);
        const size_t d1 = tensor.extent(1u);
        const size_t d2 = tensor.extent(2u);
        const size_t d3 = min(tensor.extent(3u), it3->second);
        auto ptr = make_shared<Tensor>(d0,d1,d2,d3);
        ptr->randomize(-1.0E-8, 1.0E-8); //maybe not a good idea
        for (size_t i3 = 0; i3 != d3; ++i3) {
          for (size_t i2 = 0; i2 != d2; ++i2) {
            for (size_t i1 = 0; i1 != d1; ++i1) {
              copy_n(tensor.element_ptr(0u,i1,i2,i3), d0, ptr->element_ptr(0u,i1,i2,i3));
            }
          }
        }
        out.data().emplace(it->first, move(ptr));
      }
    }
    o.dump();
    out.flush();
    data_->data().push_back(move(out));
  }
  return;
}
