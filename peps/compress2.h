#ifndef _COMPRESS2_H_
#define _COMPRESS2_H_

#include "mps2.h"
#include "tmatrix2.h"

class Compress2 { //Only Compress A Normalized MPS2
  private:
    std::shared_ptr<const MPS2> ref_;
    std::map<Qnum, size_t> qlist_;
    bool mute_ = true;    

    std::shared_ptr<MPS2> data_;
    std::unique_ptr<TMatrix2> left_;
    std::unique_ptr<TMatrix2> right_;
    std::vector<double> sweep_err_;
    double last_average_err_;
    bool conv_;
  public:
    virtual ~Compress2()=default;
    Compress2()=default;
    Compress2(const std::shared_ptr<const MPS2>&, const std::vector<Qnum>&, const std::vector<size_t>&, const bool m = true);  
    Compress2(const std::shared_ptr<const MPS2>&, const std::map<Qnum, size_t>&, const bool m = true);

    Compress2(const Compress2&)=delete;
    void operator = (const Compress2&)=delete;

    const std::map<Qnum, size_t>& qlist() const { return qlist_; }
    std::map<Qnum, size_t>& qlist() { return qlist_; }
    

    bool iteration(const double thresh, const size_t max_iter);
    void left_to_right_iter(const size_t);
    void right_to_left_iter(const size_t);
    double update_one_site(const size_t);
    double check_convergence();
    std::shared_ptr<MPS2> output() const;
  private:
    void init();
    void truncate_copy();
};
#endif
