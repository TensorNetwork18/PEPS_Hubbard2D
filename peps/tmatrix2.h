#ifndef _TMATRIX2_H_
#define _TMATRIX2_H_

#include <vector>
#include <memory>
#include "mps2.h"
#include "etensor.h"


class TMatrix2 { //NonCopyable
  public:
    enum Type{ L, R };
  private:
    Type type_;
    std::shared_ptr<const MPS2> uptr_;
    std::shared_ptr<const MPS2> dptr_;
    std::vector<std::shared_ptr<const Etensor>> pmpo_;
    bool sgn_ = true;
    std::vector<Qtensor> data_;
    std::shared_ptr<Qtensor> base_;
  public:
    virtual ~TMatrix2()=default;
    TMatrix2()=default;
    TMatrix2(Type, const std::shared_ptr<const MPS2>&, const std::shared_ptr<const MPS2>&, const std::vector<std::shared_ptr<const Etensor>>&, const bool f=true);
    TMatrix2(Type, const std::shared_ptr<const MPS2>&, const std::vector<std::shared_ptr<const Etensor>>&, const std::shared_ptr<const MPS2>&, const bool f=true);

    TMatrix2(const TMatrix2&)=delete;
    void operator = (const TMatrix2&)=delete;

    size_t length() const { return data_.size(); }
    Type type() const { return type_; }
    bool sgn() const { return sgn_; }
    const std::vector<std::shared_ptr<const Etensor>>& pmpo() const { return pmpo_; }
    const std::shared_ptr<const MPS2>& uptr() const { return uptr_; }
    const std::shared_ptr<const MPS2>& dptr() const { return dptr_; }
    const std::vector<Qtensor>& data() const { return data_; }
    std::vector<Qtensor>& data() { return data_; }

    const Qtensor& site(const size_t i) const { assert(i < length()); return data_[i]; }
    const Qtensor& operator () (const size_t i) const { return site(i); }
    Qtensor& site(const size_t i) { assert(i < length()); return data_[i]; }
    Qtensor& operator () (const size_t i) { return site(i); }
    const Qtensor& base() const { return *base_; }
    Qtensor& base() { return *base_; }

    Qtensor& evaluate(const size_t);
    Qtensor& evaluate_from_left(const size_t);
    Qtensor& evaluate_from_right(const size_t);    
  private:
    static bool sgn_permute(const std::vector<size_t>& fidx);
};
#endif
