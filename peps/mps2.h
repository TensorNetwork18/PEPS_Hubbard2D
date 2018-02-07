#ifndef _MPS2_
#define _MPS2_

#include <memory>
#include "../sparse_tensor/qtensor.h"
#include "etensor.h"

class MPS2 : public std::enable_shared_from_this<MPS2> {
  public:
    enum Type { U, D };
  private:
    Type type_;
    bool normalized_ = false;
    size_t center_;
    std::vector<Qtensor> data_;
  public:
    virtual ~MPS2()=default;

    MPS2()=default;
    MPS2(Type);
    MPS2(Type, const std::vector<std::shared_ptr<const Etensor>>&, bool n=false);

    MPS2(const MPS2&);
    MPS2(MPS2&&);
    MPS2& operator = (const MPS2&);
    MPS2& operator = (MPS2&&);

    std::shared_ptr<MPS2> copy() const { return std::make_shared<MPS2>(*this); }
    std::shared_ptr<const MPS2> share() const { return shared_from_this(); }
    std::shared_ptr<MPS2> share() { return shared_from_this(); }

    size_t length() const { return data_.size(); }
    Type type() const { return type_; }
    bool normalized() const { return normalized_; }
    size_t center() const { return center_; }
    const std::vector<Qtensor>& data() const { return data_; }
    std::vector<Qtensor>& data() { return data_; }

    const Qtensor& site(const size_t i) const { assert(i < length()); return data_[i]; }
    const Qtensor& operator () (const size_t i) const { return site(i); }
    Qtensor& site(const size_t i) { assert(i < length()); return data_[i]; }
    Qtensor& operator () (const size_t i) { return site(i); }

    void assign(Type o1, std::vector<Qtensor>&& o2) { 
      type_ = o1;
      normalized_ = true;
      center_ = o2.size() -1;
      data_ = std::move(o2); 
    }

    void flip_type();
    double left_normalize(const bool sgn = true);
    //double right_normalize(const bool sgn = true);
    void mixed_normalize(const size_t, const bool sgn = true);
    std::shared_ptr<MPS2> product(const std::vector<std::shared_ptr<const Etensor>>&) const;
    std::shared_ptr<MPS2> compress(const std::map<Qnum,size_t>&, const double thresh = 1.0E-7, const size_t max_iter = 100, const bool mute = true);

    friend std::ostream& operator << (std::ostream& os, const MPS2& o);
    void broadcast(int source);
  private:
    static bool sgn_permute(const std::vector<size_t>& fidx);
};
#endif
