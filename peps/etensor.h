#ifndef _ETENSOR_
#define _ETENSOR_

#include "../sparse_tensor/qtensor.h"

class Etensor { //NonCopyable
  private:
    const Qtensor& wfnRef_;
    Qtensor op_;
    Qtensor conj_;
    Qtensor data_;
  public:
    virtual ~Etensor()=default;
    Etensor()=default;
    Etensor(const Qtensor&);
    Etensor(const Qtensor&, const Qtensor&);
    Etensor(const Etensor&)=delete;
    Etensor(Etensor&&)=delete;
    void operator = (const Etensor&)=delete;
    void operator = (Etensor&&)=delete;

    const Qtensor& wfn() const { return wfnRef_; }
    const Qtensor& op() const { return op_; }
    const Qtensor& conj() const { return conj_; }
    const Qtensor& data() const { return data_; }
    Qtensor& data() { return data_; }
    void update();
};
#endif
