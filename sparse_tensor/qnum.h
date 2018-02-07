#ifndef QNUM_H_
#define QNUM_H_

#include <iostream>


class Qnum {
  public:
    enum Parity { EVEN = 0, ODD = 1 };
  private:
    int u1_;
    Parity parity_; 
  public:
    virtual ~Qnum()=default;
    Qnum(const int u1=0, const Parity parity=EVEN);
    Qnum(const int u1, const bool fermion);

    Qnum(const Qnum&);
    Qnum(Qnum&&);
    Qnum& operator = (const Qnum&);
    Qnum& operator = (Qnum&&);
    bool operator == (const Qnum& o) const { return (u1() == o.u1()) && (parity() == o.parity()); }    
    bool operator != (const Qnum& o) const { return (u1() != o.u1()) || (parity() != o.parity());}
    bool operator < (const Qnum& o) const { return (u1() < o.u1()) || (u1() == o.u1() && parity() < o.parity()); }
    Qnum& operator += (const Qnum& o) { u1_ += o.u1(); parity_ = (Parity)(parity() ^ o.parity()); return *this; }
    Qnum operator + (const Qnum& o) const { Qnum out(*this); out+=o; return out; }
    Qnum& operator -= (const Qnum& o) { u1_ -= o.u1(); parity_ = (Parity)(parity() ^ o.parity()); return *this; }
    Qnum operator - (const Qnum& o) const { Qnum out(*this); out-=o; return out; }
    Qnum operator - () const { Qnum out; out.u1_ = -u1(); out.parity_ = parity(); return out; }

    int u1() const { return u1_; }
    Parity parity() const { return parity_; }
    
    void assign(int, Parity);
    
    static Qnum zero() { Qnum out(0, EVEN); return out; }
    friend std::ostream& operator << (std::ostream& os, const Qnum&);

    void broadcast(int source);
};
#endif
