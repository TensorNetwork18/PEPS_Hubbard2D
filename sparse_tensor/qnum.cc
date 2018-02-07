#include <cassert>
#include "../mpi/mpi_interface.h"
#include "qnum.h"


using namespace std;


Qnum::Qnum(const int u1, const Parity parity) : u1_(u1), parity_(parity) {
}

Qnum::Qnum(const int u1, const bool fermion) : u1_(u1) {
  if (fermion && u1%2 != 0) {
    parity_ = ODD;
  } else{
    parity_ = EVEN;
  }
}

Qnum::Qnum(const Qnum& o) : u1_(o.u1()), parity_(o.parity()) {
}

Qnum::Qnum(Qnum&& o) : u1_(o.u1()), parity_(o.parity()) {
  o.u1_ = 0;
  o.parity_ = EVEN;
}

Qnum& Qnum::operator = (const Qnum& o) {
  if (this != &o) {
    u1_ = o.u1();
    parity_ = o.parity();
  }
  return *this;
}

Qnum& Qnum::operator = (Qnum&& o) {
  if (this != &o) {
    u1_ = o.u1();
    parity_ = o.parity();
    o.u1_ = 0;
    o.parity_ = EVEN;
  }
  return *this;
}

void Qnum::assign(int u1, Parity parity) {
  u1_ = u1;
  parity_ = parity;
}

ostream& operator << (ostream& os, const Qnum& o) {
  os << "Qnum: U1 = " << o.u1() << ", Parity = ";
  if (o.parity() == Qnum::EVEN) {
    os << "EVEN ";
  } else{
    os << "ODD  ";
  }
  return os;
}

void Qnum::broadcast(int source) {
  if (mpi__->size() != 1) {
    mpi__->broadcast(&u1_, 1u, source);
    if (mpi__->rank() != source)
      parity_ = (u1_%2 == 0) ? EVEN : ODD;
  }
  return;
}
