#include "etensor.h"

using namespace std;


Etensor::Etensor(const Qtensor& o) : wfnRef_(o) {
  update();
}

Etensor::Etensor(const Qtensor& o1, const Qtensor& o2) : wfnRef_(o1), op_(o2) {
  update();
}

void Etensor::update() {
  conj_ = wfnRef_.permute({4,3,2,1,0},false);
  for (auto&& k : conj_.index()) {
    k = k.reverse();
  }
  if (op_.rank() != 0) {
    conj_.contract(1.0, conj_, {2}, op_, {0});
    data_.contract(1.0, conj_, {4}, wfnRef_, {2});
  } else {
    conj_ = conj_.permute({0,1,3,4,2});
    data_.contract(1.0, conj_, {4}, wfnRef_, {2});
  }
  data_ = data_.permute({3,4,2,5,1,6,0,7}); //{l'lu'ud'dr'r}
}
