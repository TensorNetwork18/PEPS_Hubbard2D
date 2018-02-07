#include <memory>
#include "vec.h"
#include "mkl_interface.h"

using namespace std;


ostream& operator << (std::ostream& os, const Vec& o) {
  if (o.pdata()) {
    for (size_t i = 0; i != o.size(); ++i) {
      os << setw(15) << fixed << setprecision(8) << o.element(i) << " ";
    }
    os << "\n";
  } else {
    os << "It's a Null Vec" << "\n"; 
  }
  return os;
}
