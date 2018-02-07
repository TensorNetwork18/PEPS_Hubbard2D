#include <cfenv>
#include <thread>
#include "mpi_init.h"
#include "mpi_interface.h"

using namespace std;


MPI_Interface* mpi__;
std::unique_ptr<MPI_Interface> mpi;

void mpi_init() {
  // setup MPI interface. It does nothing for serial runs
  mpi = unique_ptr<MPI_Interface>(new MPI_Interface());
  mpi__ = mpi.get();
  // rounding mode in std::rint, std::lrint, and std::llrint
  fesetround(FE_TONEAREST);
}
