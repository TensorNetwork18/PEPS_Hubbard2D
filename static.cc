#include <cfenv>
#include "mpi/mpi_interface.h"
#include "shared_tensor/tag_generator.h"

using namespace std;

// GLOBAL //
TAG_Generator* tag__;
MPI_Interface* mpi__;
std::unique_ptr<TAG_Generator> tag;
std::unique_ptr<MPI_Interface> mpi;

void static_variables() {
  mpi = unique_ptr<MPI_Interface>(new MPI_Interface());
  mpi__ = mpi.get();

  tag = unique_ptr<TAG_Generator>(new TAG_Generator());
  tag__ = tag.get();

  // rounding mode in std::rint, std::lrint, and std::llrint
  fesetround(FE_TONEAREST);
}
