#if HAVE_MPI_H
  #include "../parallel/mpi_init.h"
  #include "../parallel/mpi_interface.h"
#endif
#include <gtest/gtest.h>

#ifdef TIMEING
  #include "../peps/extern.h"
#endif

#ifdef TIMEING
Stopwatch<> TIME_SVD;
Stopwatch<> TIME_PERMUTE;
Stopwatch<> TIME_CONTRACT;
Stopwatch<> TIME_ALLREDUCE;
#endif

int main(int argc, char **argv) {
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
#if HAVE_MPI_H
  mpi_init();
#endif
  
  result = RUN_ALL_TESTS();
#ifdef TIMEING
  std::cout << "TIMEING::PERMUTE " << TIME_PERMUTE.elapsed().count() *1.0E-9 << std::endl;
  std::cout << "TIMEING::CONTRACT " << TIME_CONTRACT.elapsed().count() *1.0E-9 << std::endl;
  std::cout << "TIMEING::SVD " << TIME_SVD.elapsed().count() *1.0E-9 << std::endl;
  std::cout << "TIMEING::ALLREDUCE " << TIME_ALLREDUCE.elapsed().count() *1.0E-9 << std::endl;
#endif
  return result;
}
