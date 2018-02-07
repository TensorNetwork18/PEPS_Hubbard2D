#include <cfenv>
#include <thread>
#include <string>
#include <iostream>
#include <gtest/gtest.h>
#include "mpi_init.h"
#include "mpi_interface.h"

using namespace std;


int main (int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  mpi_init();
  try {
    cout << "single ? " << endl;
  } catch (...) {
  }
  return 0;
}
