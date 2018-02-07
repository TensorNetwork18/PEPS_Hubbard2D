#ifndef _MPI_INIT_H_
#define _MPI_INIT_H_

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <cstdlib>

void mpi_init();

template<typename T>
std::string getenv_multiple(const T& head) {
  char const* val = std::getenv(head);
  return val ? std::string(val) : "";
}


template<typename T, typename ...args>
std::string getenv_multiple(const T& head, const args&... tail) {
  char const* val = std::getenv(head);
  return val ? std::string(val) : getenv_multiple(tail...);
}

#endif
