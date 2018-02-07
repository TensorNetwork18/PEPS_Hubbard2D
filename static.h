#ifndef _STATIC_H
#define _STATIC_H

#include <iostream>
#include <cstdlib>

void static_variables();

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
