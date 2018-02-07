#ifndef _BINARY_IO 
#define _BINARY_IO

#include <iostream>

template<typename T>
std::ostream& binary_write(std::ostream& stream, const T* value, const size_t s = 1u) {
    return stream.write(reinterpret_cast<const char*>(value), sizeof(T)*s);
}

template<typename T>
std::istream & binary_read(std::istream& stream, T* value, const size_t s = 1u){
    return stream.read(reinterpret_cast<char*>(value), sizeof(T)*s);
}

#endif
