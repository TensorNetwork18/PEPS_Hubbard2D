#ifndef _TAG_GENERATOR_H_
#define _TAG_GENERATOR_H_

#include <sys/stat.h>
#include <stdio.h>
#include <stdexcept>
#include <iostream>
#include <string>
#include "../mpi/mpi_interface.h"

class TAG_Generator {
  private:
    std::string prefix_ = "./T" + std::to_string(mpi__->rank()) + ".";
    size_t count_ = 0;
  public:
    TAG_Generator()=default;
    ~TAG_Generator()=default;
    //~TAG_Generator() {
    //  reset();
    //}

    void reset() {
      for (size_t i = 0; i != count_; ++i) {
        struct stat buf;
        const std::string filename = prefix_ + std::to_string(i);
        if (stat(filename.c_str(), &buf) == 0) {
          if( remove( filename.c_str() ) != 0 )
           throw std::runtime_error("TAG_Generator: Error deleting " + filename);
        }
      }
      count_ = 0;
      return;
    }

    std::string get() {
      if (!(count_ < 4294967295)) //maxmum number of files
        throw std::runtime_error("TAG_Generator: Something is very wrong");
      return prefix_ + std::to_string(count_++);
    }
};

extern TAG_Generator* tag__;
#endif
