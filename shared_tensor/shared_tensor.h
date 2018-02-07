#ifndef _SHARED_TENSOR_H_
#define _SHARED_TENSOR_H_

#include <memory>
#include <string>
#include <fstream>
#include "../binary_io.h"
#include "../math/tensor.h"
#include "tag_generator.h"

class Shared_Tensor {
  private:
    mutable std::shared_ptr<Tensor> ptensor_;
    mutable bool modified_ = false;
    bool disk_ = false;
    std::string file_; 


    void write() const {
      ofstream output;
      output.open(file_, ios::binary);
      if (!output.good())
        throw runtime_error("Shared_Tensor::write fails " + file_); 

      //Rank
      const size_t rank = ptensor_->rank();
      binary_write(output, &rank, 1);

      //Extent
      const std::vector<size_t>& extent = ptensor_->extent();
      binary_write(output, &extent[0], rank);

      //Element
      const double* pdata = ptensor_->pdata();
      binary_write(output, pdata, ptensor_->size());

      output.close();
      return;
    };

    void read() const {
      ifstream input;
      input.open(file_, ios::binary);
      if (!input.good())
        throw runtime_error("Shared_Tensor::read fails " + file_); 
      
      size_t rank;
      binary_read(input, &rank, 1);
     
      std::vector<size_t> extent(rank);
      binary_read(input, &extent[0], rank);
      
      ptensor_ = make_shared<Tensor>(extent);
      double* pdata = ptensor_->pdata();
      binary_read(input, pdata, ptensor_->size()); 

      input.close();
      return;
    }

  public:
    virtual ~Shared_Tensor()=default;
    Shared_Tensor()=default;
    Shared_Tensor(const std::shared_ptr<Tensor>& ptensor_in, const bool disk_in) : ptensor_(ptensor_in), modified_(true), disk_(disk_in) {
      if (disk_) {
        assert(ptensor_ != nullptr);
        file_ = tag__->get();
      }
      dump();
    }

    Shared_Tensor(const Shared_Tensor& obj) : ptensor_(obj.ptensor_), disk_(obj.disk_), file_(obj.file_) {
      dump();
    }

    Shared_Tensor(Shared_Tensor&& obj) : ptensor_(std::move(obj.ptensor_)), disk_(obj.disk_), file_(std::move(obj.file_)) { 
       dump();
       obj.disk_ = false;
    }
    
    Shared_Tensor& operator = (const Shared_Tensor& obj) {
      if (this != &obj) {
        ptensor_ = obj.ptensor_;
        disk_ = obj.disk_; 
        file_ = obj.file_;
        dump();
      }
      return *this;
    }

    Shared_Tensor& operator = (Shared_Tensor&& obj) {
      if (this != &obj) {
        ptensor_ = std::move(obj.ptensor_);
        disk_ = obj.disk_; 
        file_ = std::move(obj.file_);
        dump();
        obj.disk_ = false;
      } 
      return *this;
    }
 
    Shared_Tensor copy(const bool disk_in) const {
      if (disk_ && ptensor_ == nullptr) {
        assert(file_.length() != 0);
        read();
      }
      if (ptensor_ != nullptr) {
        auto out = ptensor_->copy();
        if (disk_)
          ptensor_ = nullptr;
        return Shared_Tensor(out, disk_in);
      }
      return Shared_Tensor(); 
    }

    int use_count() const {
      return ptensor_.use_count();
    }

    void delete_file() const {
      if (file_.length() != 0) {
        struct stat buf;
        if (stat(file_.c_str(), &buf) == 0) {
          if( remove( file_.c_str() ) != 0 )
           throw std::runtime_error("Shared_Tensor: Error deleting " + file_);
        }
      }
      return;
    }

    const Tensor& access() const {
     assert(ptensor_ != nullptr);
     return *ptensor_;
    }

    Tensor& access() {
     assert(ptensor_ != nullptr);
     modified_ = true;
     return *ptensor_;
    }

    void load() const {
      if (disk_ && ptensor_ == nullptr) {
        assert(file_.length()!=0);
        read();
      }
      return;
    }

    void dump() const {
      if (disk_) {
        assert(file_.length()!=0);
        if (modified_) {
          write();
          modified_ = false;     
        }
        ptensor_ = nullptr;
      }
      return;    
    }

    void flush() {
      if (disk_) {
        dump();
      } else {
        if (ptensor_ != nullptr) {
          file_ = tag__->get();
          write();
          disk_ = true;
          modified_ = false;
          ptensor_ = nullptr;
        }
      }
      return;
    }
};
#endif
