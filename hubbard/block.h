#ifndef _BLOCK_H_
#define _BLOCK_H_

#include <vector>
#include "../mpi/mpi_interface.h"
#include "../peps/mps2.h"
#include "etensor_pool.h"


namespace Hubbard2 {
struct CanonicalMPS {
    enum Type {
      Ham, AlphaC, AlphaA, BetaC, BetaA, Identity
    };
  
    Type type;
    double coeff;
    std::shared_ptr<MPS2> mps_ptr;
    int position; //Position of the creator or annihilator
  
    CanonicalMPS()=default;
    CanonicalMPS(const Type o1, const double o2, const std::shared_ptr<MPS2>& o3, const int o4) : type(o1), coeff(o2), mps_ptr(o3), position(o4) {}

    void broadcast(int source) {
      if (mpi__->size() != 1) {
        const int iproc = mpi__->world_rank();
        
        { // Type //
          int t = static_cast<int>(type);
          mpi__->broadcast(&t, 1u, source);
          if (iproc != source)
            type = static_cast<Type>(t);
        }

        if (type != Ham) { //only shared if the CanonicalMPS type is not Ham
          if (!mps_ptr)
            mps_ptr = make_shared<MPS2>();
          // Coeff
          mpi__->broadcast(&coeff, 1u, source);
          // MPS
          mps_ptr->broadcast(source);
          // Position
          mpi__->broadcast(&position, 1u, source);
        }
      }
      return;
    }

    bool exist() const {
      return (mps_ptr != nullptr);
    }
};


class Block {
  private:
    MPS2::Type type_;
    std::vector<CanonicalMPS> blocks_;
    
  public:
    virtual ~Block()=default;
    Block()=default;
    Block(const double u, const std::shared_ptr<Block>& o1, const EtensorPool& o2, const size_t o3) : Block(u, o1, o2, o3, {}) { }
    Block(const double, const std::shared_ptr<Block>&, const EtensorPool&, const size_t, const std::map<Qnum, size_t>&);
    Block(MPS2::Type type, const EtensorPool&, const size_t);

    const MPS2::Type& type() const { return type_; }
    const std::vector<CanonicalMPS>& blocks() const { return blocks_; }
    const CanonicalMPS& identity() const { 
      auto it = std::find_if(blocks_.begin(), blocks_.end(), [](const CanonicalMPS& i){ return i.type == CanonicalMPS::Identity; });
      if (it != blocks_.end())
        return *it;
      else
        throw std::logic_error("Fail @ Hubbard2::Block::identity");
    }
};
}
#endif
