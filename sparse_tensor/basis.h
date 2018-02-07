#ifndef BASIS_H_
#define BASIS_H_

#include <iostream>
#include <vector>
#include <map>
#include <iterator>
#include <numeric>
#include "qnum.h"

class Basis {
  public:
    enum Type { OUT, IN };
  private:
    Type type_;
    std::map<Qnum, size_t> basis_;
  public:
    virtual ~Basis()=default;
    Basis(const Type type = IN);
    Basis(const Type, const size_t); //No Symmetry
    Basis(const Type, const std::vector<Qnum>&);

    Basis(const Basis&);
    Basis(Basis&&);
    Basis& operator = (const Basis&);
    Basis& operator = (Basis&&);
    bool operator == (const Basis&) const;

    Type type() const { return type_; }
    const std::map<Qnum, size_t>& basis() const { return basis_; }
    std::map<Qnum, size_t>& basis() { return basis_; }
    
    size_t size() const { return basis_.size(); }
    size_t qdim(const Qnum& o) const { auto it = find(o); return (it != basis_.end()) ? it->second : 0; }
    size_t total_dim() const { return std::accumulate(basis().begin(), basis().end(), 0, [](size_t value, const std::map<Qnum, size_t>::value_type& e){ return value + e.second; }); }

    std::map<Qnum, size_t>::const_iterator find(const Qnum& o) const { return basis_.find(o); }
    std::map<Qnum, size_t>::iterator find(const Qnum& o) { return basis_.find(o); }
    bool has_basis(const Qnum& o) const { return (basis_.count(o) == 1); }
    bool insert(const Qnum qn, const size_t qd) { return basis_.emplace(qn, qd).second; }
    bool remove(const Qnum qn) { return (basis_.erase(qn) == 1); }
    void assign(Type, const size_t);
    void assign(Type, const std::vector<Qnum>&);

    Basis reverse() const;
    Basis inverse() const;
    Basis fuse(const Basis&) const;   

    friend std::ostream& operator << (std::ostream& os, const Basis&);
    void broadcast(int source);
  private:
    void init(const std::vector<Qnum>&);
};
#endif
