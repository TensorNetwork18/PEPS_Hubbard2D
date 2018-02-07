#include <iostream>
#include <functional>
#include <utility>
#include <vector>
#include <set>
#include <numeric>
#include "basis.h"
#include "qtensor.h"
#include "../math/sort_arr.h"
#ifdef HAVE_TBB
  #include "../thread/tbb_interface.h"
#endif

using namespace std; 


Qtensor Qtensor::permute(const vector<size_t>& pmt, const bool sgn) const {
  const size_t nrank = rank();
  assert(pmt.size() == nrank);
  
  // Create output Qtensor //
  Qtensor out;
  out.charge_ = charge();

  // Permute indices //
  out.index_.resize(nrank);
  for (size_t i = 0; i != nrank; ++i)
    out.index_[i] = index(pmt[i]);

#ifdef HAVE_TBB
  const size_t nblocks = block_size();
  vector<pair<vector<Qnum>, shared_ptr<Tensor>>> targets(nblocks);

  vector<function<void()>> tasks;
  size_t icount = 0;
  const bool& disk = disk_;
  for (auto&& block : data_) {
    tasks.emplace_back( [&pmt, &sgn, &nrank, &block, &targets, icount, &disk]() {
      const vector<Qnum>& key = block.first;
      const Tensor& ref = *block.second;

      vector<Qnum> new_key(nrank);
      vector<size_t> new_extent(nrank);
      vector<size_t> swap_gate;
      for (size_t k = 0; k != nrank; ++k) {
        const size_t p = pmt[k];
        const Qnum q = key[p];
        new_key[k] = q;
        new_extent[k] = ref.extent(p);
        if (sgn && q.parity() == Qnum::ODD)
          swap_gate.push_back(p);
      }

      auto ptr = make_shared<Tensor>(new_extent);
      sort_arr<>(ref.pdata(), ptr->pdata(), pmt, ref.extent());

      if (disk) {//manually delete
        assert(block.second.use_count() == 1);
        block.second = nullptr;
      }

      if (sgn && sgn_permute(swap_gate))
        ptr->scale(-1.0);

      auto&& target = targets[icount];
      target.first = move(new_key);
      target.second = move(ptr);
    });
    ++icount;
  }
  assert (icount == nblocks);

  load();
  parallel_for<>(tasks);
  dump();

  for (auto&& target : targets)
    out.data_.emplace(move(target.first), move(target.second));

#else

  load();
  for (auto&& block : data_) {
    const Tensor& ref = *block.second;

    //Permute coordinates
    vector<Qnum> new_key(nrank);
    vector<size_t> new_extent(nrank);
    vector<size_t> swap_gate;
    for (size_t k = 0; k != nrank; ++k) { 
      new_key[k] = block.first[pmt[k]];
      new_extent[k] = ref.extent(pmt[k]);
      if (sgn && new_key[k].parity() == Qnum::ODD)
        swap_gate.push_back(pmt[k]);
    }

    //Permute block elements
    auto new_ptensor = make_shared<Tensor>(new_extent);
    sort_arr<>(ref.pdata(), new_ptensor->pdata(), pmt, ref.extent());
    
    if (disk_) {//manulally delete
      assert(block.second.use_count() == 1);
      block.second = nullptr;
    }

    if (sgn && sgn_permute(swap_gate)) //Parity changed
      new_ptensor->scale(-1.0);
    out.data_.emplace(move(new_key), move(new_ptensor));
  }
  dump();
#endif
  return out;
}
