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


void Qtensor::contract(const double alpha, const Qtensor& A, const vector<size_t>& aA, const Qtensor& B, const vector<size_t>& aB, const bool sgn) {
  const size_t rank_common = aA.size();
  const size_t rankA = A.rank();
  const size_t rankB = B.rank();
  const size_t rankC = rankA + rankB - 2 * rank_common;

#ifndef NDEBUG
  assert(rank_common == aB.size());
  for (size_t i = 0; i != rank_common; ++i) {
    assert(A.index(aA[i]).type() != B.index(aB[i]).type());
    assert(A.index(aA[i]).basis() == B.index(aB[i]).basis());
  }
  if (rankC == 0)
    assert(A.charge_ + B.charge_ == Qnum::zero());
#endif

  Qtensor out;
  out.charge_ = A.charge_ + B.charge_;
 
  // Create output Qtensor //
  vector<size_t> removeA(aA);
  vector<size_t> removeB(aB);
  sort(removeA.begin(), removeA.end());
  sort(removeB.begin(), removeB.end());
  if (rankC != 0) {
    out.index_ = A.index();
    vector<Basis> lB(B.index());
    for (size_t i = rank_common; i-- > 0;) {
      out.index_.erase(out.index_.begin()+removeA[i]);
      lB.erase(lB.begin()+removeB[i]);
    }
    out.index_.insert(out.index_.end(), lB.begin(), lB.end());
  } else { //scalar
    Basis li(Basis::IN, 1u);
    Basis lo(Basis::OUT, 1u);
    out.index_ = {li, lo};
  }

#ifdef HAVE_TBB
  const size_t nblocksA = A.block_size();
  const size_t nblocksB = B.block_size();
  std::vector<std::pair<vector<Qnum>, shared_ptr<Tensor>>> targetsA(nblocksA);
  std::vector<std::pair<vector<Qnum>, shared_ptr<Tensor>>> targetsB(nblocksB);

  if (rank_common != 0) {
    vector<size_t> pmtA(A.rank());
    iota(pmtA.begin(), pmtA.end(), 0); 
    vector<size_t> pmtB(B.rank());
    iota(pmtB.begin(), pmtB.end(), 0); 
    for (size_t i = rank_common; i-- > 0;) {
      pmtA.erase(pmtA.begin()+removeA[i]);
      pmtB.erase(pmtB.begin()+removeB[i]);
    }   
    pmtA.insert(pmtA.end(), aA.begin(), aA.end());
    pmtB.insert(pmtB.begin(), aB.begin(), aB.end());

    //Permute//
    vector<function<void()>> tasks;
    tasks.reserve(nblocksA + nblocksB);

    // Permute A//
    const vector<Basis>& indexA = A.index_;
    size_t icount = 0;
    const bool& diskA = A.disk_;
    for (auto&& block : A.data_) {
      tasks.emplace_back( [&pmtA, &sgn, &rankA, &rank_common, &indexA, &block, &targetsA, icount, &diskA]() {
        const vector<Qnum>& key = block.first;
        const Tensor& ref = *block.second;

        vector<Qnum> new_key(rankA);
        vector<size_t> new_extent(rankA);
        vector<size_t> swap_gate;
        for (size_t k = 0; k != rankA - rank_common; ++k) {
          const size_t p = pmtA[k];
          const Qnum q = key[p];
          new_key[k]  = q;
          new_extent[k] = ref.extent(p);
          if (sgn && q.parity() == Qnum::ODD)
            swap_gate.push_back(p);
        }

        bool change_sgn = false;
        for (size_t k = rankA; k-- > rankA - rank_common;) { //In reversed order
          const size_t p = pmtA[k];
          const Qnum q = key[p];
          new_key[k]  = q;
          new_extent[k] = ref.extent(p);
          if (sgn && q.parity() == Qnum::ODD) {
            swap_gate.push_back(p);
            if (indexA[p].type() == Basis::OUT)
              change_sgn ^= true;
          }
        }

        auto ptr = make_shared<Tensor>(new_extent);
        sort_arr<>(ref.pdata(), ptr->pdata(), pmtA, ref.extent());

        if (diskA) {
          assert(block.second.use_count() == 1);
          block.second = nullptr;          
        }

        if (sgn && (change_sgn ^ sgn_permute(swap_gate)))
          ptr->scale(-1.0);

        auto&& target = targetsA[icount];
        target.first = move(new_key);
        target.second = move(ptr);
      });
      ++icount;
    }
    assert(icount == nblocksA);

    A.load();
    parallel_for<>(tasks);
    A.dump();
    tasks.clear();    

    // Permute B //
    icount = 0;
    const bool& diskB = B.disk_;
    for (auto&& block : B.data_) {
      tasks.emplace_back( [&pmtB, &sgn, &rankB, &block, &targetsB, icount, &diskB]() {
        const vector<Qnum>& key = block.first;
        const Tensor& ref = *block.second;
    
        vector<Qnum> new_key(rankB);
        vector<size_t> new_extent(rankB);
        vector<size_t> swap_gate;
        for (size_t k = 0; k != rankB; ++k) {
          const size_t p = pmtB[k];
          const Qnum q = key[p];
          new_key[k] = q;
          new_extent[k] = ref.extent(p);
          if (sgn && q.parity() == Qnum::ODD)
            swap_gate.push_back(p);
        }   

        auto ptr = make_shared<Tensor>(new_extent);
        sort_arr<>(ref.pdata(), ptr->pdata(), pmtB, ref.extent());

        if (diskB) {
          assert(block.second.use_count() == 1);
          block.second = nullptr;     
        }

        if (sgn && sgn_permute(swap_gate))
          ptr->scale(-1.0);

        auto&& target = targetsB[icount];
        target.first = move(new_key);
        target.second = move(ptr);
      }); 
      ++icount;
    }   
    assert(icount == nblocksB);

    B.load();
    parallel_for<>(tasks);
    B.dump();

  } else {

    size_t icount = 0;
    A.load();
    for (const auto& block : A.data_) {
      targetsA[icount].first = block.first;
      targetsA[icount].second = block.second->copy();
      ++icount;
    }
    assert(icount == A.block_size());
    A.dump();

    icount = 0;
    B.load();
    for (const auto& block : B.data_) {
      targetsB[icount].first = block.first;
      targetsB[icount].second = block.second->copy();
      ++icount;
    }
    B.dump();
    assert(icount == B.block_size());
  }

  // Build contract_map //
  map<vector<Qnum>, vector<pair<shared_ptr<Tensor>, shared_ptr<Tensor>>>> contract_map;
  if (rank_common != 0) {
    // Sort A //
    map<vector<Qnum>, vector<pair<vector<Qnum>, shared_ptr<Tensor>>>> sortA;
    for (auto&& targetA : targetsA) {
      vector<Qnum> new_key(targetA.first.end() - rank_common, targetA.first.end());
      auto it = sortA.find(new_key);
      if (it != sortA.end())
        it->second.push_back(targetA);
      else
        sortA.emplace(move(new_key), vector<pair<vector<Qnum>, shared_ptr<Tensor>>> {targetA});
    }

    // Sort B //
    map<vector<Qnum>, vector<pair<vector<Qnum>, shared_ptr<Tensor>>>> sortB;
    for (auto&& targetB : targetsB) {
      vector<Qnum> new_key(targetB.first.begin(), targetB.first.begin() + rank_common);
      auto it = sortB.find(new_key);
      if (it != sortB.end())
        it->second.push_back(targetB);
      else
        sortB.emplace(move(new_key), vector<pair<vector<Qnum>, shared_ptr<Tensor>>> {targetB});
    }

    if (rankC != 0) {
      const size_t lengthA = rankA - rank_common;
      const size_t lengthB = rankB - rank_common;
      for (const auto& i_A : sortA) {
        auto it_B = sortB.find(i_A.first);
        if (it_B != sortB.end()) {
          for (const auto& targetA : i_A.second) {
            for (const auto& targetB : it_B->second) {
              vector<Qnum> keyC(rankC);
              copy_n(targetA.first.begin(), lengthA, keyC.begin());
              copy_n(targetB.first.begin() + rank_common, lengthB, keyC.begin() + lengthA);

              auto it = contract_map.find(keyC);
              if (it != contract_map.end())
                it->second.emplace_back(targetA.second, targetB.second);
              else
                contract_map.emplace( move(keyC), vector<pair<shared_ptr<Tensor>, shared_ptr<Tensor>>> { make_pair(targetA.second, targetB.second) } );
            }
          }
          sortB.erase(it_B);
        }
      }

    } else {

      vector<Qnum> keyC(2, Qnum::zero());
      vector<pair<shared_ptr<Tensor>, shared_ptr<Tensor>>> tmp;
      for (const auto& i_A : sortA) {
        auto it_B = sortB.find(i_A.first);
        if (it_B != sortB.end()) {
          assert(i_A.second.size() == 1);
          assert(it_B->second.size() == 1);
          assert(i_A.second[0].first == it_B->second[0].first);
          tmp.push_back( make_pair(i_A.second[0].second, it_B->second[0].second) );
          sortB.erase(it_B);
        }
      }
      contract_map.emplace(move(keyC), move(tmp));
    }  
 
  } else {
   
    for (auto&& targetA : targetsA) {
      for (auto&& targetB : targetsB) {
        vector<Qnum> keyC(rankC);
        copy_n(targetA.first.begin(), rankA, keyC.begin());
        copy_n(targetB.first.begin(), rankB, keyC.begin() + rankA);

        assert(contract_map.count(keyC) == 0);
        contract_map.emplace(move(keyC), vector<pair<shared_ptr<Tensor>, shared_ptr<Tensor>>> { make_pair(targetA.second, targetB.second) } );
      }
    }
  }

  // Contract //
  if (rankC != 0) {
   const size_t nblocks = contract_map.size();
    vector<function<void()>> tasks;
    tasks.reserve(nblocks);
    vector<pair<vector<Qnum>, shared_ptr<Tensor>>> targets(nblocks);

    const size_t lengthA = rankA - rank_common;
    const size_t lengthB = rankB - rank_common;
    size_t icount = 0;
    //tbb::spin_mutex scopedMutex;

    for (auto&& ref : contract_map) {
      tasks.emplace_back( [&alpha, &rank_common, &rankC, &lengthA, &lengthB, &ref, &targets, icount]() {
        const size_t njobs = ref.second.size();
        
        //creat new target tensor
        const Tensor a0 = *ref.second[0].first;
        const Tensor b0 = *ref.second[0].second;
        vector<size_t> new_extent(rankC);
        copy_n(a0.extent().begin(), lengthA, new_extent.begin());
        copy_n(b0.extent().begin() + rank_common, lengthB, new_extent.begin()+lengthA);
        auto new_value = make_shared<Tensor>(new_extent);

        //dgemm
        size_t m = 1; for_each( a0.extent().begin(), a0.extent().end() - rank_common, [&](size_t i) {m *= i;} );
        size_t n = 1; for_each( b0.extent().begin() + rank_common, b0.extent().end(), [&](size_t i) {n *= i;} );

        for (size_t j = 0; j != njobs; ++j) {
          const Tensor& aj = *ref.second[j].first;
          const Tensor& bj = *ref.second[j].second;
          size_t k = 1; for_each( bj.extent().begin(), bj.extent().begin() + rank_common, [&](size_t i) {k *= i;} );
          const double beta = 1.0;
          dgemm("N", "N", m, n, k, alpha, aj, bj, beta, *new_value);
        }

        //{ I don't think it's thread-safe
        //  tbb::spin_mutex::scoped_lock lock(scopedMutex);
        //  ref.second.clear();
        //}

        auto&& target = targets[icount];
        target.first = ref.first;
        target.second = move(new_value);
      });
      ++icount;
    }
    assert(icount == nblocks);

    parallel_for<>(tasks);

    for (auto&& target : targets)
      out.data_.emplace(move(target.first), move(target.second));

  } else {

    //only one possible block
    assert(contract_map.size() == 1);
    const auto& key = contract_map.begin()->first;
    const auto& value = contract_map.begin()->second;
    const size_t njobs = value.size();

    vector<size_t> new_extent = {1,1};
    auto new_value = make_shared<Tensor>(new_extent);

    for (size_t j = 0; j != njobs; ++j) {
      const Tensor& aj = *value[j].first;
      const Tensor& bj = *value[j].second;
      new_value->element(0,0) += ddot(aj, bj);
    }
 
    out.data_.emplace(key, move(new_value));
  }

#else

  // Permute //
  Qtensor refA, refB;
  if (rank_common != 0) {
    vector<size_t> pmtA(A.rank());
    iota(pmtA.begin(), pmtA.end(), 0);
    vector<size_t> pmtB(B.rank());
    iota(pmtB.begin(), pmtB.end(), 0);
    for (size_t i = rank_common; i-- > 0;) {
      pmtA.erase(pmtA.begin()+removeA[i]);
      pmtB.erase(pmtB.begin()+removeB[i]);
    }
    pmtA.insert(pmtA.end(), aA.begin(), aA.end());
    pmtB.insert(pmtB.begin(), aB.begin(), aB.end());

    refA = A.permute(pmtA, sgn);
    refB = B.permute(pmtB, sgn);
  } else {
    refA = A;
    refB = B;
  }

  const size_t nA = refA.rank();
  for (const auto& ia : refA.data_) {
    //Change of sgn//
    bool change_sgn = false;
    if (sgn) {
      vector<size_t> fswaps;
      for (size_t k = 0; k != rank_common; ++k) {
        if (ia.first[nA-k-1].parity() == Qnum::ODD) { //Trace
          fswaps.push_back(nA-k-1);
          if (refA.index(nA-k-1).type() == Basis::OUT) //Contraction order is reversed
            change_sgn ^= true;
        }
      }
      change_sgn ^= sgn_permute(fswaps);
    }
    const Tensor& tensorA = *ia.second;
    
    for (const auto& ib : refB.data_) {
      const bool match = equal(ib.first.begin(), ib.first.begin() + rank_common, ia.first.end() - rank_common);
      if (match) {
        vector<Qnum> new_key;
        vector<size_t> new_extent;
     
        const Tensor& tensorB = *ib.second;
        if (rankC != 0) {
          new_key.reserve(rankC);
          new_extent.reserve(rankC);
          new_key.insert(new_key.end(), ia.first.begin(), ia.first.end() - rank_common);
          new_key.insert(new_key.end(), ib.first.begin() + rank_common, ib.first.end());
          new_extent.insert(new_extent.end(), tensorA.extent().begin(), tensorA.extent().end() - rank_common);
          new_extent.insert(new_extent.end(), tensorB.extent().begin() + rank_common, tensorB.extent().end());
        } else {
          Qnum q0 = Qnum::zero();
          assert(out.charge_ == q0);
          new_key = {q0, q0};
          new_extent = {1,1};
        }
        
        assert(new_key.size() != 0);
        assert(new_extent.size() != 0);
        auto new_value = make_shared<Tensor>(new_extent);
        assert(new_value != nullptr);
        assert(new_value->size() != 0);
        size_t m = 1; for_each( tensorA.extent().begin(), tensorA.extent().end() - rank_common, [&](size_t i) {m *= i;} );
        size_t n = 1; for_each( tensorB.extent().begin() + rank_common, tensorB.extent().end(), [&](size_t i) {n *= i;} );
        size_t k = 1; for_each( tensorB.extent().begin(), tensorB.extent().begin() + rank_common, [&](size_t i) {k *= i;} );
        const double fa = (change_sgn) ? -alpha : alpha;
        const double beta = 0.0;
        dgemm("N", "N", m, n, k, fa, tensorA, tensorB, beta, *new_value);

        auto it = out.data_.find(new_key);
        if (it != out.data_.end()) {
          *it->second += *new_value;
        } else {
          out.data_.emplace(move(new_key), move(new_value));
        }
      }
    }
  }
#endif
  assert(out.block_size() != 0);

  *this = move(out);
  return;
}
