#include <iostream>
#include <functional>
#include <utility>
#include <vector>
#include <set>
#include <numeric>
#include "basis.h"
#include "qtensor.h"
#ifdef HAVE_TBB
  #include "../thread/tbb_interface.h"
#endif

using namespace std;

vector<Qtensor> Qtensor::svd(const size_t split_rank, const Basis::Type new_basis_type, const char* shift_charge) const {
  //TODO remove blocks which are zero
  assert(split_rank > 0 && split_rank < rank());
  
  const size_t nrank = rank();
  // Construct output //
  Qtensor u, s, vt;
  if (shift_charge[0] == 'u' || shift_charge[0] == 'U') {
    u.charge_ = charge_;
  } else if (shift_charge[0] == 'v' || shift_charge[0] == 'V') {
    vt.charge_ = charge_;
  } else {
    throw logic_error("Qtensor::svd fails @ shift_charge");
  }

  Basis new_basis(new_basis_type); 
  //u: i_0, i_1, ..., i_split_rank-1, reversed_new_basis
  u.index_.resize(split_rank+1);
  copy_n(index().begin(), split_rank, u.index_.begin());
  u.index_.back() = new_basis.reverse();
  //s: new_basis, reversed_new_basis
  s.index_ = {new_basis, new_basis.reverse()};
  //v: new_basis, i_split_rank, ..., i_nrank-1
  vt.index_.resize(nrank - split_rank + 1);
  vt.index_.front() = new_basis;
  copy_n(index_.begin() + split_rank, nrank - split_rank, vt.index_.begin()+1);

  map<Qnum, vector<pair<vector<Qnum>, shared_ptr<Tensor>>>> svd_map;
  load();
  for (const auto& block : data_) {
    const vector<Qnum>& key = block.first;
    Qnum qsum = u.charge();
    for (size_t i = 0; i != split_rank; ++i)
      qsum = (index_[i].type() == Basis::IN) ? qsum + key[i] : qsum - key[i];
    if (new_basis_type == Basis::OUT)
      qsum = -qsum;

    auto it = svd_map.find(qsum);
    if (it != svd_map.end()) {
      it->second.push_back(block);
    } else {
      vector<pair<vector<Qnum>, shared_ptr<Tensor>>> tmp = {block};
      svd_map.emplace(move(qsum), move(tmp));
    }
  }
  dump();

#ifdef HAVE_TBB
  vector<function<void()>> tasks;
  tasks.reserve(svd_map.size());
  tbb::queuing_mutex scopedMutex;
#endif

  for (auto&& block_svd : svd_map) {
#ifdef HAVE_TBB
    tasks.emplace_back( [&split_rank, &u, &s, &vt, &scopedMutex, &block_svd]() {
#endif

    const Qnum& key  = block_svd.first;
    const vector<pair<vector<Qnum>, shared_ptr<Tensor>>>& value = block_svd.second;
    const size_t nblocks = value.size();

    // Map into a matrix //
    // Build quantum labels & tensor extents
    vector<vector<Qnum>> row_labels, col_labels;
    vector<vector<size_t>> row_extents, col_extents;
    vector<size_t> row_positions(nblocks), col_positions(nblocks);

    size_t iblock = 0;
    for (const auto& ref : value) {
      const auto& ref_T = *ref.second;
      { //row
        vector<Qnum> rlabel(ref.first.begin(), ref.first.begin() + split_rank);
        auto it = find(row_labels.begin(), row_labels.end(), rlabel);
        if (it != row_labels.end()) {
          row_positions[iblock] = distance(row_labels.begin(), it);
        } else {
          row_positions[iblock] = row_labels.size();
          row_labels.push_back(move(rlabel));
          vector<size_t> rextent(ref_T.extent().begin(), ref_T.extent().begin() + split_rank);
          row_extents.push_back(move(rextent));
        }
      }
      { //col
        vector<Qnum> clabel(ref.first.begin() + split_rank, ref.first.end());
        auto it = find(col_labels.begin(), col_labels.end(), clabel);
        if (it != col_labels.end()) {
          col_positions[iblock] = distance(col_labels.begin(), it);
        } else {
          col_positions[iblock] = col_labels.size();
          col_labels.push_back(move(clabel));
          vector<size_t> cextent(ref_T.extent().begin() + split_rank, ref_T.extent().end());
          col_extents.push_back(move(cextent));
        }
      }
      ++iblock;
    }

    // Compute row & col offsets
    vector<size_t> row_offsets(row_labels.size());
    for (size_t i = 0, sum = 0; i != row_labels.size(); ++i) {
      sum += accumulate( row_extents[i].begin(), row_extents[i].end(), 1, multiplies<size_t>());     
      row_offsets[i] = sum;
    }
    vector<size_t> col_offsets(col_labels.size());
    for (size_t i = 0, sum = 0; i != col_labels.size(); ++i) {
      sum += accumulate( col_extents[i].begin(), col_extents[i].end(), 1, multiplies<size_t>());     
      col_offsets[i] = sum;
    }

    vector<shared_ptr<Matrix>> ans;
    {
      // Construct the Matrix & SVD
      Matrix mat(row_offsets.back(), col_offsets.back());
      iblock = 0;
      for (const auto& ref : value) {
        const auto& ref_T = *ref.second;
        const size_t ir = row_positions[iblock];
        const size_t ic = col_positions[iblock];
        size_t row_start = (ir == 0) ? 0 : row_offsets[ir-1];
        size_t row_end = row_offsets[ir];
        size_t col_start = (ic == 0) ? 0 : col_offsets[ic-1];
        size_t col_end = col_offsets[ic];
        assert(ref_T.size() == (row_end - row_start) * (col_end - col_start));

        mat.copy_block(row_start, col_start, row_end - row_start, col_end - col_start, ref_T.pdata());
        ++iblock;
      }
      ans = mat.jsvd();
    }
    block_svd.second.clear();

    // Map matrices to tensors
    const size_t new_dim = ans[1]->nrows();
    // prepare u
    vector<pair<vector<Qnum>, shared_ptr<Tensor>>> targets_u(row_labels.size());
    {
      shared_ptr<const Matrix> mat_u = move(ans[0]);
      size_t icount = 0;
      for (auto&& row_label : row_labels) {
        size_t row_start = (icount == 0) ? 0 : row_offsets[icount-1];
        size_t row_end = row_offsets[icount];
        const shared_ptr<Matrix>& sub_matrix = mat_u->cut_row(row_start, row_end);

        vector<size_t> u_extent(move(row_extents[icount]));
        u_extent.push_back(new_dim);
        auto u_tensor = make_shared<Tensor>(u_extent);
        assert(u_tensor->size() == sub_matrix->size());
        u_tensor->set_elements(sub_matrix->size(), sub_matrix->pdata());

        vector<Qnum> u_label(move(row_label));
        u_label.push_back(key);

        targets_u[icount] = make_pair(move(u_label), move(u_tensor));
        ++icount;
      }
    } 

    // prepare s
    vector<pair<vector<Qnum>, shared_ptr<Tensor>>> targets_s(1);
    {
      shared_ptr<const Matrix> mat_s = move(ans[1]);
      auto s_tensor = make_shared<Tensor>(new_dim, new_dim);
      for (size_t i = 0; i != new_dim; ++i)
        s_tensor->element(i,i) = mat_s->element(i);
      vector<Qnum> s_label(2, key);
      targets_s[0] = make_pair(move(s_label), move(s_tensor));
    }

    // prepare vt
    vector<pair<vector<Qnum>, shared_ptr<Tensor>>> targets_vt(col_labels.size());
    {
      shared_ptr<const Matrix> mat_vt = move(ans[2]);
      size_t icount = 0;
      for (auto&& col_label : col_labels) {
        size_t col_start = (icount == 0) ? 0 : col_offsets[icount-1];
        size_t col_end = col_offsets[icount];
        const shared_ptr<Matrix>& sub_matrix = mat_vt->cut_col(col_start, col_end);
        
        vector<size_t> vt_extent(move(col_extents[icount]));
        vt_extent.insert(vt_extent.begin(), new_dim);
        auto vt_tensor = make_shared<Tensor>(vt_extent);
			  assert(vt_tensor->size() == sub_matrix->size());
        vt_tensor->set_elements(sub_matrix->size(), sub_matrix->pdata());

        vector<Qnum> vt_label(move(col_label));
        vt_label.insert(vt_label.begin(), key);

        targets_vt[icount] = make_pair(move(vt_label), move(vt_tensor));
        ++icount;
      }
    }

    {
#ifdef HAVE_TBB
      tbb::queuing_mutex::scoped_lock lock(scopedMutex);
#endif
      u.index_.back().insert(key, new_dim);
      for(auto&& i : targets_u)
        u.data_.emplace(move(i.first), move(i.second));
      
      s.index_[0].insert(key, new_dim);
      s.index_[1].insert(key, new_dim);
      for(auto&& i : targets_s)
        s.data_.emplace(move(i.first), move(i.second));

      vt.index_.front().insert(key, new_dim);
      for(auto&& i : targets_vt)
        vt.data_.emplace(move(i.first), move(i.second));
    }
#ifdef HAVE_TBB
    });
#endif
  }

#ifdef HAVE_TBB
  parallel_for<>(tasks);
#endif

  return {move(u), move(s), move(vt)};
}
