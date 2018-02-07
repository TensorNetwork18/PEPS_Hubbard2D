#include <fstream>
#include "../timer.h"
#include "minimizeE.h"

Hubbard2::MinimizeE::MinimizeE(const double u, const shared_ptr<PEPS>& psi, const std::map<Qnum, size_t>& cmap) : u_term_(u), psi_(psi), etensor_pool_(psi), compress_map_(cmap) {
}

void Hubbard2::MinimizeE::compute() {
  const size_t nr = psi_->nrows();
  const size_t nc = psi_->ncols();  

  if (!sweep_ene_.empty())
    sweep_ene_.clear();

  {
    Timer time_init;
    //Construct Boundary Block
    pBlocksU_.emplace_back(make_shared<Block>(MPS2::U, etensor_pool_, 0u));
    pBlocksD_.emplace_back(make_shared<Block>(MPS2::D, etensor_pool_, nr-1u));

    //Fill up 
    for (size_t i = nr; --i > 0;) {
      cout << "  Compute Block " << i << endl;
      pBlocksD_.emplace_back(make_shared<Block>(u_term_, pBlocksD_.back(), etensor_pool_, i, compress_map_));
    }
    cout << "  Time for Initialization:   " << time_init.tick() << endl << endl;
  }

  { //up to down
    for (size_t i = 0; i != nr; ++i) {
      pEffTN_ = make_shared<EffectiveTN>(u_term_, pBlocksU_.back(), pBlocksD_.back(), etensor_pool_, i);
      //left to right
      for (size_t j = 0; j != nc; ++j) {
        update_one_site(i, j, false);
      }
      if ( i != nr-1 ) {
        pBlocksU_.emplace_back(make_shared<Block>(u_term_, pBlocksU_.back(), etensor_pool_, i, compress_map_));
        pBlocksD_.pop_back();
      }
    }
    psi_->save();
  }

  if (mpi__->rank() == 0) {
    const double average_ene = accumulate(sweep_ene_.begin(), sweep_ene_.end(), 0.0) / sweep_ene_.size();
    cout << "  Mean_Energy = " << setw(15) << fixed << setprecision(8) << average_ene << endl;
  
    vector<double> diff(sweep_ene_.size());
    transform(sweep_ene_.begin(), sweep_ene_.end(), diff.begin(), [average_ene](double e) { return e - average_ene; });
    double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double dev = sqrt(sq_sum / sweep_ene_.size());
    cout << "  Deviation = " << setw(15) << fixed << setprecision(8) << dev << endl;

    const double min_ene = *min_element(sweep_ene_.begin(),sweep_ene_.end());
    cout << "  Minimum_Energy = " << setw(15) << fixed << setprecision(8) << min_ene << endl;
  }  

  return;
}

pair<bool, double> Hubbard2::MinimizeE::iteration(const double perturb, const double thresh, const size_t max_iter) {
  noise_ = perturb;
  last_average_ene_ = 0.0;
  conv_ = false;
  if (!sweep_ene_.empty())
    sweep_ene_.clear();

  const size_t nr = psi_->nrows();
  const double noise_min = thresh; // if (noise_ < noise_min) noise_ = 0.0; 

  {
    Timer time_init;
    //Construct Boundary Block
    pBlocksU_.emplace_back(make_shared<Block>(MPS2::U, etensor_pool_, 0u));
    pBlocksD_.emplace_back(make_shared<Block>(MPS2::D, etensor_pool_, nr-1u));
 

    //Fill up 
    for (size_t i = nr; --i > 0;) {
      cout << "  Compute Block " << i << endl;
      pBlocksD_.emplace_back(make_shared<Block>(u_term_, pBlocksD_.back(), etensor_pool_, i, compress_map_));
    }
    cout << "  Time for Initialization:   " << time_init.tick() << endl << endl;
  }

  for (size_t iter = 0; iter != max_iter; ++iter) {
    Timer time_iter;
    up_to_down_iter(0);
    down_to_up_iter(nr-1);
    
    //Dump PEPS
    psi_->save();

    //Check Convergence
    if (mpi__->rank() == 0) {
      conv_ = (noise_ < noise_min);
      const double energy_diff = check_convergence();
      conv_ &= energy_diff < thresh;
      const bool drop_noise =  energy_diff < max(1.0E-4*noise_, thresh); //probably not a good idea
      if (noise_ != 0.0) {
        if (drop_noise) {
          noise_ *= 0.1;
          if (noise_ < noise_min) {
            cout << "perturbative noise turned off" << endl;
            noise_ = 0.0;
          } else {
            cout << "perturbative noise lowered to: " << noise_ << endl;
          }   
        }   
      }
    }

    if (mpi__->size() != 1){
      const int iproc = mpi__->world_rank();
      int t = static_cast<int>(conv_);
      mpi__->broadcast(&t, 1u, 0);
      if (iproc != 0)
        conv_ = static_cast<bool>(t);
    }
    
    cout << "  Time for Iteration " << iter << ":    " << time_iter.tick() << endl << endl; 
    if (conv_) break;
  }
  return make_pair(conv_, noise_);
}

void Hubbard2::MinimizeE::up_to_down_iter(const size_t irow) {
  const size_t nr = psi_->nrows();
  for (size_t i = irow; i != nr; ++i) {
    pEffTN_ = make_shared<EffectiveTN>(u_term_, pBlocksU_.back(), pBlocksD_.back(), etensor_pool_, i);
    (i%2 == 0) ? left_to_right_iter(i) : right_to_left_iter(i);
    if ( i != nr-1 ) { 
      pBlocksU_.emplace_back(make_shared<Block>(u_term_, pBlocksU_.back(), etensor_pool_, i, compress_map_));
      pBlocksD_.pop_back();
    }
  }
  return;
}

void Hubbard2::MinimizeE::down_to_up_iter(const size_t irow) {
  for (size_t i = irow+1; i-- > 0; ) { 
    pEffTN_ = make_shared<EffectiveTN>(u_term_, pBlocksU_.back(), pBlocksD_.back(), etensor_pool_, i);
    (i%2 == 0) ? right_to_left_iter(i) : left_to_right_iter(i);
    if ( i != 0) {
      pBlocksU_.pop_back();
      pBlocksD_.push_back(make_shared<Block>(u_term_, pBlocksD_.back(), etensor_pool_, i, compress_map_));
    }
  }
  return;
}

void Hubbard2::MinimizeE::left_to_right_iter(const size_t irow) {
  const size_t nc = psi_->ncols();
  for (size_t j = 0; j != nc; ++j) {
    update_one_site(irow, j);
    etensor_pool_.update(irow,j);
  }
  return;
}

void Hubbard2::MinimizeE::right_to_left_iter(const size_t irow) {
  const size_t nc = psi_->ncols();
  const size_t jstart = (irow == nc-1 && nc%2 != 0) ? nc-1u : nc; 
  const size_t jend = (irow == 0u || (irow == nc-1 && nc%2 == 0u)) ? 1u : 0u;  
  for (size_t j = jstart; j-- > jend;) {
    update_one_site(irow, j);
    etensor_pool_.update(irow,j);
  }
  return;
}

void Hubbard2::MinimizeE::update_one_site(const size_t irow, const size_t jcol, const bool update) {
  const pair<Qtensor, Qtensor>& tmp = pEffTN_->evaluate(jcol);
  if (mpi__->rank() == 0) {
    if (update) {
#ifndef NDEBUG
      Qtensor conj = psi_->site(irow,jcol).permute({4,3,2,1,0}, false);
      for (auto&& k : conj.index()) {
        k = k.reverse();
      }
      Qtensor chkE;
      chkE.contract(1.0, conj, {0,1,2,3,4}, tmp.first, {8,6,4,2,0});
      chkE.contract(1.0, chkE, {0,1,2,3,4}, psi_->site(irow,jcol), {0,1,2,3,4});
      cout << " E before: " << chkE.get_block()->element(0,0) << endl;
  
      Qtensor chkN;
      chkN.contract(1.0, conj, {0,1,2,3,4}, tmp.second, {8,6,4,2,0});
      chkN.contract(1.0, chkN, {0,1,2,3,4}, psi_->site(irow,jcol), {0,1,2,3,4});
      cout << " N before: " << chkN.get_block()->element(0,0) << endl;
#endif

      //Cast to Matrix form
      Qtensor tensorH = tmp.first;
      for (auto&& i : tensorH.data()) { //TODO this should be threaded
        Qnum q = i.first[6] + i.first[7] + i.first[8] + i.first[9];
        if (q.parity() == Qnum::ODD) {
          i.second->scale(-1.0);
        }
      }
      tensorH = tensorH.permute({0,2,4,6,8,9,7,5,3,1});
      tensorH = tensorH.permute({0,1,2,3,4,9,8,7,6,5}, false);

      Qtensor tensorN = tmp.second;
      for (auto&& i : tensorN.data()) { //TODO this should be threaded
        Qnum q = i.first[6] + i.first[7] + i.first[8] + i.first[9];
        if (q.parity() == Qnum::ODD) {
          i.second->scale(-1.0);
        }
      }
      tensorN = tensorN.permute({0,2,4,6,8,9,7,5,3,1});
      tensorN = tensorN.permute({0,1,2,3,4,9,8,7,6,5}, false);

      // Diagonalize //
      vector<double> eigenvalues;
      vector<Qtensor> wfns = tensorH.geigs(eigenvalues, tensorN, 1);
      Qtensor out = move(wfns[0]);
      cout << "site: (" << irow << ", "<< jcol << ") :: " << setw(15) << fixed << setprecision(8) << eigenvalues[0] << endl;
      sweep_ene_.push_back(eigenvalues[0]);

#ifndef NDEBUG
      Qtensor conj2 = out.permute({4,3,2,1,0}, false);
      for (auto&& k : conj2.index())
        k = k.reverse();
      chkE.contract(1.0, conj2, {0,1,2,3,4}, tmp.first, {8,6,4,2,0});
      chkE.contract(1.0, chkE, {0,1,2,3,4}, out, {0,1,2,3,4});
      cout << " E after: " << chkE.get_block()->element(0,0) << endl;

      chkN.contract(1.0, conj2, {0,1,2,3,4}, tmp.second, {8,6,4,2,0});
      chkN.contract(1.0, chkN, {0,1,2,3,4}, out, {0,1,2,3,4});
      cout << " N after: " << chkN.get_block()->element(0,0) << endl;
#endif


      // Add Perturbation //
      if (noise_ != 0.0) {
        Qtensor conjugate = out.permute({4,3,2,1,0}, false);
        for (auto&& k : conjugate.index())
          k = k.reverse();
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<double> dis(-noise_, noise_);
        //Ref: A. Gendiar et al.,Prog. Theor. Phys. 110, No. 4 (2003) 691
        Qtensor a, b, c;
        a.contract(1.0, tmp.second, {1,3,5,7,9}, out, {0,1,2,3,4});
        b.contract(1.0, tmp.second, {1,3,5,7,9}, a, {0,1,2,3,4});
        b.contract(1.0, conjugate, {4,3,2,1,0}, b, {0,1,2,3,4});
        c = out;
        c.scale(b.get_block()->element(0,0));
        a -= c;

        //Normalize the perturbed vector
        Qtensor conja = a.permute({4,3,2,1,0}, false);
        for (auto&& k : conja.index()) {
          k = k.reverse();
        }
        Qtensor d;
        d.contract(1.0, conja, {0,1,2,3,4}, tmp.second, {8,6,4,2,0});
        d.contract(1.0, d, {0,1,2,3,4}, a, {0,1,2,3,4});
        a.scale(dis(gen)/sqrt(d.get_block()->element(0,0)));
      
        //Add perturbation
        out += a;
      } 
    psi_->site(irow,jcol) = move(out);
    
    } else {
      
      Qtensor conj = psi_->site(irow,jcol).permute({4,3,2,1,0}, false);
      for (auto&& k : conj.index()) {
        k = k.reverse();
      }   
      Qtensor chkE;
      chkE.contract(1.0, conj, {0,1,2,3,4}, tmp.first, {8,6,4,2,0});
      chkE.contract(1.0, chkE, {0,1,2,3,4}, psi_->site(irow,jcol), {0,1,2,3,4});
  
      Qtensor chkN;
      chkN.contract(1.0, conj, {0,1,2,3,4}, tmp.second, {8,6,4,2,0});
      chkN.contract(1.0, chkN, {0,1,2,3,4}, psi_->site(irow,jcol), {0,1,2,3,4});
   
      const double energy = chkE.get_block()->element(0,0) / chkN.get_block()->element(0,0);
      sweep_ene_.push_back(energy);
      cout << "site: (" << irow << ", "<< jcol << ") :: " << setw(15) << fixed << setprecision(8) << energy << endl;

      // Compute RDM //
      {
        ofstream outRDM;
        outRDM.open("RDM.dump", std::ofstream::out | std::ofstream::app);
      
        if (outRDM.is_open()) {
          const Qtensor& env = pEffTN_->evaluate2(jcol);
          Qtensor rdm;
          rdm.contract(1.0, conj, {0,1,3,4}, env, {6,4,2,0});
          rdm.contract(1.0,  rdm, {1,2,3,4}, psi_->site(irow,jcol), {0,1,3,4});
      
          Matrix dmatrix(4,4);
          Qnum q0(0, true);
          Qnum q1(1, true);
          Qnum q2(2, true);
          dmatrix.copy_block(0,0,1,1, rdm.get_block({q0,q0})->pdata());
          dmatrix.copy_block(1,1,2,2, rdm.get_block({q1,q1})->pdata());
          dmatrix.copy_block(3,3,1,1, rdm.get_block({q2,q2})->pdata());
      
          outRDM << "site: (" << irow << ", "<< jcol << ") \n" << dmatrix << endl;
        } else {
          cout << "Error opening RDM.dump";
        }
        outRDM.close();
      }
    }
  }

  if (update)
    psi_->site(irow,jcol).broadcast_pack(0);
}

double Hubbard2::MinimizeE::check_convergence() {
  const double average_ene = accumulate(sweep_ene_.begin(), sweep_ene_.end(), 0.0) / sweep_ene_.size();
  const double out = abs(last_average_ene_ - average_ene);
  last_average_ene_ = average_ene;
  cout << "  Mean_Energy = " << setw(15) << fixed << setprecision(8) << last_average_ene_ << endl;
  sweep_ene_.clear();
  return out;
}

shared_ptr<PEPS> Hubbard2::MinimizeE::output() const {
  if (conv_) {
    cout << "  *** PEPS MinimizeE Converged ***" << endl;
    cout << "  Average_Energy = " << setw(15) << fixed << setprecision(8) << last_average_ene_ << endl << endl << endl;
    return psi_;
  } else{
    cout << "  *** Failure Convergence PEPS MinimizeE  ***" << endl;
    cout << "  Average_Energy = " << setw(15) << fixed << setprecision(8) << last_average_ene_ << endl << endl << endl;
    return psi_;
  }
}
