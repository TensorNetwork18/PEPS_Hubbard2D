#include <gtest/gtest.h>
#include "../peps/peps.h"
#include "../peps/etensor.h"
#include "../peps/gauge.h"
#include "../peps/hubbard2.h"

TEST(Gauge, RDM_l) {
  size_t Nr = 3;
  size_t Nc = 3;
  size_t Dim = 3;

  std::vector<Qnum> Qs;
  for(size_t i=0; i<=Nc*Nr; ++i){
    Qnum q(i, true);
    Qs.push_back(q);
  }

  //basis
  Basis s(Basis::OUT, {Qs[0], Qs[1]});
  Basis b0(Basis::OUT, {Qs[0]});
  //vertical Basis
  Basis bv(Basis::OUT);
  bv.insert(Qs[0], Dim);
  bv.insert(-Qs[1], Dim);
  bv.insert(Qs[1], Dim);
  //horizontal Basis
  std::vector<Basis> bhs(Nc+1, Basis(Basis::OUT));
  bhs[0].insert(Qs[0], 1);
  for(size_t i=1; i!=Nc; ++i){
    bhs[i].insert(Qs[i+1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[2], Dim);
  }
  bhs[Nc].insert(Qs[0],1);
  //boudary Basis
  std::vector<Basis> bBs(Nr+1, Basis(Basis::OUT));
  bBs[0].insert(Qs[0], 1);
  for(size_t i=1; i!=Nr; ++i){
    bBs[i].insert(Qs[i*Nc]+Qs[1], Dim);
    bBs[i].insert(Qs[i*Nc], Dim);
    bBs[i].insert(Qs[i*Nc]-Qs[1], Dim);
  }
  bBs[Nr].insert(Qs[Nr*Nc-1], 1);

  //Qtensors
  std::vector<Qtensor> Qts;
  //First Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], b0, s, bv.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = 0;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }
  //Middle Row
  for(size_t i=1; i!=Nr-1; ++i){
    for(size_t j=0; j!=Nc-1; ++j){
      Qts.push_back(Qtensor({bhs[j], bv, s, bv.reverse(), bhs[j+1].reverse()}, true));
    }{ //Boundary
      size_t j = Nc-1;
      Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
    }
  }
  //Last Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], bv, s, b0.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = Nr-1;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }

  PEPS ket(Nr, Nc, Qts);
  ket.randomize();
  std::shared_ptr<Etensor> Norm[Nr][Nc];
  for(size_t i=0; i!=Nr; ++i){
    for(size_t j=0; j!=Nc; ++j){
      Norm[i][j] = std::make_shared<Etensor>(ket.site(i,j));
    }
  }

  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }

 //Site 11
 const Qtensor& A11 = Norm[1][1]->data();
 //Trace Left
 Qtensor tr_l({A11.index(0).reverse(), A11.index(1).reverse()}, true);
 for (auto&& e : tr_l.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Up
 Qtensor tr_u({A11.index(2).reverse(), A11.index(3).reverse()}, true);
 for (auto&& e : tr_u.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Dn
 Qtensor tr_d({A11.index(4).reverse(), A11.index(5).reverse()}, true);
 for (auto&& e : tr_d.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Right
 Qtensor tr_r({A11.index(6).reverse(), A11.index(7).reverse()}, true);
 for (auto&& e : tr_r.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }

 //Left DMatrix
 Qtensor DM_l;
 DM_l.contract(1.0, tr_u, {0,1}, A11, {2,3});
 DM_l.contract(1.0, tr_d, {0,1}, DM_l, {2,3});
 DM_l.contract(1.0, tr_r, {0,1}, DM_l, {2,3});
 //std::cout << "DM_l\n" << DM_l << std::endl;

 const auto& ans = DM_l.split(4);
 Qtensor tmp;
 tmp.contract(1.0, ket.site(1,0), {4}, ans[0], {0});
 ket.site(1,0) = move(tmp);
 tmp.contract(1.0, ans[1], {1}, ket.site(1,1), {0});
 ket.site(1,1) = move(tmp);

  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }
}

TEST(Gauge, RDM_u) {
  size_t Nr = 3;
  size_t Nc = 3;
  size_t Dim = 3;

  std::vector<Qnum> Qs;
  for(size_t i=0; i<=Nc*Nr; ++i){
    Qnum q(i, true);
    Qs.push_back(q);
  }

  //basis
  Basis s(Basis::OUT, {Qs[0], Qs[1]});
  Basis b0(Basis::OUT, {Qs[0]});
  //vertical Basis
  Basis bv(Basis::OUT);
  bv.insert(Qs[0], Dim);
  bv.insert(-Qs[1], Dim);
  bv.insert(Qs[1], Dim);
  //horizontal Basis
  std::vector<Basis> bhs(Nc+1, Basis(Basis::OUT));
  bhs[0].insert(Qs[0], 1);
  for(size_t i=1; i!=Nc; ++i){
    bhs[i].insert(Qs[i+1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[2], Dim);
  }
  bhs[Nc].insert(Qs[0],1);
  //boudary Basis
  std::vector<Basis> bBs(Nr+1, Basis(Basis::OUT));
  bBs[0].insert(Qs[0], 1);
  for(size_t i=1; i!=Nr; ++i){
    bBs[i].insert(Qs[i*Nc]+Qs[1], Dim);
    bBs[i].insert(Qs[i*Nc], Dim);
    bBs[i].insert(Qs[i*Nc]-Qs[1], Dim);
  }
  bBs[Nr].insert(Qs[Nr*Nc-1], 1);

  //Qtensors
  std::vector<Qtensor> Qts;
  //First Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], b0, s, bv.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = 0;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }
  //Middle Row
  for(size_t i=1; i!=Nr-1; ++i){
    for(size_t j=0; j!=Nc-1; ++j){
      Qts.push_back(Qtensor({bhs[j], bv, s, bv.reverse(), bhs[j+1].reverse()}, true));
    }{ //Boundary
      size_t j = Nc-1;
      Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
    }
  }
  //Last Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], bv, s, b0.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = Nr-1;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }

  PEPS ket(Nr, Nc, Qts);
  ket.randomize();
  std::shared_ptr<Etensor> Norm[Nr][Nc];
  for(size_t i=0; i!=Nr; ++i){
    for(size_t j=0; j!=Nc; ++j){
      Norm[i][j] = std::make_shared<Etensor>(ket.site(i,j));
    }
  }

  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }


 //Site 11
 const Qtensor& A11 = Norm[1][1]->data();
 //Trace Left
 Qtensor tr_l({A11.index(0).reverse(), A11.index(1).reverse()}, true);
 for (auto&& e : tr_l.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Up
 Qtensor tr_u({A11.index(2).reverse(), A11.index(3).reverse()}, true);
 for (auto&& e : tr_u.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Dn
 Qtensor tr_d({A11.index(4).reverse(), A11.index(5).reverse()}, true);
 for (auto&& e : tr_d.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Right
 Qtensor tr_r({A11.index(6).reverse(), A11.index(7).reverse()}, true);
 for (auto&& e : tr_r.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }

 //Up DMatrix
 Qtensor DM_u;
 DM_u.contract(1.0, tr_l, {0,1}, A11, {0,1});
 DM_u.contract(1.0, tr_d, {0,1}, DM_u, {2,3});
 DM_u.contract(1.0, tr_r, {0,1}, DM_u, {2,3});
 //std::cout << "DM_u\n" << DM_u << std::endl;

 const auto& ans = DM_u.split(4);
 Qtensor x1 = ans[0];
 Qtensor x2 = ans[1];
 Qtensor tmp;
 tmp.contract(1.0, ket.site(0,1), {3}, ans[0], {0});
 ket.site(0,1) = tmp.permute({0,1,2,4,3});
 tmp.contract(1.0, ans[1], {1}, ket.site(1,1), {1});
 ket.site(1,1) = tmp.permute({1,0,2,3,4});
  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }
}


TEST(Gauge, RDM_d) {
  size_t Nr = 3;
  size_t Nc = 3;
  size_t Dim = 3;

  std::vector<Qnum> Qs;
  for(size_t i=0; i<=Nc*Nr; ++i){
    Qnum q(i, true);
    Qs.push_back(q);
  }

  //basis
  Basis s(Basis::OUT, {Qs[0], Qs[1]});
  Basis b0(Basis::OUT, {Qs[0]});
  //vertical Basis
  Basis bv(Basis::OUT);
  bv.insert(Qs[0], Dim);
  bv.insert(-Qs[1], Dim);
  bv.insert(Qs[1], Dim);
  //horizontal Basis
  std::vector<Basis> bhs(Nc+1, Basis(Basis::OUT));
  bhs[0].insert(Qs[0], 1);
  for(size_t i=1; i!=Nc; ++i){
    bhs[i].insert(Qs[i+1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[2], Dim);
  }
  bhs[Nc].insert(Qs[0],1);
  //boudary Basis
  std::vector<Basis> bBs(Nr+1, Basis(Basis::OUT));
  bBs[0].insert(Qs[0], 1);
  for(size_t i=1; i!=Nr; ++i){
    bBs[i].insert(Qs[i*Nc]+Qs[1], Dim);
    bBs[i].insert(Qs[i*Nc], Dim);
    bBs[i].insert(Qs[i*Nc]-Qs[1], Dim);
  }
  bBs[Nr].insert(Qs[Nr*Nc-1], 1);

  //Qtensors
  std::vector<Qtensor> Qts;
  //First Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], b0, s, bv.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = 0;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }
  //Middle Row
  for(size_t i=1; i!=Nr-1; ++i){
    for(size_t j=0; j!=Nc-1; ++j){
      Qts.push_back(Qtensor({bhs[j], bv, s, bv.reverse(), bhs[j+1].reverse()}, true));
    }{ //Boundary
      size_t j = Nc-1;
      Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
    }
  }
  //Last Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], bv, s, b0.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = Nr-1;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }

  PEPS ket(Nr, Nc, Qts);
  ket.randomize();
  std::shared_ptr<Etensor> Norm[Nr][Nc];
  for(size_t i=0; i!=Nr; ++i){
    for(size_t j=0; j!=Nc; ++j){
      Norm[i][j] = std::make_shared<Etensor>(ket.site(i,j));
    }
  }

  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }


 //Site 11
 const Qtensor& A11 = Norm[1][1]->data();
 //Trace Left
 Qtensor tr_l({A11.index(0).reverse(), A11.index(1).reverse()}, true);
 for (auto&& e : tr_l.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Up
 Qtensor tr_u({A11.index(2).reverse(), A11.index(3).reverse()}, true);
 for (auto&& e : tr_u.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Dn
 Qtensor tr_d({A11.index(4).reverse(), A11.index(5).reverse()}, true);
 for (auto&& e : tr_d.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Right
 Qtensor tr_r({A11.index(6).reverse(), A11.index(7).reverse()}, true);
 for (auto&& e : tr_r.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }

 //Down DMatrix
 Qtensor DM_d;
 DM_d.contract(1.0, tr_l, {0,1}, A11, {0,1});
 DM_d.contract(1.0, tr_u, {0,1}, DM_d, {0,1});
 DM_d.contract(1.0, tr_r, {0,1}, DM_d, {2,3});

 const auto& ans = DM_d.split(4);
 Qtensor x1 = ans[0];
 Qtensor x2 = ans[1];
 Qtensor tmp;
 tmp.contract(1.0, ket.site(2,1), {1}, ans[0], {0});
 ket.site(2,1) = tmp.permute({0,4,1,2,3});
 tmp.contract(1.0, ans[1], {1}, ket.site(1,1), {3});
 ket.site(1,1) = tmp.permute({1,2,3,0,4});
  
  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }
}

TEST(Gauge, RDM_r) {
  size_t Nr = 3;
  size_t Nc = 3;
  size_t Dim = 3;

  std::vector<Qnum> Qs;
  for(size_t i=0; i<=Nc*Nr; ++i){
    Qnum q(i, true);
    Qs.push_back(q);
  }

  //basis
  Basis s(Basis::OUT, {Qs[0], Qs[1]});
  Basis b0(Basis::OUT, {Qs[0]});
  //vertical Basis
  Basis bv(Basis::OUT);
  bv.insert(Qs[0], Dim);
  bv.insert(-Qs[1], Dim);
  bv.insert(Qs[1], Dim);
  //horizontal Basis
  std::vector<Basis> bhs(Nc+1, Basis(Basis::OUT));
  bhs[0].insert(Qs[0], 1);
  for(size_t i=1; i!=Nc; ++i){
    bhs[i].insert(Qs[i+1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[2], Dim);
  }
  bhs[Nc].insert(Qs[0],1);
  //boudary Basis
  std::vector<Basis> bBs(Nr+1, Basis(Basis::OUT));
  bBs[0].insert(Qs[0], 1);
  for(size_t i=1; i!=Nr; ++i){
    bBs[i].insert(Qs[i*Nc]+Qs[1], Dim);
    bBs[i].insert(Qs[i*Nc], Dim);
    bBs[i].insert(Qs[i*Nc]-Qs[1], Dim);
  }
  bBs[Nr].insert(Qs[Nr*Nc-1], 1);

  //Qtensors
  std::vector<Qtensor> Qts;
  //First Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], b0, s, bv.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = 0;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }
  //Middle Row
  for(size_t i=1; i!=Nr-1; ++i){
    for(size_t j=0; j!=Nc-1; ++j){
      Qts.push_back(Qtensor({bhs[j], bv, s, bv.reverse(), bhs[j+1].reverse()}, true));
    }{ //Boundary
      size_t j = Nc-1;
      Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
    }
  }
  //Last Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], bv, s, b0.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = Nr-1;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }

  PEPS ket(Nr, Nc, Qts);
  ket.randomize();
  std::shared_ptr<Etensor> Norm[Nr][Nc];
  for(size_t i=0; i!=Nr; ++i){
    for(size_t j=0; j!=Nc; ++j){
      Norm[i][j] = std::make_shared<Etensor>(ket.site(i,j));
    }
  }

  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }


 //Site 11
 const Qtensor& A11 = Norm[1][1]->data();
 //Trace Left
 Qtensor tr_l({A11.index(0).reverse(), A11.index(1).reverse()}, true);
 for (auto&& e : tr_l.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Up
 Qtensor tr_u({A11.index(2).reverse(), A11.index(3).reverse()}, true);
 for (auto&& e : tr_u.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Dn
 Qtensor tr_d({A11.index(4).reverse(), A11.index(5).reverse()}, true);
 for (auto&& e : tr_d.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }
 //Trace Right
 Qtensor tr_r({A11.index(6).reverse(), A11.index(7).reverse()}, true);
 for (auto&& e : tr_r.data()) {
   assert(e.first[0] == e.first[1]);
   assert(e.second->extent(0) == e.second->extent(1));
   size_t dim = e.second->extent(0);
   for (size_t i = 0; i != dim; ++i) {
     e.second->element(i,i) = 1.0;
   }
 }

 //Down DMatrix
 Qtensor DM_r;
 DM_r.contract(1.0, tr_l, {0,1}, A11, {0,1});
 DM_r.contract(1.0, tr_u, {0,1}, DM_r, {0,1});
 DM_r.contract(1.0, tr_d, {0,1}, DM_r, {0,1});

 const auto& ans = DM_r.split(4);
 Qtensor x1 = ans[0];
 Qtensor x2 = ans[1];
 Qtensor tmp;
 tmp.contract(1.0, ans[0], {0}, ket.site(1,2), {0});
 ket.site(1,2) = move(tmp);
 tmp.contract(1.0, ket.site(1,1), {4}, ans[1], {1});
 ket.site(1,1) = move(tmp);
  
  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }
}

TEST(Gauge, LR) {
  size_t Nr = 3;
  size_t Nc = 3;
  size_t Dim = 3;

  std::vector<Qnum> Qs;
  for(size_t i=0; i<=Nc*Nr; ++i){
    Qnum q(i, true);
    Qs.push_back(q);
  }

  //basis
  Basis s(Basis::OUT, {Qs[0], Qs[1]});
  Basis b0(Basis::OUT, {Qs[0]});
  //vertical Basis
  Basis bv(Basis::OUT);
  bv.insert(Qs[0], Dim);
  bv.insert(-Qs[1], Dim);
  bv.insert(Qs[1], Dim);
  //horizontal Basis
  std::vector<Basis> bhs(Nc+1, Basis(Basis::OUT));
  bhs[0].insert(Qs[0], 1); 
  for(size_t i=1; i!=Nc; ++i){
    bhs[i].insert(Qs[i+1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[2], Dim);
  }
  bhs[Nc].insert(Qs[0],1);
  //boudary Basis
  std::vector<Basis> bBs(Nr+1, Basis(Basis::OUT));
  bBs[0].insert(Qs[0], 1); 
  for(size_t i=1; i!=Nr; ++i){
    bBs[i].insert(Qs[i*Nc]+Qs[1], Dim);
    bBs[i].insert(Qs[i*Nc], Dim);
    bBs[i].insert(Qs[i*Nc]-Qs[1], Dim);
  }
  bBs[Nr].insert(Qs[Nr*Nc-1], 1);
  //Qtensors
  std::vector<Qtensor> Qts;
  //First Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], b0, s, bv.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = 0;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }
  //Middle Row
  for(size_t i=1; i!=Nr-1; ++i){
    for(size_t j=0; j!=Nc-1; ++j){
      Qts.push_back(Qtensor({bhs[j], bv, s, bv.reverse(), bhs[j+1].reverse()}, true));
    }{ //Boundary
      size_t j = Nc-1;
      Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
    }
  }
  //Last Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], bv, s, b0.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = Nr-1;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }

  PEPS ket(Nr, Nc, Qts);
  ket.randomize();
  std::shared_ptr<Etensor> Norm[Nr][Nc];
  for(size_t i=0; i!=Nr; ++i){
    for(size_t j=0; j!=Nc; ++j){
      Norm[i][j] = std::make_shared<Etensor>(ket.site(i,j));
    }   
  }

  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }

  //site(0,1) site(11)
  const auto& riR = ket.rdmR(1,0).split(4);
  const auto& riL = ket.rdmL(1,1).split(4);
  Qtensor tmp;
  tmp.contract(1.0, riR[0], {0u}, riL[0], {0u});
  std::vector<double> sgv;
  const auto& ans = tmp.sqrt2(sgv);
  Qtensor gR;
  gR.contract(1.0, riR[1], {0}, ans[0], {0});
  Qtensor gL;
  gL.contract(1.0, ans[1], {1}, riL[1], {0});
  tmp.contract(1.0, ket.site(1,0), {4}, gR, {0});
  ket.site(1,0) = move(tmp);
  tmp.contract(1.0, gL, {1}, ket.site(1,1), {0});
  ket.site(1,1) = move(tmp);
  { //wfn
    std::cout << "AFTER" << std::endl;
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }
}

TEST(Gauge, UD) {
  size_t Nr = 3;
  size_t Nc = 3;
  size_t Dim = 3;

  std::vector<Qnum> Qs;
  for(size_t i=0; i<=Nc*Nr; ++i){
    Qnum q(i, true);
    Qs.push_back(q);
  }

  //basis
  Basis s(Basis::OUT, {Qs[0], Qs[1]});
  Basis b0(Basis::OUT, {Qs[0]});
  //vertical Basis
  Basis bv(Basis::OUT);
  bv.insert(Qs[0], Dim);
  bv.insert(-Qs[1], Dim);
  bv.insert(Qs[1], Dim);
  //horizontal Basis
  std::vector<Basis> bhs(Nc+1, Basis(Basis::OUT));
  bhs[0].insert(Qs[0], 1); 
  for(size_t i=1; i!=Nc; ++i){
    bhs[i].insert(Qs[i+1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[2], Dim);
  }
  bhs[Nc].insert(Qs[0],1);
  //boudary Basis
  std::vector<Basis> bBs(Nr+1, Basis(Basis::OUT));
  bBs[0].insert(Qs[0], 1); 
  for(size_t i=1; i!=Nr; ++i){
    bBs[i].insert(Qs[i*Nc]+Qs[1], Dim);
    bBs[i].insert(Qs[i*Nc], Dim);
    bBs[i].insert(Qs[i*Nc]-Qs[1], Dim);
  }
  bBs[Nr].insert(Qs[Nr*Nc-1], 1);
  //Qtensors
  std::vector<Qtensor> Qts;
  //First Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], b0, s, bv.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = 0;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }
  //Middle Row
  for(size_t i=1; i!=Nr-1; ++i){
    for(size_t j=0; j!=Nc-1; ++j){
      Qts.push_back(Qtensor({bhs[j], bv, s, bv.reverse(), bhs[j+1].reverse()}, true));
    }{ //Boundary
      size_t j = Nc-1;
      Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
    }
  }
  //Last Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], bv, s, b0.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = Nr-1;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }

  PEPS ket(Nr, Nc, Qts);
  ket.randomize();
  std::shared_ptr<Etensor> Norm[Nr][Nc];
  for(size_t i=0; i!=Nr; ++i){
    for(size_t j=0; j!=Nc; ++j){
      Norm[i][j] = std::make_shared<Etensor>(ket.site(i,j));
    }   
  }

  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }

  //site(1,1) site(2,1)
  const auto& riD = ket.rdmD(1,1).split(4);
  const auto& riU = ket.rdmU(2,1).split(4);
  Qtensor tmp;
  tmp.contract(1.0, riD[0], {0u}, riU[0], {0u});
  std::vector<double> sgv;
  const auto& ans = tmp.sqrt2(sgv);
  Qtensor gD;
  gD.contract(1.0, riD[1], {0}, ans[0], {0});
  Qtensor gU;
  gU.contract(1.0, ans[1], {1}, riU[1], {0});

  tmp.contract(1.0, ket.site(1,1), {3}, gD, {0});
  ket.site(1,1) = tmp.permute({0,1,2,4,3});
  tmp.contract(1.0, gU, {1}, ket.site(2,1), {1});
  ket.site(2,1) = tmp.permute({1,0,2,3,4});

  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket.site(0,0), {4}, ket.site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket.site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket.site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket.site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket.site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket.site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket.site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket.site(2,2), {0,1});
    std::cout << psi << std::endl;
  }
}

TEST(Gauge, TEST) {
  size_t Nr = 3;
  size_t Nc = 3;
  size_t Dim = 2;

  std::vector<Qnum> Qs; 
  for(size_t i=0; i<=Nc*Nr; ++i){
    Qnum q(i, true);
    Qs.push_back(q);
  }

  //basis
  Basis s(Basis::OUT, {Qs[0], Qs[1]});
  Basis b0(Basis::OUT, {Qs[0]});
  //vertical Basis
  Basis bv(Basis::OUT);
  bv.insert(Qs[0], Dim);
  bv.insert(-Qs[1], Dim);
  bv.insert(Qs[1], Dim);
  //horizontal Basis
  std::vector<Basis> bhs(Nc+1, Basis(Basis::OUT));
  bhs[0].insert(Qs[0], 1); 
  for(size_t i=1; i!=Nc; ++i){
    bhs[i].insert(Qs[i+1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[2], Dim);
  }
  bhs[Nc].insert(Qs[0],1);
  //boudary Basis
  std::vector<Basis> bBs(Nr+1, Basis(Basis::OUT));
  bBs[0].insert(Qs[0], 1); 
  for(size_t i=1; i!=Nr; ++i){
    bBs[i].insert(Qs[i*Nc]+Qs[1], Dim);
    bBs[i].insert(Qs[i*Nc], Dim);
    bBs[i].insert(Qs[i*Nc]-Qs[1], Dim);
  }
  bBs[Nr].insert(Qs[Nr*Nc-1], 1);

  //Qtensors
  std::vector<Qtensor> Qts;
  //First Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], b0, s, bv.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = 0;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }
  //Middle Row
  for(size_t i=1; i!=Nr-1; ++i){
    for(size_t j=0; j!=Nc-1; ++j){
      Qts.push_back(Qtensor({bhs[j], bv, s, bv.reverse(), bhs[j+1].reverse()}, true));
    }{ //Boundary
      size_t j = Nc-1;
      Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
    }
  }
  //Last Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], bv, s, b0.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = Nr-1;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }

  auto ket = std::make_shared<PEPS>(Nr, Nc, Qts);
  ket->randomize();
  std::shared_ptr<Etensor> Norm[Nr][Nc];
  for(size_t i=0; i!=Nr; ++i){
    for(size_t j=0; j!=Nc; ++j){
      Norm[i][j] = std::make_shared<Etensor>(ket->site(i,j));
    }
  }

  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket->site(0,0), {4}, ket->site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket->site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket->site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket->site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket->site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket->site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket->site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket->site(2,2), {0,1});
    std::cout << psi << std::endl;
  }

  Gauge obj(ket);
  obj.iteration();
  ket = obj.output();

  { //wfn
    Qtensor psi;
    psi.contract(1.0, ket->site(0,0), {4}, ket->site(0,1), {0});
    psi.contract(1.0, psi, {7}, ket->site(0,2), {0});
    psi.contract(1.0, psi, {3}, ket->site(1,0), {1});
    psi.contract(1.0, psi, {13,5}, ket->site(1,1), {0,1});
    psi.contract(1.0, psi, {14,7}, ket->site(1,2), {0,1});
    psi.contract(1.0, psi, {10}, ket->site(2,0), {1});
    psi.contract(1.0, psi, {18,11}, ket->site(2,1), {0,1});
    psi.contract(1.0, psi, {19,12}, ket->site(2,2), {0,1});
    std::cout << psi << std::endl;
  } 
}

TEST(Gauge, TEST2) {
  size_t Nr = 4;
  size_t Nc = 4;
  size_t Dim = 2;

  std::vector<Qnum> Qs; 
  for(size_t i=0; i<=Nc*Nr; ++i){
    Qnum q(i, true);
    Qs.push_back(q);
  }

  //basis
  Basis s(Basis::OUT, {Qs[0], Qs[1], Qs[1], Qs[2]});
  Basis b0(Basis::OUT, {Qs[0]});
  //vertical Basis
  Basis bv(Basis::OUT);
  bv.insert(Qs[0], Dim);
  bv.insert(-Qs[1], Dim);
  bv.insert(Qs[1], Dim);
  //horizontal Basis
  std::vector<Basis> bhs(Nc+1, Basis(Basis::OUT));
  bhs[0].insert(Qs[0], 1); 
  for(size_t i=1; i!=Nc; ++i){
    bhs[i].insert(Qs[i+1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[1], Dim);
    bhs[i].insert(Qs[i+1]-Qs[2], Dim);
  }
  bhs[Nc].insert(Qs[0],1);
  //boudary Basis
  std::vector<Basis> bBs(Nr+1, Basis(Basis::OUT));
  bBs[0].insert(Qs[0], 1); 
  for(size_t i=1; i!=Nr; ++i){
    bBs[i].insert(Qs[i*Nc]+Qs[1], Dim);
    bBs[i].insert(Qs[i*Nc], Dim);
    bBs[i].insert(Qs[i*Nc]-Qs[1], Dim);
  }
  bBs[Nr].insert(Qs[Nr*Nc], 1);

  //Qtensors
  std::vector<Qtensor> Qts;
  //First Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], b0, s, bv.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = 0;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }
  //Middle Row
  for(size_t i=1; i!=Nr-1; ++i){
    for(size_t j=0; j!=Nc-1; ++j){
      Qts.push_back(Qtensor({bhs[j], bv, s, bv.reverse(), bhs[j+1].reverse()}, true));
    }{ //Boundary
      size_t j = Nc-1;
      Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
    }
  }
  //Last Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], bv, s, b0.reverse(), bhs[j+1].reverse()}, true));
  }{ //Boundary
    size_t i = Nr-1;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }

  auto ket = std::make_shared<PEPS>(Nr, Nc, Qts);
  ket->randomize();
  std::shared_ptr<Etensor> Norm[Nr][Nc];
  for(size_t i=0; i!=Nr; ++i){
    for(size_t j=0; j!=Nc; ++j){
      Norm[i][j] = std::make_shared<Etensor>(ket->site(i,j));
    }
  }

  Gauge obj(ket);
  obj.iteration();
  ket = obj.output();
  
  std::vector<Qnum> qd = {-Qs[2], -Qs[1], Qs[0], Qs[1], Qs[2]};
  std::vector<size_t> qm = {2, 4, 6, 4, 2};
  Hubbard2::MinimizeE opt(ket, 0.0, qd, qm); 
  opt.iteration();
  opt.output();
}
