#include <gtest/gtest.h>
#include "../peps/peps.h"
#include "../peps/etensor.h"
#include "../peps/mps2.h"
#include "../peps/compress2.h"

TEST(MPS2, Up) { //5x3 N = 11
  size_t Nr = 3;
  size_t Nc = 5;
  size_t Dim = 1;

  std::vector<Qnum> Qs;
  for(size_t i=0; i<=Nc*Nr; ++i){
    if(i%2==0){
      Qnum q(i,Qnum::EVEN);
      Qs.push_back(q);
    } else {
      Qnum q(i,Qnum::ODD);
      Qs.push_back(q);
    }   
  }
  //basis
  Basis s(Basis::OUT, {Qs[0], Qs[1], Qs[2]});
  Basis b0(Basis::OUT, {Qs[0]});
  //vertical Basis
  Basis bv(Basis::OUT);
  bv.insert(Qs[0], Dim);
  bv.insert(-Qs[1], Dim);
  bv.insert(Qs[1], Dim);
  bv.insert(Qs[2], Dim);
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
    bBs[i].insert(Qs[i*Nc]-Qs[i]+Qs[1], Dim);
    bBs[i].insert(Qs[i*Nc]-Qs[i], Dim);
    bBs[i].insert(Qs[i*Nc]-Qs[i]-Qs[1], Dim);
  }
  //bBs[Nr].insert(Qs[Nr*Nc], 1);
  bBs[Nr].insert(Qs[Nr*Nc]-Qs[1], 1);
  
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

  MPS2 R_1(MPS2::U, {Norm[0][0], Norm[0][1], Norm[0][2], Norm[0][3], Norm[0][4]});
  auto R0 = R_1.product({Norm[0][0], Norm[0][1], Norm[0][2], Norm[0][3], Norm[0][4]});
  double scal = R0->left_normalize();
  ASSERT_DOUBLE_EQ(1.0, R0->left_normalize());

  Qtensor row0;
  row0.contract(1.0, Norm[0][0]->data(), {6,7}, Norm[0][1]->data(), {0,1});
  row0.contract(1.0, row0, {10,11}, Norm[0][2]->data(), {0,1});
  row0.contract(1.0, row0, {14,15}, Norm[0][3]->data(), {0,1});
  row0.contract(1.0, row0, {18,19}, Norm[0][4]->data(), {0,1});
  std::cout << "ans1 " << row0 << std::endl;
  
  Qtensor Row0;
  Row0.contract(1.0, R0->site(0), {3}, R0->site(1), {0});
  Row0.contract(1.0, Row0, {5}, R0->site(2), {0});
  Row0.contract(1.0, Row0, {7}, R0->site(3), {0});
  Row0.contract(1.0, Row0, {9}, R0->site(4), {0});
  Row0.scale(scal);
  std::cout << "ans2 " << Row0 << std::endl;
}


TEST(MPS2, Dn) { //5x3 N = 12
  size_t Nr = 3;
  size_t Nc = 5;
  size_t Dim = 1;

  std::vector<Qnum> Qs;
  for(size_t i=0; i<=Nc*Nr; ++i){
    if(i%2==0){
      Qnum q(i,Qnum::EVEN);
      Qs.push_back(q);
    } else {
      Qnum q(i,Qnum::ODD);
      Qs.push_back(q);
    }   
  }
  //basis
  Basis s(Basis::OUT, {Qs[0], Qs[1], Qs[2]});
  Basis b0(Basis::OUT, {Qs[0]});
  //vertical Basis
  Basis bv(Basis::OUT);
  bv.insert(Qs[0], Dim);
  bv.insert(-Qs[1], Dim);
  bv.insert(Qs[1], Dim);
  bv.insert(Qs[2], Dim);
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
  //bBs[Nr].insert(Qs[Nr*Nc], 1);
  bBs[Nr].insert(Qs[Nr*Nc]-Qs[3], 1);
  
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

  MPS2 R_1(MPS2::D, {Norm[2][0], Norm[2][1], Norm[2][2], Norm[2][3], Norm[2][4]});
  auto R0 = R_1.product({Norm[2][0], Norm[2][1], Norm[2][2], Norm[2][3], Norm[2][4]});
  auto scal = R0->left_normalize();
  ASSERT_DOUBLE_EQ(1.0, R0->left_normalize());

  Qtensor row0;
  row0.contract(1.0, Norm[2][0]->data(), {6,7}, Norm[2][1]->data(), {0,1});
  row0.contract(1.0, row0, {10,11}, Norm[2][2]->data(), {0,1});
  row0.contract(1.0, row0, {14,15}, Norm[2][3]->data(), {0,1});
  row0.contract(1.0, row0, {18,19}, Norm[2][4]->data(), {0,1});
  std::cout << "ans1 " << row0 << std::endl;
  Qtensor Row0;
  Row0.contract(1.0, R0->site(0), {3}, R0->site(1), {0});
  Row0.contract(1.0, Row0, {5}, R0->site(2), {0});
  Row0.contract(1.0, Row0, {7}, R0->site(3), {0});
  Row0.contract(1.0, Row0, {9}, R0->site(4), {0});
  Row0.scale(scal);
  std::cout << "ans2 " << Row0 << std::endl;
}
