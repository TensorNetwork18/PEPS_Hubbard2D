#include <gtest/gtest.h>
#include "../peps/peps.h"
#include "../peps/tmatrix2.h"


TEST(TMatrix2, TEST) {
  size_t Nr = 3;
  size_t Nc = 3;
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
    bBs[i].insert(Qs[i*Nc]-Qs[2], Dim);
  }
  bBs[Nr].insert(Qs[Nr*Nc]-Qs[2], 1);

  //Qtensors
  std::vector<Qtensor> Qts;
  //First Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], b0, s, bv.reverse(), bhs[j+1].reverse()}, true));
  }
  { //Boundary
    size_t i = 0;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }

  //Middle Row
  for(size_t i=1; i!=Nr-1; ++i){
    for(size_t j=0; j!=Nc-1; ++j){
      Qts.push_back(Qtensor({bhs[j], bv, s, bv.reverse(), bhs[j+1].reverse()}, true));
    }
    { //Boundary
      size_t j = Nc-1;
      Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
    }
  }

  //Last Row
  for(size_t j=0; j!=Nc-1; ++j){
    Qts.push_back(Qtensor({bhs[j], bv, s, b0.reverse(), bhs[j+1].reverse()}, true));
  }
  { //Boundary
    size_t i = Nr-1;
    size_t j = Nc-1;
    Qts.push_back(Qtensor({bhs[j], bBs[i], s, bBs[i+1].reverse(), b0.reverse()}, true));
  }

  PEPS ket(Nr, Nc, Qts);
  ket.randomize();

  std::shared_ptr<Etensor> Norm[Nr][Nc];
  Norm[0][0] = std::make_shared<Etensor>(ket.site(0,0));
  Norm[0][1] = std::make_shared<Etensor>(ket.site(0,1));
  Norm[0][2] = std::make_shared<Etensor>(ket.site(0,2));
  Norm[1][0] = std::make_shared<Etensor>(ket.site(1,0));
  Norm[1][1] = std::make_shared<Etensor>(ket.site(1,1));
  Norm[1][2] = std::make_shared<Etensor>(ket.site(1,2));
  Norm[2][0] = std::make_shared<Etensor>(ket.site(2,0));
  Norm[2][1] = std::make_shared<Etensor>(ket.site(2,1));
  Norm[2][2] = std::make_shared<Etensor>(ket.site(2,2));


  auto R0 = std::make_shared<MPS2>(MPS2(MPS2::U, {Norm[0][0], Norm[0][1], Norm[0][2]}));
  R0 = R0->product({Norm[0][0], Norm[0][1], Norm[0][2]});
  const double scal0 = R0->left_normalize();
  auto R2 = std::make_shared<MPS2>(MPS2(MPS2::D, {Norm[2][0], Norm[2][1], Norm[2][2]}));
  R2 = R2->product({Norm[2][0], Norm[2][1], Norm[2][2]});
  const double scal2 = R2->left_normalize();
  std::vector<std::shared_ptr<const Etensor>> R1 = {Norm[1][0], Norm[1][1], Norm[1][2]};
  TMatrix2 tml(TMatrix2::L, R0, R2, R1);
  TMatrix2 tmr(TMatrix2::R, R0, R2, R1);

  {
    Qtensor Row0;
    Row0.contract(1.0, Norm[0][0]->data(), {6,7}, Norm[0][1]->data(), {0,1});
    Row0.contract(1.0, Row0, {10,11}, Norm[0][2]->data(), {0,1});
    Qtensor Row1;
    Row1.contract(1.0, Norm[1][0]->data(), {6,7}, Norm[1][1]->data(), {0,1});
    Row1.contract(1.0, Row1, {10,11}, Norm[1][2]->data(), {0,1});
    Qtensor Row2;
    Row2.contract(1.0, Norm[2][0]->data(), {6,7}, Norm[2][1]->data(), {0,1});
    Row2.contract(1.0, Row2, {10,11}, Norm[2][2]->data(), {0,1});

    Qtensor ans;
    ans.contract(1.0, Row0, {4,5,8,9,12,13}, Row1, {2,3,6,7,10,11});
    ans.contract(1.0, ans, {12,13,14,15,16,17}, Row2, {2,3,6,7,10,11});
    std::cout << ans << std::endl;
  }

  {
    Qtensor Col0;
    Col0.contract(1.0, Norm[0][0]->data(), {4,5}, Norm[1][0]->data(), {2,3});
    Col0.contract(1.0, Col0, {8,9}, Norm[2][0]->data(), {2,3});
    Qtensor Col1;
    Col1.contract(1.0, Norm[0][1]->data(), {4,5}, Norm[1][1]->data(), {2,3});
    Col1.contract(1.0, Col1, {8,9}, Norm[2][1]->data(), {2,3});
    Qtensor Col2;
    Col2.contract(1.0, Norm[0][2]->data(), {4,5}, Norm[1][2]->data(), {2,3});
    Col2.contract(1.0, Col2, {8,9}, Norm[2][2]->data(), {2,3});

    Qtensor ans;
    ans.contract(1.0, Col0, {4,5,8,9,14,15}, Col1, {0,1,6,7,10,11});
    ans.contract(1.0, ans, {12,13,14,15,18,19}, Col2, {0,1,6,7,10,11});
    std::cout << ans << std::endl;
  }

  // -1 0
  Qtensor norm_l = tml.evaluate(3);
  Qtensor norm_r = tmr.evaluate(0);
  Qtensor norm; 
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;
  // 0 1
  norm_l = tml.evaluate(0);
  norm_r = tmr.evaluate(1);
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;
  // 1 2
  norm_l = tml.evaluate(1);
  norm_r = tmr.evaluate(2);
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;  
  // 2 3
  norm_l = tml.evaluate(2);
  norm_r = tmr.evaluate(3);
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;
}


TEST(TMatrix2, SelfIntegral) { //3x4
  size_t Nr = 3;
  size_t Nc = 4;
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

  {
    Qtensor Row0;
    Row0.contract(1.0, Norm[0][0]->data(), {6,7}, Norm[0][1]->data(), {0,1});
    Row0.contract(1.0, Row0, {10,11}, Norm[0][2]->data(), {0,1});
    Row0.contract(1.0, Row0, {14,15}, Norm[0][3]->data(), {0,1});
    Qtensor Row1;
    Row1.contract(1.0, Norm[1][0]->data(), {6,7}, Norm[1][1]->data(), {0,1});
    Row1.contract(1.0, Row1, {10,11}, Norm[1][2]->data(), {0,1});
    Row1.contract(1.0, Row1, {14,15}, Norm[1][3]->data(), {0,1});
    Qtensor Row2;
    Row2.contract(1.0, Norm[2][0]->data(), {6,7}, Norm[2][1]->data(), {0,1});
    Row2.contract(1.0, Row2, {10,11}, Norm[2][2]->data(), {0,1});
    Row2.contract(1.0, Row2, {14,15}, Norm[2][3]->data(), {0,1});

    Qtensor ans;
    ans.contract(1.0, Row0, {4,5,8,9,12,13,16,17}, Row1, {2,3,6,7,10,11,14,15});
    ans.contract(1.0, ans, {14,15,16,17,18,19,20,21}, Row2, {2,3,6,7,10,11,14,15});
    std::cout << ans << std::endl;
  }

  auto R0 = std::make_shared<MPS2>(MPS2(MPS2::U, {Norm[0][0], Norm[0][1], Norm[0][2], Norm[0][3]}));
  R0 = R0->product({Norm[0][0], Norm[0][1], Norm[0][2], Norm[0][3]});
  const double scal0 = R0->left_normalize();  
  auto R2 = std::make_shared<MPS2>(MPS2(MPS2::D, {Norm[2][0], Norm[2][1], Norm[2][2], Norm[2][3]}));
  R2 = R2->product({Norm[2][0], Norm[2][1], Norm[2][2], Norm[2][3]});
  const double scal2 = R2->left_normalize();
  std::vector<std::shared_ptr<const Etensor>> R1 = {Norm[1][0], Norm[1][1], Norm[1][2], Norm[1][3]};
  TMatrix2 tml(TMatrix2::L, R0, R2, R1);
  TMatrix2 tmr(TMatrix2::R, R0, R2, R1);

  // -1 0
  Qtensor norm_l = tml.evaluate(4);
  Qtensor norm_r = tmr.evaluate(0);
  Qtensor norm;
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;
  // 0 1
  norm_l = tml.evaluate(0);
  norm_r = tmr.evaluate(1);
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;
  // 1 2
  norm_l = tml.evaluate(1);
  norm_r = tmr.evaluate(2);
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;  
  // 2 3
  norm_l = tml.evaluate(2);
  norm_r = tmr.evaluate(3);
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;
  // 3 4
  norm_l = tml.evaluate(3);
  norm_r = tmr.evaluate(4);
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;
}


TEST(TMatrix2, Integral) { //3x3
  size_t Nr = 3;
  size_t Nc = 3;
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

  Qtensor up_dag(Qs[1], {s, s.reverse()}, true);
  up_dag.get_block({Qs[1],Qs[0]})->element(0,0) = 1;
  Qtensor up(-Qs[1], {s, s.reverse()}, true);
  up.get_block({Qs[0],Qs[1]})->element(0,0) = 1;

  Norm[2][0] = std::make_shared<Etensor>(ket.site(2,0), up);
  Norm[1][1] = std::make_shared<Etensor>(ket.site(1,1), up_dag);

  {
    Qtensor Row0;
    Row0.contract(1.0, Norm[0][0]->data(), {6,7}, Norm[0][1]->data(), {0,1});
    Row0.contract(1.0, Row0, {10,11}, Norm[0][2]->data(), {0,1});
    Qtensor Row1;
    Row1.contract(1.0, Norm[1][0]->data(), {6,7}, Norm[1][1]->data(), {0,1});
    Row1.contract(1.0, Row1, {10,11}, Norm[1][2]->data(), {0,1});
    Qtensor Row2;
    Row2.contract(1.0, Norm[2][0]->data(), {6,7}, Norm[2][1]->data(), {0,1});
    Row2.contract(1.0, Row2, {10,11}, Norm[2][2]->data(), {0,1});

    Qtensor ans;
    ans.contract(1.0, Row0, {4,5,8,9,12,13}, Row1, {2,3,6,7,10,11});
    ans.contract(1.0, ans, {12,13,14,15,16,17}, Row2, {2,3,6,7,10,11});
    std::cout << ans << std::endl;
  }

  auto R0 = std::make_shared<MPS2>(MPS2(MPS2::U, {Norm[0][0], Norm[0][1], Norm[0][2]}));
  R0 = R0->product({Norm[0][0], Norm[0][1], Norm[0][2]});
  const double scal0 = R0->left_normalize();
  auto R2 = std::make_shared<MPS2>(MPS2(MPS2::D, {Norm[2][0], Norm[2][1], Norm[2][2]}));
  R2 = R2->product({Norm[2][0], Norm[2][1], Norm[2][2]});
  const double scal2 = R2->left_normalize();
  std::vector<std::shared_ptr<const Etensor>> R1 = {Norm[1][0], Norm[1][1], Norm[1][2]};
  TMatrix2 tml(TMatrix2::L, R0, R2, R1);
  TMatrix2 tmr(TMatrix2::R, R0, R2, R1);

  // -1 0
  Qtensor norm_l = tml.evaluate(3);
  Qtensor norm_r = tmr.evaluate(0);
  Qtensor norm; 
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;
  // 0 1
  norm_l = tml.evaluate(0);
  norm_r = tmr.evaluate(1);
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;
  // 1 2
  norm_l = tml.evaluate(1);
  norm_r = tmr.evaluate(2);
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;  
  // 2 3
  norm_l = tml.evaluate(2);
  norm_r = tmr.evaluate(3);
  norm.contract(scal0*scal2, norm_l, {0,1,2,3}, norm_r, {0,1,2,3});
  std::cout << norm << std::endl;
}
