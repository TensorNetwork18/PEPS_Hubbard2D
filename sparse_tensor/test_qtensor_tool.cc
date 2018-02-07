#include <gtest/gtest.h>
#include "qtensor.h"

TEST(QtensorTool, Plus) {
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  
  Basis legI1(Basis::IN, {q0, q1});
  Basis legI2(Basis::IN, {q0, q1});
  Basis legO1(Basis::OUT, {q0, q1, q2});

  Qtensor obj1({legI1, legI2, legO1}, false);
  auto ptr = std::make_shared<Tensor>(1u, 1u, 1u);
  ptr->fill(-1.0);
  obj1.insert_block({q0,q0,q0}, *ptr);

  Qtensor obj2({legI1, legI2, legO1}, true); 
  obj2.fill(1.0);

  obj1 += obj2;
  ASSERT_DOUBLE_EQ(0.0, obj1.get_block({q0,q0,q0})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, obj1.get_block({q0,q1,q1})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, obj1.get_block({q1,q0,q1})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, obj1.get_block({q1,q1,q2})->element(0,0,0));
}


TEST(QtensorTool, Plus2) {
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  
  Basis legI1(Basis::IN, {q0, q1});
  Basis legI2(Basis::IN, {q0, q1});
  Basis legO1(Basis::OUT, {q0, q1, q2});

  Qtensor obj1({legI1, legI2, legO1}, false);
  auto ptr = std::make_shared<Tensor>(1u, 1u, 1u);
  ptr->fill(-1.0);
  obj1.insert_block({q0,q0,q0}, *ptr);

  Qtensor obj2({legI1, legI2, legO1}, true); 
  obj2.fill(1.0);

  auto out = obj1 + obj2;
  ASSERT_DOUBLE_EQ(0.0, out.get_block({q0,q0,q0})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, out.get_block({q0,q1,q1})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, out.get_block({q1,q0,q1})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, out.get_block({q1,q1,q2})->element(0,0,0));
}

TEST(QtensorTool, Minus) {
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  
  Basis legI1(Basis::IN, {q0, q1});
  Basis legI2(Basis::IN, {q0, q1});
  Basis legO1(Basis::OUT, {q0, q1, q2});

  Qtensor obj1({legI1, legI2, legO1}, false);
  auto ptr = std::make_shared<Tensor>(1u, 1u, 1u);
  ptr->fill(-1.0);
  obj1.insert_block({q0,q0,q0}, *ptr);

  Qtensor obj2({legI1, legI2, legO1}, true); 
  obj2.fill(-1.0);

  obj1 -= obj2;
  ASSERT_DOUBLE_EQ(0.0, obj1.get_block({q0,q0,q0})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, obj1.get_block({q0,q1,q1})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, obj1.get_block({q1,q0,q1})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, obj1.get_block({q1,q1,q2})->element(0,0,0));
}

TEST(QtensorTool, Minus2) {
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  
  Basis legI1(Basis::IN, {q0, q1});
  Basis legI2(Basis::IN, {q0, q1});
  Basis legO1(Basis::OUT, {q0, q1, q2});

  Qtensor obj1({legI1, legI2, legO1}, false);
  auto ptr = std::make_shared<Tensor>(1u, 1u, 1u);
  ptr->fill(-1.0);
  obj1.insert_block({q0,q0,q0}, *ptr);

  Qtensor obj2({legI1, legI2, legO1}, true);
  obj2.fill(-1.0);

  auto out = obj1 - obj2;
  ASSERT_DOUBLE_EQ(0.0, out.get_block({q0,q0,q0})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, out.get_block({q0,q1,q1})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, out.get_block({q1,q0,q1})->element(0,0,0));
  ASSERT_DOUBLE_EQ(1.0, out.get_block({q1,q1,q2})->element(0,0,0));
}


TEST(QtensorTool, Permute) {
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  Basis legI1(Basis::IN, {q0, q1, q1});
  Basis legI2(Basis::IN, {q0, q1});
  Basis legO1(Basis::OUT, {q0, q1});
  Basis legO2(Basis::OUT, {q0, q1, q2, q2});
  Qtensor obj({legI2, legO1, legI1, legO2}); 
  obj.fill(1.0);
  std::cout << obj << std::endl;
  auto out = obj.permute({2,0,1,3}); //(legI1, legI2, legO1, legO2)
  std::cout << out << std::endl;
}

TEST(QtensorTool, SVD) { 
//IN: {0,1}, {0,1} {0,1}
//OUT: {0,1}, {0,1,2}
//Convert T{i0,i1,i2,i3,i4} to M{i0*i1, i2*i3*i4}
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  Basis legI(Basis::IN, {q0, q0, q1, q1});
  Basis legO1(Basis::OUT, {q0, q1});
  Basis legO2(Basis::OUT, {q0, q1, q2});
  //Qtensor obj({legI, legI, legI, legO1, legO2}); 
  //Qtensor obj({legI, legO1, legI, legI, legO2}); 
  Qtensor obj({legI, legI, legO1, legO2}); 
  obj.randomize();
  std::cout << obj << std::endl;

  auto ans = obj.svd(2, Basis::IN, "U");
std::cout << "PASS SVD" << std::endl;
  std::cout << "u:" << std::endl;
  std::cout << ans[0] << std::endl;
  std::cout << "s:" << std::endl;
  std::cout << ans[1] << std::endl;
  std::cout << "vt:" << std::endl;
  std::cout << ans[2] << std::endl;

  //check contract back//
  auto out = make_shared<Qtensor>();
  out->contract(1.0, ans[0], {2}, ans[1], {0});
  out->contract(1.0, *out, {2}, ans[2], {0});
  std::cout << *out << std::endl;
}

TEST(QtensorTool, SVD2) { 
//IN: {0,1}, {0,1} {0,1}
//OUT: {0,1}, {0,1,2}
//Convert T{i0,i1,i2,i3,i4} to M{i0*i1, i2*i3*i4}
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  Basis legI(Basis::IN, {q0, q0, q1, q1});
  Basis legO1(Basis::OUT, {q0, q1});
  Basis legO2(Basis::OUT, {q0, q1, q2});
  //Qtensor obj({legI, legI, legI, legO1, legO2}); 
  //Qtensor obj({legI, legO1, legI, legI, legO2}); 
  Qtensor obj(q1,{legI, legI, legO1, legO2}); 
  obj.randomize();
  std::cout << obj << std::endl;

  auto ans = obj.svd(2, Basis::OUT, "U");
  std::cout << "u:" << std::endl;
  std::cout << ans[0] << std::endl;
  std::cout << "s:" << std::endl;
  std::cout << ans[1] << std::endl;
  std::cout << "vt:" << std::endl;
  std::cout << ans[2] << std::endl;

  //check contract back//
  auto out = make_shared<Qtensor>();
  out->contract(1.0, ans[0], {2}, ans[1], {0});
  out->contract(1.0, *out, {2}, ans[2], {0});
  std::cout << *out << std::endl;
}

TEST(QtensorTool, Contract) {
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  Qnum q3(3, Qnum::ODD);
  Basis legI1(Basis::IN, {q0, q1});
  Basis legI2(Basis::IN, {q0, q1, q1, q2});
  Basis legO1(Basis::OUT, {q0, q1});
  Basis legO2(Basis::OUT, {q0, q1});
  Basis legO3(Basis::OUT, {q0, q1, q2});
  Qtensor obj1({legI1, legI2, legO1, legO2, legO3});
  obj1.fill(1.5);
  std::cout << "obj1:\n" << obj1 << std::endl;

  Basis lI1(Basis::IN, {q0, q1, q2});
  Basis lI2(Basis::IN, {q0, q1});
  Basis lI3(Basis::IN, {q0, q1});
  Basis lO1(Basis::OUT, {q0, q1, q2, q3});
  Qtensor obj2({lI1, lI2, lI3, lO1});
  obj2.fill(0.5);
  std::cout << "obj2:\n" << obj2 << std::endl;

  Qtensor out;
  out.contract(1.0, obj1, {2,3,4}, obj2, {1,2,0});
  std::cout << "out:\n" << out << std::endl;
}

TEST(QtensorTool, Contract2) {
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  Qnum q3(3, Qnum::ODD);
  Basis legI1(Basis::IN, {q0, q1});
  Basis legI2(Basis::IN, {q0, q1, q1, q2});
  Basis legO1(Basis::OUT, {q0, q1});
  Basis legO2(Basis::OUT, {q0, q1});
  Basis legO3(Basis::OUT, {q0, q1, q2});
  Qtensor obj1({legI1, legI2, legO1, legO2, legO3});
  obj1.fill(0.3);
  std::cout << "obj1:\n" << obj1 << std::endl;

  Basis lI1(Basis::IN, {q0, q1, q2});
  Basis lI2(Basis::IN, {q0, q1});
  Basis lI3(Basis::IN, {q0, q1});
  Basis lO1(Basis::OUT, {q0, q1, q2, q3});
  Qtensor obj2({lI1, lI2, lI3, lO1});
  obj2.fill(1.3);
  std::cout << "obj2:\n" << obj2 << std::endl;

  Qtensor out;
  out.contract(1.0, obj1, {}, obj2, {});
  std::cout << "out:\n" << out << std::endl;
}

TEST(QtensorTool, Contract3) {
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  Qnum q3(3, Qnum::ODD);

  Basis legI1(Basis::IN, {q0, q1});
  Basis legI2(Basis::IN, {q0, q1, q1, q2});
  Basis legO1(Basis::OUT, {q0, q1});
  Basis legO2(Basis::OUT, {q0, q1, q1, q2});
  Qtensor obj1({legI1, legI2, legO1, legO2});
  obj1.fill(0.39);
  std::cout << "obj1:\n" << obj1 << std::endl;
  Qtensor out;
  out.contract(1.0, obj1, {0,1,2,3}, obj1, {2,3,0,1});
  std::cout << "out:\n" << out << std::endl;
}
