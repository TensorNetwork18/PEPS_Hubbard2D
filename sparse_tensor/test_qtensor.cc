#include <gtest/gtest.h>
#include "qtensor.h"

TEST(Qtensor, DefaultCtr) {
  Qtensor obj;
  ASSERT_EQ(obj.rank(), 0u);
  ASSERT_EQ(obj.block_size(), 0u);
  ASSERT_EQ(obj.total_element(), 0u);
  ASSERT_EQ(obj.total_size(), 0u);
}

TEST(Qtensor, CtrNosym) {
  //Vector like
  Basis leg0(Basis::IN, 10u);
  std::vector<Basis> lv = {leg0};
  Qtensor obj(lv);
  ASSERT_EQ(obj.rank(), 1u);
  ASSERT_EQ(obj.total_element(), 10u);
  ASSERT_EQ(obj.total_size(), 10u);
}

TEST(Qtensor, CtrNosym2) {
  //Matrix like
  Basis leg0(Basis::IN, 10u);
  Basis leg1(Basis::OUT, 10u);
  Qtensor obj({leg0, leg1});
  ASSERT_EQ(obj.rank(), 2u);
  ASSERT_EQ(obj.total_element(), 100u);
  ASSERT_EQ(obj.total_size(), 100u);
}

TEST(Qtensor, ConstructorSym) {
  //Block Tensor Non-degenerate
  // In (0,1)X(0,1) OUT (0,1,2) 
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  Basis legI(Basis::IN, {q0, q1});
  Basis legO(Basis::OUT, {q0, q1, q2});
  Qtensor obj({legI, legI, legO});
  ASSERT_EQ(obj.rank(), 3u);
  ASSERT_EQ(obj.total_element(), 4u);
  ASSERT_EQ(obj.total_size(), 12u);
}

TEST(Qtensor, ConstructorSym2) {
  //Block Tensor with degeneracy
  // In (0,1)X(0,1) OUT (0,1,1,2) 
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  Basis legI(Basis::IN, {q0, q1});
  Basis legO(Basis::OUT, {q0, q1, q1, q2});
  Qtensor obj({legI, legI, legO});
  ASSERT_EQ(obj.rank(), 3u);
  ASSERT_EQ(obj.total_element(), 6u);
  ASSERT_EQ(obj.total_size(), 16u);
}
