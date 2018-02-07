#include <gtest/gtest.h>
#include "basis.h"


TEST(Basis, DefaultCtr) {
  Basis obj;
  ASSERT_EQ(obj.type(), Basis::IN);
  ASSERT_EQ(obj.total_dim(), 0u);
  ASSERT_TRUE(obj.basis().empty());
}

TEST(Basis, CtrNosym) {
  Basis obj(Basis::IN, 10u);
  Qnum q0 = Qnum::zero();
  ASSERT_EQ(obj.type(), Basis::IN);
  ASSERT_EQ(obj.size(), 1u);
  ASSERT_EQ(obj.total_dim(), 10u);
  ASSERT_EQ(obj.qdim(q0), 10u);
}

TEST(Basis, Ctr) {
  Qnum q_1(-1, Qnum::ODD);
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);

  std::vector<Qnum> qlist {q_1, q_1, q1, q0, q1, q_1};
  Basis obj(Basis::IN, qlist);
  ASSERT_EQ(obj.type(), Basis::IN);
  ASSERT_EQ(obj.total_dim(), 6u);
  ASSERT_EQ(obj.size(), 3u);
  ASSERT_EQ(obj.qdim(q_1), 3u);
  ASSERT_EQ(obj.qdim( q1), 2u);
  ASSERT_EQ(obj.qdim( q0), 1u);
}

TEST(Basis, CopyCtr) {
  Qnum q_1(-1, Qnum::ODD);
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  std::vector<Qnum> qlist {q_1, q_1, q1, q0, q1, q_1};
  Basis obj1(Basis::IN, qlist);
  Basis obj2(obj1);
  ASSERT_EQ(obj1.type(), obj2.type());
  ASSERT_EQ(obj1.basis(), obj2.basis());
}

TEST(Basis, MoveCtr) {
  Qnum q_1(-1, Qnum::ODD);
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  std::vector<Qnum> qlist {q_1, q_1, q1, q0, q1, q_1};
  Basis obj1(Basis::IN, qlist);
  Basis obj2(std::move(obj1));
  ASSERT_EQ(obj1.size(), 0u);
  ASSERT_EQ(obj2.type(), Basis::IN);
  ASSERT_EQ(obj2.total_dim(), 6u);
  ASSERT_EQ(obj2.size(), 3u);
  ASSERT_EQ(obj2.qdim(q_1), 3u);
  ASSERT_EQ(obj2.qdim( q1), 2u);
  ASSERT_EQ(obj2.qdim( q0), 1u);  
}

TEST(Basis, CopyAssign) {
  Qnum q_1(-1, Qnum::ODD);
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  std::vector<Qnum> qlist {q_1, q_1, q1, q0, q1, q_1};
  Basis obj1(Basis::IN, qlist);
  Basis obj2 = obj1;
  ASSERT_EQ(obj1.type(), obj2.type());
  ASSERT_EQ(obj1.basis(), obj2.basis());
}

TEST(Basis, MoveAssign) {
  Qnum q_1(-1, Qnum::ODD);
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  std::vector<Qnum> qlist {q_1, q_1, q1, q0, q1, q_1};
  Basis obj1(Basis::IN, qlist);
  Basis obj2 = std::move(obj1);
  ASSERT_EQ(obj1.size(), 0u);
  ASSERT_EQ(obj2.type(), Basis::IN);
  ASSERT_EQ(obj2.total_dim(), 6u);
  ASSERT_EQ(obj2.size(), 3u);
  ASSERT_EQ(obj2.qdim(q_1), 3u);
  ASSERT_EQ(obj2.qdim( q1), 2u);
  ASSERT_EQ(obj2.qdim( q0), 1u);  
}

TEST(Basis, Reverse) {
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);

  Basis obj(Basis::OUT, {q0,q1,q2});
  obj = obj.reverse();
  ASSERT_EQ(obj.type(), Basis::IN);
  obj = obj.reverse();
  ASSERT_EQ(obj.type(), Basis::OUT);
}

TEST(Basis, Inverse) {
  Qnum q_1(-1, Qnum::ODD);
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  std::vector<Qnum> qlist {q_1, q_1, q1, q0, q1, q_1};
  Basis obj(Basis::IN, qlist);
  auto out = obj.inverse();
  ASSERT_EQ(out.type(), Basis::IN);
  ASSERT_EQ(out.total_dim(), 6u);
  ASSERT_EQ(out.size(), 3u);
  ASSERT_EQ(out.qdim( q1), 3u);
  ASSERT_EQ(out.qdim(q_1), 2u);
  ASSERT_EQ(out.qdim( q0), 1u);
}

TEST(Basis, Fuse) {
  Qnum q_2(-2, Qnum::EVEN);
  Qnum q_1(-1, Qnum::ODD);
  Qnum q0(0, Qnum::EVEN);
  Qnum q1(1, Qnum::ODD);
  Qnum q2(2, Qnum::EVEN);
  std::vector<Qnum> qlist {q_1, q0, q1};
  Basis obj1(Basis::IN, qlist);
  Basis obj2(Basis::IN, qlist);
  auto out = obj1.fuse(obj2);
  ASSERT_EQ(out.type(), Basis::IN); 
  ASSERT_EQ(out.total_dim(), 9u);
  ASSERT_EQ(out.size(), 5u);
  ASSERT_EQ(out.qdim(q_2), 1u);
  ASSERT_EQ(out.qdim(q_1), 2u);
  ASSERT_EQ(out.qdim( q0), 3u);
  ASSERT_EQ(out.qdim( q1), 2u);
  ASSERT_EQ(out.qdim( q2), 1u); 
}
