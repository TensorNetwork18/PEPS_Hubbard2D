#include <gtest/gtest.h>
#include "qnum.h"

TEST(Qnum, DefaultConstructor) {
  Qnum obj;
  ASSERT_EQ(0,obj.u1());
  ASSERT_EQ(Qnum::EVEN, obj.parity());
}

TEST(Qnum, Constructor) {
  Qnum obj(100,Qnum::ODD);
  ASSERT_EQ(100,obj.u1());
  ASSERT_EQ(Qnum::ODD, obj.parity());
}

TEST(Qnum, CopyConstructor) {
  Qnum obj1(20,Qnum::ODD);
  Qnum obj2(obj1);
  ASSERT_EQ(obj1.u1(), obj2.u1());
  ASSERT_EQ(obj1.parity(), obj2.parity());
  ASSERT_EQ(obj1, obj2);
}

TEST(Qnum, MoveConstructor) {
  Qnum obj1(30,Qnum::ODD);
  Qnum obj2(std::move(obj1));
  ASSERT_EQ(obj1.u1(),0);
  ASSERT_EQ(obj1.parity(),Qnum::EVEN);
  ASSERT_EQ(obj2.u1(),30);
  ASSERT_EQ(obj2.parity(),Qnum::ODD);
}

TEST(Qnum, CopyAssignment) {
  Qnum obj1(20,Qnum::ODD);
  Qnum obj2 = obj1;
  ASSERT_EQ(obj1.u1(), obj2.u1());
  ASSERT_EQ(obj1.parity(), obj2.parity());
  ASSERT_EQ(obj1, obj2);
}

TEST(Qnum, MoveAssignment) {
  Qnum obj1(30,Qnum::ODD);
  Qnum obj2 = std::move(obj1);
  ASSERT_EQ(obj1.u1(),0);
  ASSERT_EQ(obj1.parity(),Qnum::EVEN);
  ASSERT_EQ(obj2.u1(),30);
  ASSERT_EQ(obj2.parity(),Qnum::ODD);
}

TEST(Qnum, Assign) {
  Qnum obj(30, Qnum::ODD);
  obj.assign(10,Qnum::EVEN);
  ASSERT_EQ(obj.u1(),10);
  ASSERT_EQ(obj.parity(),Qnum::EVEN);
}

TEST(Qnum, Comparison) {
  Qnum obj1(30,Qnum::EVEN);
  Qnum obj2(30,Qnum::ODD);
  ASSERT_FALSE((obj1==obj2));
}
