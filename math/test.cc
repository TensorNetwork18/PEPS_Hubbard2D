#include <gtest/gtest.h>
#include <cmath>
#include "tensor_base.h"
#include "tensor.h"
#include "matrix.h"
#include "vec.h"
#include "mkl_extern.h"
#include "mkl_interface.h"
#include "davidson.h"
#include "sort_arr.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST(Math, Eig_sym) {
  size_t N = 3;
  Matrix A(N,N);
  for(size_t i=0; i<N; ++i){
    for(size_t j=0; j<N; ++j){
      A.element(i,j) = (double)(std::max(i,j)+1);
    }
  }
  Vec ene(N);
  A.eig_sym(ene);
  ASSERT_LE(abs(-1.17761678 - ene(0)), 1.0E-8);
  ASSERT_LE(abs(-0.33892171 - ene(1)), 1.0E-8);
  ASSERT_LE(abs( 7.51653849 - ene(2)), 1.0E-8);
}

TEST(Math, Split) {
  size_t N = 3;
  Matrix A(N,N);
  for(size_t i=0; i<N; ++i){
    for(size_t j=0; j<N; ++j){
      A.element(i,j) = (double)(std::max(i,j)+1);
      if(i!=j) A.element(i,j) *= 0.02;
    }   
  }
  auto out = A.split();
  Matrix C = *out[0] * *out[1];
  std::cout << C << std::endl;
}

TEST(Math, Davidson) {
  size_t N = 3;
  Matrix A(N,N);
  for(size_t i=0; i<N; ++i){
    for(size_t j=0; j<N; ++j){
      A.element(i,j) = (double)(std::max(i,j)+1);
    }
  } 
  auto guess = make_shared<Matrix>(N,2);
  guess->element(0,0) = 1.0;
  guess->element(1,1) = 1.0;
  auto ene = A.eig_davidson(*guess);
  ASSERT_LE(abs(-1.17761678 - ene[0]), 1.0E-8);
  ASSERT_LE(abs(-0.33892171 - ene[1]), 1.0E-8);
}

TEST(Math, Geig_sym) {
  Matrix A(4,4);
  A(0,0) =  0.24;
  A(1,0) =  0.39; A(0,1) = A(1,0);
  A(2,0) =  0.42; A(0,2) = A(2,0);
  A(1,1) = -0.11;
  A(2,1) =  0.79; A(1,2) = A(2,1);
  A(3,1) =  0.63; A(1,3) = A(3,1);
  A(2,2) = -0.25;
  A(3,2) =  0.48; A(2,3) = A(3,2);
  A(3,3) = -0.03;
  Matrix H(A);

  Matrix B(4,4);
  B(0,0) =  2.07;
  B(1,0) =  0.95; B(0,1) = B(1,0);
  B(1,1) =  1.69;
  B(2,1) = -0.29; B(1,2) = B(2,1);
  B(2,2) =  0.65;
  B(3,2) = -0.33; B(2,3) = B(3,2);
  B(3,3) =  1.17;
  Matrix S(B);

  Matrix wfn;
  auto ene = A.geig_sym(B, wfn);
  ASSERT_LE(abs(-0.8305 - ene[0]), 1.0E-4);
  ASSERT_LE(abs(-0.6401 - ene[1]), 1.0E-4);
  ASSERT_LE(abs( 0.0992 - ene[2]), 1.0E-4);
  ASSERT_LE(abs( 1.8525 - ene[3]), 1.0E-4);
}

TEST(Math, Geig_sym2) {
  Matrix A(8,8);
  Matrix B(8,8);
  A(0,0) = 6;
  A(1,1) = 5;
  A(2,2) = 4;
  A(3,3) = 3;
  A(4,4) = 2;
  A(5,5) = 1;

  A(0,6) = 1;
  A(1,7) = 1;
  A(6,0) = 1;
  A(7,1) = 1;
  
  B(0,0) = 1;
  B(1,1) = 1;
  B(2,2) = 1;
  B(3,3) = 1;
  B(4,4) = 1.0E-6;
  B(5,5) = 1.0E-6;
  B(6,6) = 1.0E-6;
  B(7,7) = 1.0E-6;

  Matrix psi;
  auto ene = A.geig_sym(B, psi, 1.0E-6);
  ASSERT_LE(abs( 3.00000 - ene[0]), 1.0E-5);
  ASSERT_LE(abs( 4.00000 - ene[1]), 1.0E-5);
}

TEST(Math, SVD) {
  Matrix A(4,6);
  A(0,0) =  7.52; A(0,1) = -0.76; A(0,2) =  5.13; A(0,3) = -4.75; A(0,4) =  1.33; A(0,5) = -2.40;
  A(1,0) = -1.10; A(1,1) =  0.62; A(1,2) =  6.62; A(1,3) =  8.52; A(1,4) =  4.91; A(1,5) = -6.77;
  A(2,0) = -7.95; A(2,1) =  9.34; A(2,2) = -5.66; A(2,3) =  5.75; A(2,4) = -5.49; A(2,5) =  2.34;
  A(3,0) =  1.08; A(3,1) = -7.10; A(3,2) =  0.87; A(3,3) =  5.30; A(3,4) = -3.52; A(3,5) =  3.95;

  auto ans = A.transpose()->svd();
  ASSERT_LE(abs(18.36597845 - ans[1]->element(0,0)), 1.0E-8);
  ASSERT_LE(abs(13.62997968 - ans[1]->element(1,0)), 1.0E-8);
  ASSERT_LE(abs(10.85333573 - ans[1]->element(2,0)), 1.0E-8);
  ASSERT_LE(abs( 4.49156909 - ans[1]->element(3,0)), 1.0E-8);
  cout << "u\n" << *ans[0] << endl; 
  cout << "vt\n" << *ans[2] << endl; 
}

/*TEST(Math, SVDD) {
  Matrix A(4,6);
  A(0,0) =  7.52; A(0,1) = -0.76; A(0,2) =  5.13; A(0,3) = -4.75; A(0,4) =  1.33; A(0,5) = -2.40;
  A(1,0) = -1.10; A(1,1) =  0.62; A(1,2) =  6.62; A(1,3) =  8.52; A(1,4) =  4.91; A(1,5) = -6.77;
  A(2,0) = -7.95; A(2,1) =  9.34; A(2,2) = -5.66; A(2,3) =  5.75; A(2,4) = -5.49; A(2,5) =  2.34;
  A(3,0) =  1.08; A(3,1) = -7.10; A(3,2) =  0.87; A(3,3) =  5.30; A(3,4) = -3.52; A(3,5) =  3.95;

  auto ans = A.svdd();
  ASSERT_LE(abs(18.36597845 - ans[1]->element(0,0)), 1.0E-8);
  ASSERT_LE(abs(13.62997968 - ans[1]->element(1,0)), 1.0E-8);
  ASSERT_LE(abs(10.85333573 - ans[1]->element(2,0)), 1.0E-8);
  ASSERT_LE(abs( 4.49156909 - ans[1]->element(3,0)), 1.0E-8);
}*/

TEST(Math, JSVD) {
  Matrix A(4,6);
  A(0,0) =  7.52; A(0,1) = -0.76; A(0,2) =  5.13; A(0,3) = -4.75; A(0,4) =  1.33; A(0,5) = -2.40;
  A(1,0) = -1.10; A(1,1) =  0.62; A(1,2) =  6.62; A(1,3) =  8.52; A(1,4) =  4.91; A(1,5) = -6.77;
  A(2,0) = -7.95; A(2,1) =  9.34; A(2,2) = -5.66; A(2,3) =  5.75; A(2,4) = -5.49; A(2,5) =  2.34;
  A(3,0) =  1.08; A(3,1) = -7.10; A(3,2) =  0.87; A(3,3) =  5.30; A(3,4) = -3.52; A(3,5) =  3.95;

  auto ans = A.jsvd();
  ASSERT_LE(abs(18.36597845 - ans[1]->element(0,0)), 1.0E-8);
  ASSERT_LE(abs(13.62997968 - ans[1]->element(1,0)), 1.0E-8);
  ASSERT_LE(abs(10.85333573 - ans[1]->element(2,0)), 1.0E-8);
  ASSERT_LE(abs( 4.49156909 - ans[1]->element(3,0)), 1.0E-8);
  cout << "u\n" << *ans[0] << endl; 
  cout << "vt\n" << *ans[2] << endl; 
}

TEST(Math, SVDJ) {
  Matrix A(4,6);
  A(0,0) =  7.52; A(0,1) = -0.76; A(0,2) =  5.13; A(0,3) = -4.75; A(0,4) =  1.33; A(0,5) = -2.40;
  A(1,0) = -1.10; A(1,1) =  0.62; A(1,2) =  6.62; A(1,3) =  8.52; A(1,4) =  4.91; A(1,5) = -6.77;
  A(2,0) = -7.95; A(2,1) =  9.34; A(2,2) = -5.66; A(2,3) =  5.75; A(2,4) = -5.49; A(2,5) =  2.34;
  A(3,0) =  1.08; A(3,1) = -7.10; A(3,2) =  0.87; A(3,3) =  5.30; A(3,4) = -3.52; A(3,5) =  3.95;

  auto ans = A.svdj();
  ASSERT_LE(abs(18.36597845 - ans[1]->element(0,0)), 1.0E-8);
  ASSERT_LE(abs(13.62997968 - ans[1]->element(1,0)), 1.0E-8);
  ASSERT_LE(abs(10.85333573 - ans[1]->element(2,0)), 1.0E-8);
  ASSERT_LE(abs( 4.49156909 - ans[1]->element(3,0)), 1.0E-8);
  cout << "u\n" << *ans[0] << endl; 
  cout << "vt\n" << *ans[2] << endl; 
}

TEST(Math, QR) {
  Matrix A(3,3);
  A(0,0) = 12; A(0,1) = -51; A(0,2) = 4;
  A(1,0) = 6; A(1,1) = 167; A(1,2) = -68;
  A(2,0) = -4; A(2,1) = 24; A(2,2) = -41;
  A.qr();
  std::cout << A << std::endl;
}


TEST(Math, QR2) {
  Matrix A(3,2);
  A(0,0) = 1; A(0,1) = 2;
  A(1,0) = 3; A(1,1) = 4;
  A(2,0) = 5; A(2,1) = 6;

  A.qr();
  std::cout << A << std::endl;
}

TEST(Math, QRFull) {
  Matrix A(3u,2u);
  A(0,0) = 1; A(0,1) = 2;
  A(1,0) = 3; A(1,1) = 4;
  A(2,0) = 5; A(2,1) = 6;

  auto out = A.qr_full();
  std::cout << *out << std::endl;
  std::cout << A << std::endl;
}

TEST(Math, Sort2) {
  const size_t N = 2;
  std::vector<size_t> pmt(N);
  std::iota(pmt.begin(), pmt.end(), 0);
  Tensor A(3u,4u);
  A.randomize();
  
  do {
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    EXPECT_TRUE(B->is_equal(C));
  } while ( std::next_permutation( pmt.begin(), pmt.end()) );
}

TEST(Math, Sort3) {
  const size_t N = 3;
  std::vector<size_t> pmt(N);
  std::iota(pmt.begin(), pmt.end(), 0); 
  Tensor A(3u,5u,4u);
  A.randomize();
  
  do {
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    EXPECT_TRUE(B->is_equal(C));
  } while ( std::next_permutation( pmt.begin(), pmt.end()) );
}

TEST(Math, Sort4) {
  const size_t N = 4;
  std::vector<size_t> pmt(N);
  std::iota(pmt.begin(), pmt.end(), 0);
  Tensor A(3u,5u,2u,4u);
  A.randomize();

  do {
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    EXPECT_TRUE(B->is_equal(C));
  } while ( std::next_permutation( pmt.begin(), pmt.end()) );
}

TEST(Math, Sort5) {
  const size_t N = 5;
  std::vector<size_t> pmt(N);
  std::iota(pmt.begin(), pmt.end(), 0);
  Tensor A(3u,5u,2u,6u,4u);
  A.randomize();

  do {
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    const bool chk = B->is_equal(C);
    EXPECT_TRUE(chk);
    if (!chk) {
      for (size_t i = 0; i != N; ++i) 
        std::cout << " " << pmt[i];
      std::cout << std::endl;
    }
  } while ( std::next_permutation( pmt.begin(), pmt.end()) );
}

TEST(Math, Sort6) {
  const size_t N = 6;
  std::vector<size_t> pmt(N);
  std::iota(pmt.begin(), pmt.end(), 0); 
  Tensor A(3u,5u,2u,4u,3u,7u);
  A.randomize();

  do {
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    const bool chk = B->is_equal(C);
    EXPECT_TRUE(chk);
    if (!chk) {
      for (size_t i = 0; i != N; ++i) 
        std::cout << " " << pmt[i];
      std::cout << std::endl;
    }   
  } while ( std::next_permutation( pmt.begin(), pmt.end()) );
}

TEST(Math, Sort7) {
  const size_t N = 7;
  std::vector<size_t> pmt(N);
  std::iota(pmt.begin(), pmt.end(), 0); 
  Tensor A(3u,5u,2u,4u,3u,7u,2u);
  A.randomize();

  do {
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    const bool chk = B->is_equal(C);
    EXPECT_TRUE(chk);
    if (!chk) {
      for (size_t i = 0; i != N; ++i) 
        std::cout << " " << pmt[i];
      std::cout << std::endl;
    }   
  } while ( std::next_permutation( pmt.begin(), pmt.end()) );
}

TEST(Math, Sort8) {
  const size_t N = 8;
  std::vector<size_t> pmt(N);
  std::iota(pmt.begin(), pmt.end(), 0); 
  Tensor A(3u,2u,4u,3u,5u,4u,2u,5u);
  A.randomize();

  do {
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    const bool chk = B->is_equal(C);
    EXPECT_TRUE(chk);
    if (!chk) {
      for (size_t i = 0; i != N; ++i) 
        std::cout << " " << pmt[i];
      std::cout << std::endl;
    }   
  } while ( std::next_permutation( pmt.begin(), pmt.end()) );
}

TEST(Math, Sort9) {
  const size_t N = 9; 
  std::vector<size_t> pmt(N);
  std::iota(pmt.begin(), pmt.end(), 0); 
  Tensor A(3u,2u,4u,3u,2u,4u,3u,5u,2u);
  A.randomize();

 {   
    pmt = {0,2,4,6,8,7,5,3,1};
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    const bool chk = B->is_equal(C);
    EXPECT_TRUE(chk);
    if (!chk) {
      for (size_t i = 0; i != N; ++i)
        std::cout << " " << pmt[i];
      std::cout << std::endl;
    }   
  }
}

TEST(Math, Sort10) {
  const size_t N = 10;
  std::vector<size_t> pmt(N);
  std::iota(pmt.begin(), pmt.end(), 0);
  Tensor A(3u,2u,4u,3u,2u,4u,3u,5u,2u,3u);
  A.randomize();

 {   
    pmt = {0,1,2,3,4,9,8,7,6,5};
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    const bool chk = B->is_equal(C);
    EXPECT_TRUE(chk);
    if (!chk) {
      for (size_t i = 0; i != N; ++i)
        std::cout << " " << pmt[i];
      std::cout << std::endl;
    }   
  }

 {   
    pmt = {0,2,4,6,8,1,3,5,7,9};
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    const bool chk = B->is_equal(C);
    EXPECT_TRUE(chk);
    if (!chk) {
      for (size_t i = 0; i != N; ++i)
        std::cout << " " << pmt[i];
      std::cout << std::endl;
    }   
  }

  { 
    pmt = {0,2,4,6,8,9,7,5,3,1};
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    const bool chk = B->is_equal(C);
    EXPECT_TRUE(chk);
    if (!chk) {
      for (size_t i = 0; i != N; ++i)
        std::cout << " " << pmt[i];
      std::cout << std::endl;
    }   
  }

  {   
    pmt = {8,6,4,2,0,1,3,5,7,9};
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    const bool chk = B->is_equal(C);
    EXPECT_TRUE(chk);
    if (!chk) {
      for (size_t i = 0; i != N; ++i)
        std::cout << " " << pmt[i];
      std::cout << std::endl;
    }
  }
  /*do {
    auto B = A.permute(pmt);
    vector<size_t> tmp(N);
    for (size_t i = 0; i != N; ++i) tmp[i] = A.extent(pmt[i]);
    Tensor C(tmp);
    sort_arr<>(A.pdata(), C.pdata(), pmt, A.extent());
    const bool chk = B->is_equal(C);
    EXPECT_TRUE(chk);
    if (!chk) {
      for (size_t i = 0; i != N; ++i)
        std::cout << " " << pmt[i];
      std::cout << std::endl;
    }
  } while ( std::next_permutation( pmt.begin(), pmt.end()) );
  */
}
