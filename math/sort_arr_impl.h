#ifndef _SORT_ARR_IMPL_H_
#define _SORT_ARR_IMPL_H_

#include <iostream>
#include <algorithm>
#include <memory>
#include <vector>
#include <array>
#include "mkl_extern.h"


//1d
template <size_t i>
void sort_arr_impl(const double* inp, double* out) {
  out[0] = inp[0];
  return;
}

//2d
template <size_t i, size_t j>
void sort_arr_impl(const double* inp, double* out, const size_t d0, const size_t d1) {
  size_t id[2];
  size_t jd[2] = {d0, d1};
  
  long icount = 0;
  for (size_t j1 = 0; j1 != d1; ++j1) {
    id[1] = j1;
    for (size_t j0 = 0; j0 != d0; ++j0, ++icount) {
      id[0] = j0;
      long dist = id[i] + jd[i] * id[j];
      out[dist] = inp[icount];
    }
  }
  return;
}

template <>
void sort_arr_impl<0,1>(const double* inp, double* out, const size_t d0, const size_t d1);

template <>
void sort_arr_impl<1,0>(const double* inp, double* out, const size_t d0, const size_t d1);

//3d
template <size_t i, size_t j, size_t k>
void sort_arr_impl(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2) {
    size_t id[3];
  size_t jd[3] = {d0,d1,d2};

  long icount = 0;
  for (size_t j2 = 0; j2 != d2; ++j2) {
    id[2] = j2; 
    for (size_t j1 = 0; j1 != d1; ++ j1) {
      id[1] = j1; 
      for (size_t j0 = 0; j0 != d0; ++ j0, ++icount) {
        id[0] = j0; 
        long dist = id[i] + jd[i] * (id[j] + jd[j] * id[k]);
        out[dist] = inp[icount];
      }   
    }   
  }
  return;
}

template<>
void sort_arr_impl<0,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2);

template<>
void sort_arr_impl<0,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2);

template<>
void sort_arr_impl<1,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2);

template<>
void sort_arr_impl<1,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2);

template<>
void sort_arr_impl<2,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2);

//4d
template <size_t i, size_t j, size_t k, size_t l>
void sort_arr_impl(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  size_t id[4];
  size_t jd[4] = {d0,d1,d2,d3};

  long icount = 0;
  for (size_t j3 = 0; j3 != d3; ++j3) {
    id[3] = j3;
    for (size_t j2 = 0; j2 != d2; ++j2) {
      id[2] = j2; 
      for (size_t j1 = 0; j1 != d1; ++ j1) {
        id[1] = j1; 
        for (size_t j0 = 0; j0 != d0; ++ j0, ++icount) {
          id[0] = j0; 
          long dist = id[i] + jd[i] * (id[j] + jd[j] * (id[k] + jd[k] * id[l]));
          out[dist] = inp[icount];
        }
      }
    }   
  }
  return;
}

template<>
void sort_arr_impl<0,1,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<0,1,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<0,2,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<0,2,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<0,3,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<0,3,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<1,0,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<1,2,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<1,2,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<2,0,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<2,3,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

template<>
void sort_arr_impl<3,0,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3);

//5d
template <size_t i, size_t j, size_t k, size_t l, size_t m>
void sort_arr_impl(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  size_t id[5];
  size_t jd[5] = {d0,d1,d2,d3,d4};

  long icount = 0;
  for (size_t j4 = 0; j4 != d4; ++j4) {
    id[4] = j4;
    for (size_t j3 = 0; j3 != d3; ++j3) {
      id[3] = j3; 
      for (size_t j2 = 0; j2 != d2; ++j2) {
        id[2] = j2; 
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          id[1] = j1; 
          for (size_t j0 = 0; j0 != d0; ++ j0, ++icount) {
            id[0] = j0; 
            long dist = id[i] + jd[i] * (id[j] + jd[j] * (id[k] + jd[k] * (id[l] + jd[l] * id[m])));
            out[dist] = inp[icount];
          }
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,1,2,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,1,3,2,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,1,3,4,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,1,4,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,1,4,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,2,1,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,2,1,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,2,3,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,2,3,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,2,4,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,2,4,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,3,1,2,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,3,1,4,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,3,2,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,3,2,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,3,4,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,3,4,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,4,1,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,4,1,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,4,2,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,4,2,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,4,3,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<0,4,3,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,0,2,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,0,2,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<1,0,3,2,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,0,3,4,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,0,4,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<1,0,4,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,2,0,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,2,0,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,2,3,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,2,3,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,2,4,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,2,4,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<1,3,2,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,3,2,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,3,4,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,3,4,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,4,0,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<1,4,0,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,4,2,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,4,2,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<1,4,3,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<1,4,3,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<2,0,1,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<2,0,1,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,0,3,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);  

template<>
void sort_arr_impl<2,0,3,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,0,4,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,0,4,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<2,1,0,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,1,0,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,1,3,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,1,3,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,1,4,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,1,4,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<2,3,0,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<2,3,0,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,3,1,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<2,3,1,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<2,3,4,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<2,3,4,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<2,4,0,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,4,0,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,4,1,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,4,1,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<2,4,3,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<2,4,3,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<3,0,1,2,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<3,0,1,4,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<3,0,2,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


//template<>
//void sort_arr_impl<3,0,2,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<3,0,4,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<3,0,4,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


//template<>
//void sort_arr_impl<3,1,0,2,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


//template<>
//void sort_arr_impl<3,1,0,4,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


template<>
void sort_arr_impl<3,1,2,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


template<>
void sort_arr_impl<3,1,2,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


//template<>
//void sort_arr_impl<3,1,4,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


//template<>
//void sort_arr_impl<3,1,4,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


template<>
void sort_arr_impl<3,2,0,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


//template<>
//void sort_arr_impl<3,2,0,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);
//   return;
//}

//template<>
//void sort_arr_impl<3,2,1,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


//template<>
//void sort_arr_impl<3,2,1,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


template<>
void sort_arr_impl<3,2,4,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<3,2,4,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<3,4,0,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<3,4,0,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<3,4,1,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<3,4,1,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<3,4,2,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<3,4,2,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,0,1,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,0,1,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,0,2,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<4,0,2,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,0,3,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<4,0,3,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,1,0,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<4,1,0,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,1,2,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,1,2,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<4,1,3,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


//template<>
//void sort_arr_impl<4,1,3,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


template<>
void sort_arr_impl<4,2,0,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);


//template<>
//void sort_arr_impl<4,2,0,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<4,2,1,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<4,2,1,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,2,3,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,2,3,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,3,0,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<4,3,0,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<4,3,1,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,3,1,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

template<>
void sort_arr_impl<4,3,2,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//template<>
//void sort_arr_impl<4,3,2,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4);

//6d
template <size_t i, size_t j, size_t k, size_t l, size_t m, size_t n>
void sort_arr_impl(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5) {
  size_t id[6];
  size_t jd[6] = {d0,d1,d2,d3,d4,d5};

  long icount = 0;
  for (size_t j5 = 0; j5 != d5; ++j5) {
    id[5] = j5;
    for (size_t j4 = 0; j4 != d4; ++j4) {
      id[4] = j4;
      for (size_t j3 = 0; j3 != d3; ++j3) {
        id[3] = j3;
        for (size_t j2 = 0; j2 != d2; ++j2) {
          id[2] = j2;
          for (size_t j1 = 0; j1 != d1; ++ j1) {
            id[1] = j1;
            for (size_t j0 = 0; j0 != d0; ++ j0, ++icount) {
              id[0] = j0;
              long dist = id[i] + jd[i] * (id[j] + jd[j] * (id[k] + jd[k] * (id[l] + jd[l] * (id[m] + jd[m] * id[n]))));
              out[dist] = inp[icount];
            }
          }
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5);

template<>
void sort_arr_impl<0,2,4,1,3,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5);

template<>
void sort_arr_impl<0,3,1,2,4,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5);

template<>
void sort_arr_impl<0,5,4,3,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5);

template<>
void sort_arr_impl<1,2,4,5,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5);

template<>
void sort_arr_impl<1,2,5,0,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5);

template<>
void sort_arr_impl<2,3,0,1,4,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5);

//7d
template <size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t o>
void sort_arr_impl(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6) {
  size_t id[7];
  size_t jd[7] = {d0,d1,d2,d3,d4,d5,d6};

  long icount = 0;
  for (size_t j6 = 0; j6 != d6; ++j6) {
    id[6] = j6; 
    for (size_t j5 = 0; j5 != d5; ++j5) {
      id[5] = j5;
      for (size_t j4 = 0; j4 != d4; ++j4) {
        id[4] = j4;
        for (size_t j3 = 0; j3 != d3; ++j3) {
          id[3] = j3;
          for (size_t j2 = 0; j2 != d2; ++j2) {
            id[2] = j2;
            for (size_t j1 = 0; j1 != d1; ++ j1) {
              id[1] = j1;
              for (size_t j0 = 0; j0 != d0; ++ j0, ++icount) {
                id[0] = j0;
                long dist = id[i] + jd[i] * (id[j] + jd[j] * (id[k] + jd[k] * (id[l] + jd[l] * (id[m] + jd[m] * (id[n] + jd[n] * id[o])))));
                out[dist] = inp[icount];
              }
            }
          }
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4,5,6>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6);

template<>
void sort_arr_impl<0,1,2,4,6,3,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6);

template<>
void sort_arr_impl<0,2,3,4,5,1,6>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6);

template<>
void sort_arr_impl<0,2,4,6,5,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6);

template<>
void sort_arr_impl<1,3,4,5,0,2,6>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6);

template<>
void sort_arr_impl<2,4,6,0,1,3,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6);

//8d
template <size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t o, size_t p>
void sort_arr_impl(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  size_t id[8];
  size_t jd[8] = {d0,d1,d2,d3,d4,d5,d6,d7};

  long icount = 0;
  for (size_t j7 = 0; j7 != d7; ++j7) {
    id[7] = j7;
    for (size_t j6 = 0; j6 != d6; ++j6) {
      id[6] = j6; 
      for (size_t j5 = 0; j5 != d5; ++j5) {
        id[5] = j5; 
        for (size_t j4 = 0; j4 != d4; ++j4) {
          id[4] = j4; 
          for (size_t j3 = 0; j3 != d3; ++j3) {
            id[3] = j3; 
            for (size_t j2 = 0; j2 != d2; ++j2) {
              id[2] = j2; 
              for (size_t j1 = 0; j1 != d1; ++ j1) {
                id[1] = j1; 
                for (size_t j0 = 0; j0 != d0; ++ j0, ++icount) {
                  id[0] = j0; 
                  long dist = id[i] + jd[i] * (id[j] + jd[j] * (id[k] + jd[k] * (id[l] + jd[l] * (id[m] + jd[m] * (id[n] + jd[n] * (id[o] + jd[o] * id[p]))))));
                  out[dist] = inp[icount];
                }
              }
            }  
          }   
        }   
      }   
    }   
  }
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4,5,6,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

template<>
void sort_arr_impl<0,1,2,3,6,7,4,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

template<>
void sort_arr_impl<0,1,3,4,6,7,2,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

template<>
void sort_arr_impl<0,1,6,2,3,4,5,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

template<>
void sort_arr_impl<0,2,3,4,5,1,6,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

template<>
void sort_arr_impl<0,4,5,3,6,1,2,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

template<>
void sort_arr_impl<2,3,0,1,4,5,6,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

//template<>
//void sort_arr_impl<2,5,3,1,6,0,7,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

template<>
void sort_arr_impl<3,4,2,5,1,6,0,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

template<>
void sort_arr_impl<5,6,0,1,2,3,4,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

template<>
void sort_arr_impl<6,4,2,0,1,3,5,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

//template<>
//void sort_arr_impl<7,5,3,1,0,2,4,6>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7);

//9d
template <size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t o, size_t p, size_t q>
void sort_arr_impl(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8) {
  size_t id[9];
  size_t jd[9] = {d0,d1,d2,d3,d4,d5,d6,d7,d8};

  long icount = 0;
  for (size_t j8 = 0; j8 != d8; ++j8) {
    id[8] = j8;
    for (size_t j7 = 0; j7 != d7; ++j7) {
      id[7] = j7;
      for (size_t j6 = 0; j6 != d6; ++j6) {
        id[6] = j6;
        for (size_t j5 = 0; j5 != d5; ++j5) {
          id[5] = j5;
          for (size_t j4 = 0; j4 != d4; ++j4) {
            id[4] = j4;
            for (size_t j3 = 0; j3 != d3; ++j3) {
              id[3] = j3;
              for (size_t j2 = 0; j2 != d2; ++j2) {
                id[2] = j2;
                for (size_t j1 = 0; j1 != d1; ++ j1) {
                  id[1] = j1;
                  for (size_t j0 = 0; j0 != d0; ++ j0, ++icount) {
                    id[0] = j0;
                    long dist = id[i] + jd[i] * (id[j] + jd[j] * (id[k] + jd[k] * (id[l] + jd[l] * (id[m] + jd[m] * (id[n] + jd[n] * (id[o] + jd[o] * (id[p] + jd[p] * id[q])))))));
                    out[dist] = inp[icount];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4,5,6,7,8>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8);

template<>
void sort_arr_impl<0,2,4,6,8,7,5,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8);

template<>
void sort_arr_impl<8,6,4,2,0,1,3,5,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8);


//10d
template <size_t i, size_t j, size_t k, size_t l, size_t m, size_t n, size_t o, size_t p, size_t q, size_t r>
void sort_arr_impl(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9) {
  size_t id[10];
  size_t jd[10] = {d0,d1,d2,d3,d4,d5,d6,d7,d8,d9};

  long icount = 0;
  for (size_t j9 = 0; j9 != d9; ++j9) {
    id[9] = j9;
    for (size_t j8 = 0; j8 != d8; ++j8) {
      id[8] = j8;
      for (size_t j7 = 0; j7 != d7; ++j7) {
        id[7] = j7;
        for (size_t j6 = 0; j6 != d6; ++j6) {
          id[6] = j6;
          for (size_t j5 = 0; j5 != d5; ++j5) {
            id[5] = j5;
            for (size_t j4 = 0; j4 != d4; ++j4) {
              id[4] = j4;
              for (size_t j3 = 0; j3 != d3; ++j3) {
                id[3] = j3;
                for (size_t j2 = 0; j2 != d2; ++j2) {
                  id[2] = j2;
                  for (size_t j1 = 0; j1 != d1; ++ j1) {
                    id[1] = j1;
                    for (size_t j0 = 0; j0 != d0; ++ j0, ++icount) {
                      id[0] = j0;
                      long dist = id[i] + jd[i] * (id[j] + jd[j] * (id[k] + jd[k] * (id[l] + jd[l] * (id[m] + jd[m] * (id[n] + jd[n] * (id[o] + jd[o] * (id[p] + jd[p] * (id[q] + jd[q] * id[r]))))))));
                      out[dist] = inp[icount];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4,5,6,7,8,9>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9);

template<>
void sort_arr_impl<0,1,2,3,4,9,8,7,6,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9);

template<>
void sort_arr_impl<0,2,4,6,8,1,3,5,7,9>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9);

template<>
void sort_arr_impl<0,2,4,6,8,9,7,5,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9);

template<>
void sort_arr_impl<8,6,4,2,0,1,3,5,7,9>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9);
#endif /* _sort_arr_impl_IMPL_H_  */
