#ifndef _SORT_ARR_H_
#define _SORT_ARR_H_

#include <memory>
#include <vector>
#include "sort_arr_impl.h"

template<typename DataType = double>
void sort_arr(const DataType* inp, DataType* out, const std::vector<size_t>& pmt, const std::vector<size_t>& dim) {
  size_t nr = pmt.size();
  assert(dim.size() == nr);

  size_t tag = 0ull;
  for (size_t i = 0; i !=nr ; ++i)
    tag = (tag << 4) + pmt[i];

  switch (tag) {
    case 0ull :
      sort_arr_impl<0>(inp, out); break;
    case 1ull :
      sort_arr_impl<0,1>(inp, out, dim[0], dim[1]); break;
    case 16ull :
      sort_arr_impl<1,0>(inp, out, dim[0], dim[1]); break;
    case 18ull :
      sort_arr_impl<0,1,2>(inp, out, dim[0], dim[1], dim[2]); break;
    case 33ull :
      sort_arr_impl<0,2,1>(inp, out, dim[0], dim[1], dim[2]); break;
    case 258ull :
      sort_arr_impl<1,0,2>(inp, out, dim[0], dim[1], dim[2]); break;
    case 288ull :
      sort_arr_impl<1,2,0>(inp, out, dim[0], dim[1], dim[2]); break;
    case 513ull :
      sort_arr_impl<2,0,1>(inp, out, dim[0], dim[1], dim[2]); break;
    case 528ull :
      sort_arr_impl<2,1,0>(inp, out, dim[0], dim[1], dim[2]); break;
    case 291ull :
      sort_arr_impl<0,1,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 306ull :
      sort_arr_impl<0,1,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 531ull :
      sort_arr_impl<0,2,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 561ull :
      sort_arr_impl<0,2,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 786ull :
      sort_arr_impl<0,3,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 801ull :
      sort_arr_impl<0,3,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 4131ull :
      sort_arr_impl<1,0,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 4146ull :
      sort_arr_impl<1,0,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 4611ull :
      sort_arr_impl<1,2,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 4656ull :
      sort_arr_impl<1,2,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 4866ull :
      sort_arr_impl<1,3,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 4896ull :
      sort_arr_impl<1,3,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 8211ull :
      sort_arr_impl<2,0,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 8241ull :
      sort_arr_impl<2,0,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 8451ull :
      sort_arr_impl<2,1,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 8496ull :
      sort_arr_impl<2,1,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 8961ull :
      sort_arr_impl<2,3,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 8976ull :
      sort_arr_impl<2,3,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 12306ull :
      sort_arr_impl<3,0,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 12321ull :
      sort_arr_impl<3,0,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 12546ull :
      sort_arr_impl<3,1,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 12576ull :
      sort_arr_impl<3,1,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 12801ull :
      sort_arr_impl<3,2,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 12816ull :
      sort_arr_impl<3,2,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3]); break; 
    case 4660ull:
      sort_arr_impl<0,1,2,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 4675ull:
      sort_arr_impl<0,1,2,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 4900ull:
      sort_arr_impl<0,1,3,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 4930ull:
      sort_arr_impl<0,1,3,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 5155ull:
      sort_arr_impl<0,1,4,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 5170ull:
      sort_arr_impl<0,1,4,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 8500ull:
      sort_arr_impl<0,2,1,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 8515ull:
      sort_arr_impl<0,2,1,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 8980ull:
      sort_arr_impl<0,2,3,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 9025ull:
      sort_arr_impl<0,2,3,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 9235ull:
      sort_arr_impl<0,2,4,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 9265ull:
      sort_arr_impl<0,2,4,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 12580ull:
      sort_arr_impl<0,3,1,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 12610ull:
      sort_arr_impl<0,3,1,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 12820ull:
      sort_arr_impl<0,3,2,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 12865ull:
      sort_arr_impl<0,3,2,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 13330ull:
      sort_arr_impl<0,3,4,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 13345ull:
      sort_arr_impl<0,3,4,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 16675ull:
      sort_arr_impl<0,4,1,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 16690ull:
      sort_arr_impl<0,4,1,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 16915ull:
      sort_arr_impl<0,4,2,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 16945ull:
      sort_arr_impl<0,4,2,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 17170ull:
      sort_arr_impl<0,4,3,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 17185ull:
      sort_arr_impl<0,4,3,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 66100ull:
      sort_arr_impl<1,0,2,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 66115ull:
      sort_arr_impl<1,0,2,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 66340ull:
      sort_arr_impl<1,0,3,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 66370ull:
      sort_arr_impl<1,0,3,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 66595ull:
      sort_arr_impl<1,0,4,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 66610ull:
      sort_arr_impl<1,0,4,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 73780ull:
      sort_arr_impl<1,2,0,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 73795ull:
      sort_arr_impl<1,2,0,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 74500ull:
      sort_arr_impl<1,2,3,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 74560ull:
      sort_arr_impl<1,2,3,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 74755ull:
      sort_arr_impl<1,2,4,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 74800ull:
      sort_arr_impl<1,2,4,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 77860ull:
      sort_arr_impl<1,3,0,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 77890ull:
      sort_arr_impl<1,3,0,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 78340ull:
      sort_arr_impl<1,3,2,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 78400ull:
      sort_arr_impl<1,3,2,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 78850ull:
      sort_arr_impl<1,3,4,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 78880ull:
      sort_arr_impl<1,3,4,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 81955ull:
      sort_arr_impl<1,4,0,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 81970ull:
      sort_arr_impl<1,4,0,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 82435ull:
      sort_arr_impl<1,4,2,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 82480ull:
      sort_arr_impl<1,4,2,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 82690ull:
      sort_arr_impl<1,4,3,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 82720ull:
      sort_arr_impl<1,4,3,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 131380ull:
      sort_arr_impl<2,0,1,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 131395ull:
      sort_arr_impl<2,0,1,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 131860ull:
      sort_arr_impl<2,0,3,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 131905ull:
      sort_arr_impl<2,0,3,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 132115ull:
      sort_arr_impl<2,0,4,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 132145ull:
      sort_arr_impl<2,0,4,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 135220ull:
      sort_arr_impl<2,1,0,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 135235ull:
      sort_arr_impl<2,1,0,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 135940ull:
      sort_arr_impl<2,1,3,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 136000ull:
      sort_arr_impl<2,1,3,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 136195ull:
      sort_arr_impl<2,1,4,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 136240ull:
      sort_arr_impl<2,1,4,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 143380ull:
      sort_arr_impl<2,3,0,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 143425ull:
      sort_arr_impl<2,3,0,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 143620ull:
      sort_arr_impl<2,3,1,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 143680ull:
      sort_arr_impl<2,3,1,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 144385ull:
      sort_arr_impl<2,3,4,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 144400ull:
      sort_arr_impl<2,3,4,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 147475ull:
      sort_arr_impl<2,4,0,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 147505ull:
      sort_arr_impl<2,4,0,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 147715ull:
      sort_arr_impl<2,4,1,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 147760ull:
      sort_arr_impl<2,4,1,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 148225ull:
      sort_arr_impl<2,4,3,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 148240ull:
      sort_arr_impl<2,4,3,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 196900ull:
      sort_arr_impl<3,0,1,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 196930ull:
      sort_arr_impl<3,0,1,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 197140ull:
      sort_arr_impl<3,0,2,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 197185ull:
      sort_arr_impl<3,0,2,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 197650ull:
      sort_arr_impl<3,0,4,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 197665ull:
      sort_arr_impl<3,0,4,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 200740ull:
      sort_arr_impl<3,1,0,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 200770ull:
      sort_arr_impl<3,1,0,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 201220ull:
      sort_arr_impl<3,1,2,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 201280ull:
      sort_arr_impl<3,1,2,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 201730ull:
      sort_arr_impl<3,1,4,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 201760ull:
      sort_arr_impl<3,1,4,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 204820ull:
      sort_arr_impl<3,2,0,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 204865ull:
      sort_arr_impl<3,2,0,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 205060ull:
      sort_arr_impl<3,2,1,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 205120ull:
      sort_arr_impl<3,2,1,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 205825ull:
      sort_arr_impl<3,2,4,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 205840ull:
      sort_arr_impl<3,2,4,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 213010ull:
      sort_arr_impl<3,4,0,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 213025ull:
      sort_arr_impl<3,4,0,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 213250ull:
      sort_arr_impl<3,4,1,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 213280ull:
      sort_arr_impl<3,4,1,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 213505ull:
      sort_arr_impl<3,4,2,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 213520ull:
      sort_arr_impl<3,4,2,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 262435ull:
      sort_arr_impl<4,0,1,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 262450ull:
      sort_arr_impl<4,0,1,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 262675ull:
      sort_arr_impl<4,0,2,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 262705ull:
      sort_arr_impl<4,0,2,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 262930ull:
      sort_arr_impl<4,0,3,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 262945ull:
      sort_arr_impl<4,0,3,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 266275ull:
      sort_arr_impl<4,1,0,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 266290ull:
      sort_arr_impl<4,1,0,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 266755ull:
      sort_arr_impl<4,1,2,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 266800ull:
      sort_arr_impl<4,1,2,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 267010ull:
      sort_arr_impl<4,1,3,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 267040ull:
      sort_arr_impl<4,1,3,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 270355ull:
      sort_arr_impl<4,2,0,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 270385ull:
      sort_arr_impl<4,2,0,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 270595ull:
      sort_arr_impl<4,2,1,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 270640ull:
      sort_arr_impl<4,2,1,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 271105ull:
      sort_arr_impl<4,2,3,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 271120ull:
      sort_arr_impl<4,2,3,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 274450ull:
      sort_arr_impl<4,3,0,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 274465ull:
      sort_arr_impl<4,3,0,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 274690ull:
      sort_arr_impl<4,3,1,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 274720ull:
      sort_arr_impl<4,3,1,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 274945ull:
      sort_arr_impl<4,3,2,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 274960ull:
      sort_arr_impl<4,3,2,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4]); break;
    case 74565ull:
      sort_arr_impl<0,1,2,3,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 74580ull:
      sort_arr_impl<0,1,2,3,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 74805ull:
      sort_arr_impl<0,1,2,4,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 74835ull:
      sort_arr_impl<0,1,2,4,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 75060ull:
      sort_arr_impl<0,1,2,5,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 75075ull:
      sort_arr_impl<0,1,2,5,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 78405ull:
      sort_arr_impl<0,1,3,2,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 78420ull:
      sort_arr_impl<0,1,3,2,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 78885ull:
      sort_arr_impl<0,1,3,4,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 78930ull:
      sort_arr_impl<0,1,3,4,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 79140ull:
      sort_arr_impl<0,1,3,5,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 79170ull:
      sort_arr_impl<0,1,3,5,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 82485ull:
      sort_arr_impl<0,1,4,2,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 82515ull:
      sort_arr_impl<0,1,4,2,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 82725ull:
      sort_arr_impl<0,1,4,3,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 82770ull:
      sort_arr_impl<0,1,4,3,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 83235ull:
      sort_arr_impl<0,1,4,5,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 83250ull:
      sort_arr_impl<0,1,4,5,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 86580ull:
      sort_arr_impl<0,1,5,2,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 86595ull:
      sort_arr_impl<0,1,5,2,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 86820ull:
      sort_arr_impl<0,1,5,3,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 86850ull:
      sort_arr_impl<0,1,5,3,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 87075ull:
      sort_arr_impl<0,1,5,4,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 87090ull:
      sort_arr_impl<0,1,5,4,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 136005ull:
      sort_arr_impl<0,2,1,3,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 136020ull:
      sort_arr_impl<0,2,1,3,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 136245ull:
      sort_arr_impl<0,2,1,4,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 136275ull:
      sort_arr_impl<0,2,1,4,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 136500ull:
      sort_arr_impl<0,2,1,5,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 136515ull:
      sort_arr_impl<0,2,1,5,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 143685ull:
      sort_arr_impl<0,2,3,1,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 143700ull:
      sort_arr_impl<0,2,3,1,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 144405ull:
      sort_arr_impl<0,2,3,4,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 144465ull:
      sort_arr_impl<0,2,3,4,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 144660ull:
      sort_arr_impl<0,2,3,5,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 144705ull:
      sort_arr_impl<0,2,3,5,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 147765ull:
      sort_arr_impl<0,2,4,1,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 147795ull:
      sort_arr_impl<0,2,4,1,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 148245ull:
      sort_arr_impl<0,2,4,3,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 148305ull:
      sort_arr_impl<0,2,4,3,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 148755ull:
      sort_arr_impl<0,2,4,5,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 148785ull:
      sort_arr_impl<0,2,4,5,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 151860ull:
      sort_arr_impl<0,2,5,1,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 151875ull:
      sort_arr_impl<0,2,5,1,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 152340ull:
      sort_arr_impl<0,2,5,3,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 152385ull:
      sort_arr_impl<0,2,5,3,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 152595ull:
      sort_arr_impl<0,2,5,4,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 152625ull:
      sort_arr_impl<0,2,5,4,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 201285ull:
      sort_arr_impl<0,3,1,2,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 201300ull:
      sort_arr_impl<0,3,1,2,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 201765ull:
      sort_arr_impl<0,3,1,4,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 201810ull:
      sort_arr_impl<0,3,1,4,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 202020ull:
      sort_arr_impl<0,3,1,5,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 202050ull:
      sort_arr_impl<0,3,1,5,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 205125ull:
      sort_arr_impl<0,3,2,1,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 205140ull:
      sort_arr_impl<0,3,2,1,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 205845ull:
      sort_arr_impl<0,3,2,4,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 205905ull:
      sort_arr_impl<0,3,2,4,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 206100ull:
      sort_arr_impl<0,3,2,5,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 206145ull:
      sort_arr_impl<0,3,2,5,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 213285ull:
      sort_arr_impl<0,3,4,1,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 213330ull:
      sort_arr_impl<0,3,4,1,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 213525ull:
      sort_arr_impl<0,3,4,2,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 213585ull:
      sort_arr_impl<0,3,4,2,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 214290ull:
      sort_arr_impl<0,3,4,5,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 214305ull:
      sort_arr_impl<0,3,4,5,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 217380ull:
      sort_arr_impl<0,3,5,1,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 217410ull:
      sort_arr_impl<0,3,5,1,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 217620ull:
      sort_arr_impl<0,3,5,2,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 217665ull:
      sort_arr_impl<0,3,5,2,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 218130ull:
      sort_arr_impl<0,3,5,4,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 218145ull:
      sort_arr_impl<0,3,5,4,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 266805ull:
      sort_arr_impl<0,4,1,2,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 266835ull:
      sort_arr_impl<0,4,1,2,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 267045ull:
      sort_arr_impl<0,4,1,3,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 267090ull:
      sort_arr_impl<0,4,1,3,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 267555ull:
      sort_arr_impl<0,4,1,5,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 267570ull:
      sort_arr_impl<0,4,1,5,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 270645ull:
      sort_arr_impl<0,4,2,1,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 270675ull:
      sort_arr_impl<0,4,2,1,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 271125ull:
      sort_arr_impl<0,4,2,3,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 271185ull:
      sort_arr_impl<0,4,2,3,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 271635ull:
      sort_arr_impl<0,4,2,5,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 271665ull:
      sort_arr_impl<0,4,2,5,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 274725ull:
      sort_arr_impl<0,4,3,1,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 274770ull:
      sort_arr_impl<0,4,3,1,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 274965ull:
      sort_arr_impl<0,4,3,2,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 275025ull:
      sort_arr_impl<0,4,3,2,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 275730ull:
      sort_arr_impl<0,4,3,5,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 275745ull:
      sort_arr_impl<0,4,3,5,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 282915ull:
      sort_arr_impl<0,4,5,1,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 282930ull:
      sort_arr_impl<0,4,5,1,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 283155ull:
      sort_arr_impl<0,4,5,2,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 283185ull:
      sort_arr_impl<0,4,5,2,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 283410ull:
      sort_arr_impl<0,4,5,3,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 283425ull:
      sort_arr_impl<0,4,5,3,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 332340ull:
      sort_arr_impl<0,5,1,2,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 332355ull:
      sort_arr_impl<0,5,1,2,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 332580ull:
      sort_arr_impl<0,5,1,3,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 332610ull:
      sort_arr_impl<0,5,1,3,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 332835ull:
      sort_arr_impl<0,5,1,4,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 332850ull:
      sort_arr_impl<0,5,1,4,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 336180ull:
      sort_arr_impl<0,5,2,1,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 336195ull:
      sort_arr_impl<0,5,2,1,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 336660ull:
      sort_arr_impl<0,5,2,3,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 336705ull:
      sort_arr_impl<0,5,2,3,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 336915ull:
      sort_arr_impl<0,5,2,4,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 336945ull:
      sort_arr_impl<0,5,2,4,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 340260ull:
      sort_arr_impl<0,5,3,1,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 340290ull:
      sort_arr_impl<0,5,3,1,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 340500ull:
      sort_arr_impl<0,5,3,2,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 340545ull:
      sort_arr_impl<0,5,3,2,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 341010ull:
      sort_arr_impl<0,5,3,4,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 341025ull:
      sort_arr_impl<0,5,3,4,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 344355ull:
      sort_arr_impl<0,5,4,1,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 344370ull:
      sort_arr_impl<0,5,4,1,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 344595ull:
      sort_arr_impl<0,5,4,2,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 344625ull:
      sort_arr_impl<0,5,4,2,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 344850ull:
      sort_arr_impl<0,5,4,3,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 344865ull:
      sort_arr_impl<0,5,4,3,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1057605ull:
      sort_arr_impl<1,0,2,3,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1057620ull:
      sort_arr_impl<1,0,2,3,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1057845ull:
      sort_arr_impl<1,0,2,4,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1057875ull:
      sort_arr_impl<1,0,2,4,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1058100ull:
      sort_arr_impl<1,0,2,5,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1058115ull:
      sort_arr_impl<1,0,2,5,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1061445ull:
      sort_arr_impl<1,0,3,2,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1061460ull:
      sort_arr_impl<1,0,3,2,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1061925ull:
      sort_arr_impl<1,0,3,4,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1061970ull:
      sort_arr_impl<1,0,3,4,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1062180ull:
      sort_arr_impl<1,0,3,5,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1062210ull:
      sort_arr_impl<1,0,3,5,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1065525ull:
      sort_arr_impl<1,0,4,2,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1065555ull:
      sort_arr_impl<1,0,4,2,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1065765ull:
      sort_arr_impl<1,0,4,3,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1065810ull:
      sort_arr_impl<1,0,4,3,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1066275ull:
      sort_arr_impl<1,0,4,5,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1066290ull:
      sort_arr_impl<1,0,4,5,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1069620ull:
      sort_arr_impl<1,0,5,2,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1069635ull:
      sort_arr_impl<1,0,5,2,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1069860ull:
      sort_arr_impl<1,0,5,3,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1069890ull:
      sort_arr_impl<1,0,5,3,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1070115ull:
      sort_arr_impl<1,0,5,4,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1070130ull:
      sort_arr_impl<1,0,5,4,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1180485ull:
      sort_arr_impl<1,2,0,3,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1180500ull:
      sort_arr_impl<1,2,0,3,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1180725ull:
      sort_arr_impl<1,2,0,4,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1180755ull:
      sort_arr_impl<1,2,0,4,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1180980ull:
      sort_arr_impl<1,2,0,5,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1180995ull:
      sort_arr_impl<1,2,0,5,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1192005ull:
      sort_arr_impl<1,2,3,0,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1192020ull:
      sort_arr_impl<1,2,3,0,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1192965ull:
      sort_arr_impl<1,2,3,4,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1193040ull:
      sort_arr_impl<1,2,3,4,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1193220ull:
      sort_arr_impl<1,2,3,5,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1193280ull:
      sort_arr_impl<1,2,3,5,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1196085ull:
      sort_arr_impl<1,2,4,0,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1196115ull:
      sort_arr_impl<1,2,4,0,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1196805ull:
      sort_arr_impl<1,2,4,3,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1196880ull:
      sort_arr_impl<1,2,4,3,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1197315ull:
      sort_arr_impl<1,2,4,5,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1197360ull:
      sort_arr_impl<1,2,4,5,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1200180ull:
      sort_arr_impl<1,2,5,0,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1200195ull:
      sort_arr_impl<1,2,5,0,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1200900ull:
      sort_arr_impl<1,2,5,3,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1200960ull:
      sort_arr_impl<1,2,5,3,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1201155ull:
      sort_arr_impl<1,2,5,4,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1201200ull:
      sort_arr_impl<1,2,5,4,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1245765ull:
      sort_arr_impl<1,3,0,2,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1245780ull:
      sort_arr_impl<1,3,0,2,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1246245ull:
      sort_arr_impl<1,3,0,4,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1246290ull:
      sort_arr_impl<1,3,0,4,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1246500ull:
      sort_arr_impl<1,3,0,5,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1246530ull:
      sort_arr_impl<1,3,0,5,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1253445ull:
      sort_arr_impl<1,3,2,0,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1253460ull:
      sort_arr_impl<1,3,2,0,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1254405ull:
      sort_arr_impl<1,3,2,4,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1254480ull:
      sort_arr_impl<1,3,2,4,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1254660ull:
      sort_arr_impl<1,3,2,5,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1254720ull:
      sort_arr_impl<1,3,2,5,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1261605ull:
      sort_arr_impl<1,3,4,0,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1261650ull:
      sort_arr_impl<1,3,4,0,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1262085ull:
      sort_arr_impl<1,3,4,2,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1262160ull:
      sort_arr_impl<1,3,4,2,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1262850ull:
      sort_arr_impl<1,3,4,5,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1262880ull:
      sort_arr_impl<1,3,4,5,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1265700ull:
      sort_arr_impl<1,3,5,0,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1265730ull:
      sort_arr_impl<1,3,5,0,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1266180ull:
      sort_arr_impl<1,3,5,2,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1266240ull:
      sort_arr_impl<1,3,5,2,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1266690ull:
      sort_arr_impl<1,3,5,4,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1266720ull:
      sort_arr_impl<1,3,5,4,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1311285ull:
      sort_arr_impl<1,4,0,2,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1311315ull:
      sort_arr_impl<1,4,0,2,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1311525ull:
      sort_arr_impl<1,4,0,3,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1311570ull:
      sort_arr_impl<1,4,0,3,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1312035ull:
      sort_arr_impl<1,4,0,5,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1312050ull:
      sort_arr_impl<1,4,0,5,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1318965ull:
      sort_arr_impl<1,4,2,0,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1318995ull:
      sort_arr_impl<1,4,2,0,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1319685ull:
      sort_arr_impl<1,4,2,3,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1319760ull:
      sort_arr_impl<1,4,2,3,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1320195ull:
      sort_arr_impl<1,4,2,5,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1320240ull:
      sort_arr_impl<1,4,2,5,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1323045ull:
      sort_arr_impl<1,4,3,0,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1323090ull:
      sort_arr_impl<1,4,3,0,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1323525ull:
      sort_arr_impl<1,4,3,2,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1323600ull:
      sort_arr_impl<1,4,3,2,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1324290ull:
      sort_arr_impl<1,4,3,5,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1324320ull:
      sort_arr_impl<1,4,3,5,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1331235ull:
      sort_arr_impl<1,4,5,0,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1331250ull:
      sort_arr_impl<1,4,5,0,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1331715ull:
      sort_arr_impl<1,4,5,2,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1331760ull:
      sort_arr_impl<1,4,5,2,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1331970ull:
      sort_arr_impl<1,4,5,3,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1332000ull:
      sort_arr_impl<1,4,5,3,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1376820ull:
      sort_arr_impl<1,5,0,2,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1376835ull:
      sort_arr_impl<1,5,0,2,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1377060ull:
      sort_arr_impl<1,5,0,3,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1377090ull:
      sort_arr_impl<1,5,0,3,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1377315ull:
      sort_arr_impl<1,5,0,4,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1377330ull:
      sort_arr_impl<1,5,0,4,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1384500ull:
      sort_arr_impl<1,5,2,0,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1384515ull:
      sort_arr_impl<1,5,2,0,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1385220ull:
      sort_arr_impl<1,5,2,3,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1385280ull:
      sort_arr_impl<1,5,2,3,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1385475ull:
      sort_arr_impl<1,5,2,4,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1385520ull:
      sort_arr_impl<1,5,2,4,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1388580ull:
      sort_arr_impl<1,5,3,0,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1388610ull:
      sort_arr_impl<1,5,3,0,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1389060ull:
      sort_arr_impl<1,5,3,2,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1389120ull:
      sort_arr_impl<1,5,3,2,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1389570ull:
      sort_arr_impl<1,5,3,4,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1389600ull:
      sort_arr_impl<1,5,3,4,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1392675ull:
      sort_arr_impl<1,5,4,0,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1392690ull:
      sort_arr_impl<1,5,4,0,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1393155ull:
      sort_arr_impl<1,5,4,2,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1393200ull:
      sort_arr_impl<1,5,4,2,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1393410ull:
      sort_arr_impl<1,5,4,3,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1393440ull:
      sort_arr_impl<1,5,4,3,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2102085ull:
      sort_arr_impl<2,0,1,3,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2102100ull:
      sort_arr_impl<2,0,1,3,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2102325ull:
      sort_arr_impl<2,0,1,4,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2102355ull:
      sort_arr_impl<2,0,1,4,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2102580ull:
      sort_arr_impl<2,0,1,5,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2102595ull:
      sort_arr_impl<2,0,1,5,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2109765ull:
      sort_arr_impl<2,0,3,1,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2109780ull:
      sort_arr_impl<2,0,3,1,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2110485ull:
      sort_arr_impl<2,0,3,4,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2110545ull:
      sort_arr_impl<2,0,3,4,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2110740ull:
      sort_arr_impl<2,0,3,5,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2110785ull:
      sort_arr_impl<2,0,3,5,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2113845ull:
      sort_arr_impl<2,0,4,1,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2113875ull:
      sort_arr_impl<2,0,4,1,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2114325ull:
      sort_arr_impl<2,0,4,3,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2114385ull:
      sort_arr_impl<2,0,4,3,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2114835ull:
      sort_arr_impl<2,0,4,5,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2114865ull:
      sort_arr_impl<2,0,4,5,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2117940ull:
      sort_arr_impl<2,0,5,1,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2117955ull:
      sort_arr_impl<2,0,5,1,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2118420ull:
      sort_arr_impl<2,0,5,3,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2118465ull:
      sort_arr_impl<2,0,5,3,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2118675ull:
      sort_arr_impl<2,0,5,4,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2118705ull:
      sort_arr_impl<2,0,5,4,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2163525ull:
      sort_arr_impl<2,1,0,3,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2163540ull:
      sort_arr_impl<2,1,0,3,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2163765ull:
      sort_arr_impl<2,1,0,4,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2163795ull:
      sort_arr_impl<2,1,0,4,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2164020ull:
      sort_arr_impl<2,1,0,5,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2164035ull:
      sort_arr_impl<2,1,0,5,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2175045ull:
      sort_arr_impl<2,1,3,0,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2175060ull:
      sort_arr_impl<2,1,3,0,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2176005ull:
      sort_arr_impl<2,1,3,4,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2176080ull:
      sort_arr_impl<2,1,3,4,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2176260ull:
      sort_arr_impl<2,1,3,5,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2176320ull:
      sort_arr_impl<2,1,3,5,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2179125ull:
      sort_arr_impl<2,1,4,0,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2179155ull:
      sort_arr_impl<2,1,4,0,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2179845ull:
      sort_arr_impl<2,1,4,3,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2179920ull:
      sort_arr_impl<2,1,4,3,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2180355ull:
      sort_arr_impl<2,1,4,5,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2180400ull:
      sort_arr_impl<2,1,4,5,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2183220ull:
      sort_arr_impl<2,1,5,0,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2183235ull:
      sort_arr_impl<2,1,5,0,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2183940ull:
      sort_arr_impl<2,1,5,3,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2184000ull:
      sort_arr_impl<2,1,5,3,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2184195ull:
      sort_arr_impl<2,1,5,4,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2184240ull:
      sort_arr_impl<2,1,5,4,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2294085ull:
      sort_arr_impl<2,3,0,1,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2294100ull:
      sort_arr_impl<2,3,0,1,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2294805ull:
      sort_arr_impl<2,3,0,4,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2294865ull:
      sort_arr_impl<2,3,0,4,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2295060ull:
      sort_arr_impl<2,3,0,5,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2295105ull:
      sort_arr_impl<2,3,0,5,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2297925ull:
      sort_arr_impl<2,3,1,0,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2297940ull:
      sort_arr_impl<2,3,1,0,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2298885ull:
      sort_arr_impl<2,3,1,4,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2298960ull:
      sort_arr_impl<2,3,1,4,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2299140ull:
      sort_arr_impl<2,3,1,5,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2299200ull:
      sort_arr_impl<2,3,1,5,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2310165ull:
      sort_arr_impl<2,3,4,0,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2310225ull:
      sort_arr_impl<2,3,4,0,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2310405ull:
      sort_arr_impl<2,3,4,1,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2310480ull:
      sort_arr_impl<2,3,4,1,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2311425ull:
      sort_arr_impl<2,3,4,5,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2311440ull:
      sort_arr_impl<2,3,4,5,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2314260ull:
      sort_arr_impl<2,3,5,0,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2314305ull:
      sort_arr_impl<2,3,5,0,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2314500ull:
      sort_arr_impl<2,3,5,1,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2314560ull:
      sort_arr_impl<2,3,5,1,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2315265ull:
      sort_arr_impl<2,3,5,4,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2315280ull:
      sort_arr_impl<2,3,5,4,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2359605ull:
      sort_arr_impl<2,4,0,1,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2359635ull:
      sort_arr_impl<2,4,0,1,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2360085ull:
      sort_arr_impl<2,4,0,3,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2360145ull:
      sort_arr_impl<2,4,0,3,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2360595ull:
      sort_arr_impl<2,4,0,5,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2360625ull:
      sort_arr_impl<2,4,0,5,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2363445ull:
      sort_arr_impl<2,4,1,0,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2363475ull:
      sort_arr_impl<2,4,1,0,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2364165ull:
      sort_arr_impl<2,4,1,3,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2364240ull:
      sort_arr_impl<2,4,1,3,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2364675ull:
      sort_arr_impl<2,4,1,5,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2364720ull:
      sort_arr_impl<2,4,1,5,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2371605ull:
      sort_arr_impl<2,4,3,0,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2371665ull:
      sort_arr_impl<2,4,3,0,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2371845ull:
      sort_arr_impl<2,4,3,1,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2371920ull:
      sort_arr_impl<2,4,3,1,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2372865ull:
      sort_arr_impl<2,4,3,5,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2372880ull:
      sort_arr_impl<2,4,3,5,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2379795ull:
      sort_arr_impl<2,4,5,0,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2379825ull:
      sort_arr_impl<2,4,5,0,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2380035ull:
      sort_arr_impl<2,4,5,1,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2380080ull:
      sort_arr_impl<2,4,5,1,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2380545ull:
      sort_arr_impl<2,4,5,3,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2380560ull:
      sort_arr_impl<2,4,5,3,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2425140ull:
      sort_arr_impl<2,5,0,1,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2425155ull:
      sort_arr_impl<2,5,0,1,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2425620ull:
      sort_arr_impl<2,5,0,3,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2425665ull:
      sort_arr_impl<2,5,0,3,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2425875ull:
      sort_arr_impl<2,5,0,4,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2425905ull:
      sort_arr_impl<2,5,0,4,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2428980ull:
      sort_arr_impl<2,5,1,0,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2428995ull:
      sort_arr_impl<2,5,1,0,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2429700ull:
      sort_arr_impl<2,5,1,3,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2429760ull:
      sort_arr_impl<2,5,1,3,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2429955ull:
      sort_arr_impl<2,5,1,4,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2430000ull:
      sort_arr_impl<2,5,1,4,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2437140ull:
      sort_arr_impl<2,5,3,0,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2437185ull:
      sort_arr_impl<2,5,3,0,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2437380ull:
      sort_arr_impl<2,5,3,1,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2437440ull:
      sort_arr_impl<2,5,3,1,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2438145ull:
      sort_arr_impl<2,5,3,4,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2438160ull:
      sort_arr_impl<2,5,3,4,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2441235ull:
      sort_arr_impl<2,5,4,0,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2441265ull:
      sort_arr_impl<2,5,4,0,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2441475ull:
      sort_arr_impl<2,5,4,1,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2441520ull:
      sort_arr_impl<2,5,4,1,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2441985ull:
      sort_arr_impl<2,5,4,3,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 2442000ull:
      sort_arr_impl<2,5,4,3,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3150405ull:
      sort_arr_impl<3,0,1,2,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3150420ull:
      sort_arr_impl<3,0,1,2,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3150885ull:
      sort_arr_impl<3,0,1,4,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3150930ull:
      sort_arr_impl<3,0,1,4,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3151140ull:
      sort_arr_impl<3,0,1,5,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3151170ull:
      sort_arr_impl<3,0,1,5,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3154245ull:
      sort_arr_impl<3,0,2,1,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3154260ull:
      sort_arr_impl<3,0,2,1,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3154965ull:
      sort_arr_impl<3,0,2,4,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3155025ull:
      sort_arr_impl<3,0,2,4,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3155220ull:
      sort_arr_impl<3,0,2,5,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3155265ull:
      sort_arr_impl<3,0,2,5,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3162405ull:
      sort_arr_impl<3,0,4,1,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3162450ull:
      sort_arr_impl<3,0,4,1,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3162645ull:
      sort_arr_impl<3,0,4,2,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3162705ull:
      sort_arr_impl<3,0,4,2,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3163410ull:
      sort_arr_impl<3,0,4,5,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3163425ull:
      sort_arr_impl<3,0,4,5,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3166500ull:
      sort_arr_impl<3,0,5,1,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3166530ull:
      sort_arr_impl<3,0,5,1,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3166740ull:
      sort_arr_impl<3,0,5,2,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3166785ull:
      sort_arr_impl<3,0,5,2,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3167250ull:
      sort_arr_impl<3,0,5,4,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3167265ull:
      sort_arr_impl<3,0,5,4,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3211845ull:
      sort_arr_impl<3,1,0,2,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3211860ull:
      sort_arr_impl<3,1,0,2,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3212325ull:
      sort_arr_impl<3,1,0,4,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3212370ull:
      sort_arr_impl<3,1,0,4,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3212580ull:
      sort_arr_impl<3,1,0,5,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3212610ull:
      sort_arr_impl<3,1,0,5,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3219525ull:
      sort_arr_impl<3,1,2,0,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3219540ull:
      sort_arr_impl<3,1,2,0,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3220485ull:
      sort_arr_impl<3,1,2,4,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3220560ull:
      sort_arr_impl<3,1,2,4,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3220740ull:
      sort_arr_impl<3,1,2,5,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3220800ull:
      sort_arr_impl<3,1,2,5,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3227685ull:
      sort_arr_impl<3,1,4,0,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3227730ull:
      sort_arr_impl<3,1,4,0,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3228165ull:
      sort_arr_impl<3,1,4,2,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3228240ull:
      sort_arr_impl<3,1,4,2,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3228930ull:
      sort_arr_impl<3,1,4,5,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3228960ull:
      sort_arr_impl<3,1,4,5,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3231780ull:
      sort_arr_impl<3,1,5,0,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3231810ull:
      sort_arr_impl<3,1,5,0,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3232260ull:
      sort_arr_impl<3,1,5,2,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3232320ull:
      sort_arr_impl<3,1,5,2,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3232770ull:
      sort_arr_impl<3,1,5,4,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3232800ull:
      sort_arr_impl<3,1,5,4,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3277125ull:
      sort_arr_impl<3,2,0,1,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3277140ull:
      sort_arr_impl<3,2,0,1,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3277845ull:
      sort_arr_impl<3,2,0,4,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3277905ull:
      sort_arr_impl<3,2,0,4,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3278100ull:
      sort_arr_impl<3,2,0,5,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3278145ull:
      sort_arr_impl<3,2,0,5,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3280965ull:
      sort_arr_impl<3,2,1,0,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3280980ull:
      sort_arr_impl<3,2,1,0,5,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3281925ull:
      sort_arr_impl<3,2,1,4,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3282000ull:
      sort_arr_impl<3,2,1,4,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3282180ull:
      sort_arr_impl<3,2,1,5,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3282240ull:
      sort_arr_impl<3,2,1,5,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3293205ull:
      sort_arr_impl<3,2,4,0,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3293265ull:
      sort_arr_impl<3,2,4,0,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3293445ull:
      sort_arr_impl<3,2,4,1,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3293520ull:
      sort_arr_impl<3,2,4,1,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3294465ull:
      sort_arr_impl<3,2,4,5,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3294480ull:
      sort_arr_impl<3,2,4,5,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3297300ull:
      sort_arr_impl<3,2,5,0,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3297345ull:
      sort_arr_impl<3,2,5,0,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3297540ull:
      sort_arr_impl<3,2,5,1,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3297600ull:
      sort_arr_impl<3,2,5,1,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3298305ull:
      sort_arr_impl<3,2,5,4,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3298320ull:
      sort_arr_impl<3,2,5,4,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3408165ull:
      sort_arr_impl<3,4,0,1,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3408210ull:
      sort_arr_impl<3,4,0,1,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3408405ull:
      sort_arr_impl<3,4,0,2,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3408465ull:
      sort_arr_impl<3,4,0,2,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3409170ull:
      sort_arr_impl<3,4,0,5,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3409185ull:
      sort_arr_impl<3,4,0,5,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3412005ull:
      sort_arr_impl<3,4,1,0,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3412050ull:
      sort_arr_impl<3,4,1,0,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3412485ull:
      sort_arr_impl<3,4,1,2,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3412560ull:
      sort_arr_impl<3,4,1,2,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3413250ull:
      sort_arr_impl<3,4,1,5,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3413280ull:
      sort_arr_impl<3,4,1,5,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3416085ull:
      sort_arr_impl<3,4,2,0,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3416145ull:
      sort_arr_impl<3,4,2,0,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3416325ull:
      sort_arr_impl<3,4,2,1,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3416400ull:
      sort_arr_impl<3,4,2,1,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3417345ull:
      sort_arr_impl<3,4,2,5,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3417360ull:
      sort_arr_impl<3,4,2,5,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3428370ull:
      sort_arr_impl<3,4,5,0,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3428385ull:
      sort_arr_impl<3,4,5,0,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3428610ull:
      sort_arr_impl<3,4,5,1,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3428640ull:
      sort_arr_impl<3,4,5,1,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3428865ull:
      sort_arr_impl<3,4,5,2,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3428880ull:
      sort_arr_impl<3,4,5,2,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3473700ull:
      sort_arr_impl<3,5,0,1,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3473730ull:
      sort_arr_impl<3,5,0,1,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3473940ull:
      sort_arr_impl<3,5,0,2,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3473985ull:
      sort_arr_impl<3,5,0,2,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3474450ull:
      sort_arr_impl<3,5,0,4,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3474465ull:
      sort_arr_impl<3,5,0,4,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3477540ull:
      sort_arr_impl<3,5,1,0,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3477570ull:
      sort_arr_impl<3,5,1,0,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3478020ull:
      sort_arr_impl<3,5,1,2,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3478080ull:
      sort_arr_impl<3,5,1,2,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3478530ull:
      sort_arr_impl<3,5,1,4,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3478560ull:
      sort_arr_impl<3,5,1,4,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3481620ull:
      sort_arr_impl<3,5,2,0,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3481665ull:
      sort_arr_impl<3,5,2,0,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3481860ull:
      sort_arr_impl<3,5,2,1,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3481920ull:
      sort_arr_impl<3,5,2,1,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3482625ull:
      sort_arr_impl<3,5,2,4,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3482640ull:
      sort_arr_impl<3,5,2,4,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3489810ull:
      sort_arr_impl<3,5,4,0,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3489825ull:
      sort_arr_impl<3,5,4,0,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3490050ull:
      sort_arr_impl<3,5,4,1,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3490080ull:
      sort_arr_impl<3,5,4,1,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3490305ull:
      sort_arr_impl<3,5,4,2,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 3490320ull:
      sort_arr_impl<3,5,4,2,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4198965ull:
      sort_arr_impl<4,0,1,2,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4198995ull:
      sort_arr_impl<4,0,1,2,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4199205ull:
      sort_arr_impl<4,0,1,3,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4199250ull:
      sort_arr_impl<4,0,1,3,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4199715ull:
      sort_arr_impl<4,0,1,5,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4199730ull:
      sort_arr_impl<4,0,1,5,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4202805ull:
      sort_arr_impl<4,0,2,1,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4202835ull:
      sort_arr_impl<4,0,2,1,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4203285ull:
      sort_arr_impl<4,0,2,3,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4203345ull:
      sort_arr_impl<4,0,2,3,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4203795ull:
      sort_arr_impl<4,0,2,5,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4203825ull:
      sort_arr_impl<4,0,2,5,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4206885ull:
      sort_arr_impl<4,0,3,1,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4206930ull:
      sort_arr_impl<4,0,3,1,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4207125ull:
      sort_arr_impl<4,0,3,2,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4207185ull:
      sort_arr_impl<4,0,3,2,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4207890ull:
      sort_arr_impl<4,0,3,5,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4207905ull:
      sort_arr_impl<4,0,3,5,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4215075ull:
      sort_arr_impl<4,0,5,1,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4215090ull:
      sort_arr_impl<4,0,5,1,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4215315ull:
      sort_arr_impl<4,0,5,2,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4215345ull:
      sort_arr_impl<4,0,5,2,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4215570ull:
      sort_arr_impl<4,0,5,3,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4215585ull:
      sort_arr_impl<4,0,5,3,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4260405ull:
      sort_arr_impl<4,1,0,2,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4260435ull:
      sort_arr_impl<4,1,0,2,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4260645ull:
      sort_arr_impl<4,1,0,3,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4260690ull:
      sort_arr_impl<4,1,0,3,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4261155ull:
      sort_arr_impl<4,1,0,5,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4261170ull:
      sort_arr_impl<4,1,0,5,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4268085ull:
      sort_arr_impl<4,1,2,0,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4268115ull:
      sort_arr_impl<4,1,2,0,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4268805ull:
      sort_arr_impl<4,1,2,3,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4268880ull:
      sort_arr_impl<4,1,2,3,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4269315ull:
      sort_arr_impl<4,1,2,5,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4269360ull:
      sort_arr_impl<4,1,2,5,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4272165ull:
      sort_arr_impl<4,1,3,0,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4272210ull:
      sort_arr_impl<4,1,3,0,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4272645ull:
      sort_arr_impl<4,1,3,2,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4272720ull:
      sort_arr_impl<4,1,3,2,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4273410ull:
      sort_arr_impl<4,1,3,5,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4273440ull:
      sort_arr_impl<4,1,3,5,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4280355ull:
      sort_arr_impl<4,1,5,0,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4280370ull:
      sort_arr_impl<4,1,5,0,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4280835ull:
      sort_arr_impl<4,1,5,2,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4280880ull:
      sort_arr_impl<4,1,5,2,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4281090ull:
      sort_arr_impl<4,1,5,3,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4281120ull:
      sort_arr_impl<4,1,5,3,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4325685ull:
      sort_arr_impl<4,2,0,1,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4325715ull:
      sort_arr_impl<4,2,0,1,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4326165ull:
      sort_arr_impl<4,2,0,3,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4326225ull:
      sort_arr_impl<4,2,0,3,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4326675ull:
      sort_arr_impl<4,2,0,5,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4326705ull:
      sort_arr_impl<4,2,0,5,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4329525ull:
      sort_arr_impl<4,2,1,0,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4329555ull:
      sort_arr_impl<4,2,1,0,5,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4330245ull:
      sort_arr_impl<4,2,1,3,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4330320ull:
      sort_arr_impl<4,2,1,3,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4330755ull:
      sort_arr_impl<4,2,1,5,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4330800ull:
      sort_arr_impl<4,2,1,5,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4337685ull:
      sort_arr_impl<4,2,3,0,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4337745ull:
      sort_arr_impl<4,2,3,0,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4337925ull:
      sort_arr_impl<4,2,3,1,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4338000ull:
      sort_arr_impl<4,2,3,1,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4338945ull:
      sort_arr_impl<4,2,3,5,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4338960ull:
      sort_arr_impl<4,2,3,5,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4345875ull:
      sort_arr_impl<4,2,5,0,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4345905ull:
      sort_arr_impl<4,2,5,0,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4346115ull:
      sort_arr_impl<4,2,5,1,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4346160ull:
      sort_arr_impl<4,2,5,1,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4346625ull:
      sort_arr_impl<4,2,5,3,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4346640ull:
      sort_arr_impl<4,2,5,3,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4391205ull:
      sort_arr_impl<4,3,0,1,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4391250ull:
      sort_arr_impl<4,3,0,1,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4391445ull:
      sort_arr_impl<4,3,0,2,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4391505ull:
      sort_arr_impl<4,3,0,2,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4392210ull:
      sort_arr_impl<4,3,0,5,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4392225ull:
      sort_arr_impl<4,3,0,5,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4395045ull:
      sort_arr_impl<4,3,1,0,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4395090ull:
      sort_arr_impl<4,3,1,0,5,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4395525ull:
      sort_arr_impl<4,3,1,2,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4395600ull:
      sort_arr_impl<4,3,1,2,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4396290ull:
      sort_arr_impl<4,3,1,5,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4396320ull:
      sort_arr_impl<4,3,1,5,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4399125ull:
      sort_arr_impl<4,3,2,0,1,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4399185ull:
      sort_arr_impl<4,3,2,0,5,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4399365ull:
      sort_arr_impl<4,3,2,1,0,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4399440ull:
      sort_arr_impl<4,3,2,1,5,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4400385ull:
      sort_arr_impl<4,3,2,5,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4400400ull:
      sort_arr_impl<4,3,2,5,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4411410ull:
      sort_arr_impl<4,3,5,0,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4411425ull:
      sort_arr_impl<4,3,5,0,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4411650ull:
      sort_arr_impl<4,3,5,1,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4411680ull:
      sort_arr_impl<4,3,5,1,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4411905ull:
      sort_arr_impl<4,3,5,2,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4411920ull:
      sort_arr_impl<4,3,5,2,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4522275ull:
      sort_arr_impl<4,5,0,1,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4522290ull:
      sort_arr_impl<4,5,0,1,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4522515ull:
      sort_arr_impl<4,5,0,2,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4522545ull:
      sort_arr_impl<4,5,0,2,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4522770ull:
      sort_arr_impl<4,5,0,3,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4522785ull:
      sort_arr_impl<4,5,0,3,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4526115ull:
      sort_arr_impl<4,5,1,0,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4526130ull:
      sort_arr_impl<4,5,1,0,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4526595ull:
      sort_arr_impl<4,5,1,2,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4526640ull:
      sort_arr_impl<4,5,1,2,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4526850ull:
      sort_arr_impl<4,5,1,3,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4526880ull:
      sort_arr_impl<4,5,1,3,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4530195ull:
      sort_arr_impl<4,5,2,0,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4530225ull:
      sort_arr_impl<4,5,2,0,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4530435ull:
      sort_arr_impl<4,5,2,1,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4530480ull:
      sort_arr_impl<4,5,2,1,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4530945ull:
      sort_arr_impl<4,5,2,3,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4530960ull:
      sort_arr_impl<4,5,2,3,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4534290ull:
      sort_arr_impl<4,5,3,0,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4534305ull:
      sort_arr_impl<4,5,3,0,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4534530ull:
      sort_arr_impl<4,5,3,1,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4534560ull:
      sort_arr_impl<4,5,3,1,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4534785ull:
      sort_arr_impl<4,5,3,2,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 4534800ull:
      sort_arr_impl<4,5,3,2,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5247540ull:
      sort_arr_impl<5,0,1,2,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5247555ull:
      sort_arr_impl<5,0,1,2,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5247780ull:
      sort_arr_impl<5,0,1,3,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5247810ull:
      sort_arr_impl<5,0,1,3,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5248035ull:
      sort_arr_impl<5,0,1,4,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5248050ull:
      sort_arr_impl<5,0,1,4,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5251380ull:
      sort_arr_impl<5,0,2,1,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5251395ull:
      sort_arr_impl<5,0,2,1,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5251860ull:
      sort_arr_impl<5,0,2,3,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5251905ull:
      sort_arr_impl<5,0,2,3,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5252115ull:
      sort_arr_impl<5,0,2,4,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5252145ull:
      sort_arr_impl<5,0,2,4,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5255460ull:
      sort_arr_impl<5,0,3,1,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5255490ull:
      sort_arr_impl<5,0,3,1,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5255700ull:
      sort_arr_impl<5,0,3,2,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5255745ull:
      sort_arr_impl<5,0,3,2,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5256210ull:
      sort_arr_impl<5,0,3,4,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5256225ull:
      sort_arr_impl<5,0,3,4,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5259555ull:
      sort_arr_impl<5,0,4,1,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5259570ull:
      sort_arr_impl<5,0,4,1,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5259795ull:
      sort_arr_impl<5,0,4,2,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5259825ull:
      sort_arr_impl<5,0,4,2,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5260050ull:
      sort_arr_impl<5,0,4,3,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5260065ull:
      sort_arr_impl<5,0,4,3,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5308980ull:
      sort_arr_impl<5,1,0,2,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5308995ull:
      sort_arr_impl<5,1,0,2,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5309220ull:
      sort_arr_impl<5,1,0,3,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5309250ull:
      sort_arr_impl<5,1,0,3,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5309475ull:
      sort_arr_impl<5,1,0,4,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5309490ull:
      sort_arr_impl<5,1,0,4,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5316660ull:
      sort_arr_impl<5,1,2,0,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5316675ull:
      sort_arr_impl<5,1,2,0,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5317380ull:
      sort_arr_impl<5,1,2,3,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5317440ull:
      sort_arr_impl<5,1,2,3,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5317635ull:
      sort_arr_impl<5,1,2,4,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5317680ull:
      sort_arr_impl<5,1,2,4,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5320740ull:
      sort_arr_impl<5,1,3,0,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5320770ull:
      sort_arr_impl<5,1,3,0,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5321220ull:
      sort_arr_impl<5,1,3,2,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5321280ull:
      sort_arr_impl<5,1,3,2,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5321730ull:
      sort_arr_impl<5,1,3,4,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5321760ull:
      sort_arr_impl<5,1,3,4,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5324835ull:
      sort_arr_impl<5,1,4,0,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5324850ull:
      sort_arr_impl<5,1,4,0,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5325315ull:
      sort_arr_impl<5,1,4,2,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5325360ull:
      sort_arr_impl<5,1,4,2,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5325570ull:
      sort_arr_impl<5,1,4,3,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5325600ull:
      sort_arr_impl<5,1,4,3,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5374260ull:
      sort_arr_impl<5,2,0,1,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5374275ull:
      sort_arr_impl<5,2,0,1,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5374740ull:
      sort_arr_impl<5,2,0,3,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5374785ull:
      sort_arr_impl<5,2,0,3,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5374995ull:
      sort_arr_impl<5,2,0,4,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5375025ull:
      sort_arr_impl<5,2,0,4,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5378100ull:
      sort_arr_impl<5,2,1,0,3,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5378115ull:
      sort_arr_impl<5,2,1,0,4,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5378820ull:
      sort_arr_impl<5,2,1,3,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5378880ull:
      sort_arr_impl<5,2,1,3,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5379075ull:
      sort_arr_impl<5,2,1,4,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5379120ull:
      sort_arr_impl<5,2,1,4,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5386260ull:
      sort_arr_impl<5,2,3,0,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5386305ull:
      sort_arr_impl<5,2,3,0,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5386500ull:
      sort_arr_impl<5,2,3,1,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5386560ull:
      sort_arr_impl<5,2,3,1,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5387265ull:
      sort_arr_impl<5,2,3,4,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5387280ull:
      sort_arr_impl<5,2,3,4,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5390355ull:
      sort_arr_impl<5,2,4,0,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5390385ull:
      sort_arr_impl<5,2,4,0,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5390595ull:
      sort_arr_impl<5,2,4,1,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5390640ull:
      sort_arr_impl<5,2,4,1,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5391105ull:
      sort_arr_impl<5,2,4,3,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5391120ull:
      sort_arr_impl<5,2,4,3,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5439780ull:
      sort_arr_impl<5,3,0,1,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5439810ull:
      sort_arr_impl<5,3,0,1,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5440020ull:
      sort_arr_impl<5,3,0,2,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5440065ull:
      sort_arr_impl<5,3,0,2,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5440530ull:
      sort_arr_impl<5,3,0,4,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5440545ull:
      sort_arr_impl<5,3,0,4,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5443620ull:
      sort_arr_impl<5,3,1,0,2,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5443650ull:
      sort_arr_impl<5,3,1,0,4,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5444100ull:
      sort_arr_impl<5,3,1,2,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5444160ull:
      sort_arr_impl<5,3,1,2,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5444610ull:
      sort_arr_impl<5,3,1,4,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5444640ull:
      sort_arr_impl<5,3,1,4,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5447700ull:
      sort_arr_impl<5,3,2,0,1,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5447745ull:
      sort_arr_impl<5,3,2,0,4,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5447940ull:
      sort_arr_impl<5,3,2,1,0,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5448000ull:
      sort_arr_impl<5,3,2,1,4,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5448705ull:
      sort_arr_impl<5,3,2,4,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5448720ull:
      sort_arr_impl<5,3,2,4,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5455890ull:
      sort_arr_impl<5,3,4,0,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5455905ull:
      sort_arr_impl<5,3,4,0,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5456130ull:
      sort_arr_impl<5,3,4,1,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5456160ull:
      sort_arr_impl<5,3,4,1,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5456385ull:
      sort_arr_impl<5,3,4,2,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5456400ull:
      sort_arr_impl<5,3,4,2,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5505315ull:
      sort_arr_impl<5,4,0,1,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5505330ull:
      sort_arr_impl<5,4,0,1,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5505555ull:
      sort_arr_impl<5,4,0,2,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5505585ull:
      sort_arr_impl<5,4,0,2,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5505810ull:
      sort_arr_impl<5,4,0,3,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5505825ull:
      sort_arr_impl<5,4,0,3,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5509155ull:
      sort_arr_impl<5,4,1,0,2,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5509170ull:
      sort_arr_impl<5,4,1,0,3,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5509635ull:
      sort_arr_impl<5,4,1,2,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5509680ull:
      sort_arr_impl<5,4,1,2,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5509890ull:
      sort_arr_impl<5,4,1,3,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5509920ull:
      sort_arr_impl<5,4,1,3,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5513235ull:
      sort_arr_impl<5,4,2,0,1,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5513265ull:
      sort_arr_impl<5,4,2,0,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5513475ull:
      sort_arr_impl<5,4,2,1,0,3>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5513520ull:
      sort_arr_impl<5,4,2,1,3,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5513985ull:
      sort_arr_impl<5,4,2,3,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5514000ull:
      sort_arr_impl<5,4,2,3,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5517330ull:
      sort_arr_impl<5,4,3,0,1,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5517345ull:
      sort_arr_impl<5,4,3,0,2,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5517570ull:
      sort_arr_impl<5,4,3,1,0,2>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5517600ull:
      sort_arr_impl<5,4,3,1,2,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5517825ull:
      sort_arr_impl<5,4,3,2,0,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 5517840ull:
      sort_arr_impl<5,4,3,2,1,0>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5]); break;
    case 1193046ull:
      sort_arr_impl<0,1,2,3,4,5,6>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6]); break;
    case 1197621ull:
      sort_arr_impl<0,1,2,4,6,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6]); break;
    case 2311446ull:
      sort_arr_impl<0,2,3,4,5,1,6>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6]); break;
    case 2385201ull:
      sort_arr_impl<0,2,4,6,5,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6]); break;
    case 20205606ull:
      sort_arr_impl<1,3,4,5,0,2,6>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6]); break;
    case 38142261ull:
      sort_arr_impl<2,4,6,0,1,3,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6]); break;
    case 19088743ull:
      sort_arr_impl<0,1,2,3,4,5,6,7>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 19097413ull:
      sort_arr_impl<0,1,2,3,6,7,4,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 20211493ull:
      sort_arr_impl<0,1,3,4,6,7,2,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 23213143ull:
      sort_arr_impl<0,1,6,2,3,4,5,7>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 36983143ull:
      sort_arr_impl<0,2,3,4,5,1,6,7>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 72573223ull:
      sort_arr_impl<0,4,5,3,6,1,2,7>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 587285863ull:
      sort_arr_impl<2,3,0,1,4,5,6,7>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 623992948ull:
      sort_arr_impl<2,5,3,1,6,0,7,4>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 874845703ull:
      sort_arr_impl<3,4,2,5,1,6,0,7>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 1442915143ull:
      sort_arr_impl<5,6,0,1,2,3,4,7>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 1679823703ull:
      sort_arr_impl<6,4,2,0,1,3,5,7>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 1966146118ull:
      sort_arr_impl<7,5,3,1,0,2,4,6>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7]); break;
    case 305419896ull:
      sort_arr_impl<0,1,2,3,4,5,6,7,8>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8]); break;
    case 610825521ull:
      sort_arr_impl<0,2,4,6,8,7,5,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8]); break;
    case 36039562071ull:
      sort_arr_impl<8,6,4,2,0,1,3,5,7>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8]); break;
    case 4886718345ull:
      sort_arr_impl<0,1,2,3,4,5,6,7,8,9>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8], dim[9]); break;
    case 4886988645ull:
      sort_arr_impl<0,1,2,3,4,9,8,7,6,5>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8], dim[9]); break;
    case 9772807545ull:
      sort_arr_impl<0,2,4,6,8,1,3,5,7,9>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8], dim[9]); break;
    case 9773348145ull:
      sort_arr_impl<0,2,4,6,8,9,7,5,3,1>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8], dim[9]); break;
    case 576632993145ull:
      sort_arr_impl<8,6,4,2,0,1,3,5,7,9>(inp, out, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7], dim[8], dim[9]); break;
    default : {
      std::vector<size_t> tmp(nr);

      std::cout << "Warning: tensor permutation is not optimized " << tag << std::endl;
      std::cout << " " << nr;
      for (size_t k = 0; k != nr; ++k)
        std::cout << " " << pmt[k];
      std::cout << std::endl;

      size_t ns = 1u;
      for_each(dim.begin(), dim.end(), [&](size_t i) {ns *= i;});
      for (size_t i = 0; i != ns; ++i) {
        size_t dist = i;
        for (size_t j = 0; j != nr; ++j) {
          tmp[j] = dist%dim[j];
          dist = dist/dim[j];
        }
        dist = tmp[pmt[nr-1]];
        for (size_t j = nr-1; j-->0;) {
          dist *= dim[pmt[j]];
          dist += tmp[pmt[j]];
        }
        out[dist] = inp[i];
      }
    }
  }
}
#endif /* _SORT_H_  */
