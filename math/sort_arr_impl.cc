#include "sort_arr_impl.h"
#include "mkl_extern.h"


inline void transpose_arr(const double alpha, const int nrows, const int ncols, const double* inp, double* out) {
  mkl_domatcopy_("C","T", &nrows, &ncols, &alpha, inp, &nrows, out, &ncols);
  return;
}

template <>
void sort_arr_impl<0,1>(const double* inp, double* out, const size_t d0, const size_t d1) {
  std::copy_n(inp, d0*d1, out);
  return;
}

template <>
void sort_arr_impl<1,0>(const double* inp, double* out, const size_t d0, const size_t d1) {
  transpose_arr(1.0, d0, d1, inp, out);
  return;
}

template<>
void sort_arr_impl<0,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2) {
  std::copy_n(inp, d0*d1*d2, out);
  return;
}

template<>
void sort_arr_impl<0,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2) {
  for (size_t i = 0; i != d2; ++i)
    for (size_t j = 0; j != d1; ++j)
      std::copy_n(inp + d0*(j + d1 * i), d0, out + d0 * (i + d2 * j));
  return;
}

template<>
void sort_arr_impl<1,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2) {
  for (size_t i = 0, off = 0; i != d2; ++i, off += d0*d1)
    sort_arr_impl<1,0>(inp+off, out+off, d0, d1);
  return;
}

template<>
void sort_arr_impl<1,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2) {
  transpose_arr(1.0, d0, d1*d2, inp, out);
  return;
}

template<>
void sort_arr_impl<2,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2) {
  transpose_arr(1.0, d0*d1, d2, inp, out);
  return;
}

template<>
void sort_arr_impl<0,1,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  std::copy_n(inp, d0*d1*d2*d3, out); 
  return;
}

template<>
void sort_arr_impl<0,1,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  sort_arr_impl<0,2,1>(inp, out, d0*d1, d2, d3);
  return;
}

template<>
void sort_arr_impl<0,2,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  for (size_t j3 = 0; j3 != d3; ++j3)
    for (size_t j2 = 0; j2 != d2; ++j2)
      for (size_t j1 = 0; j1 != d1; ++j1)
        std::copy_n(inp + d0*(j1+d1*(j2+d2*j3)), d0, out + d0*(j2+d2*(j1+d1*j3)));
  return;
}

template<>
void sort_arr_impl<0,2,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  sort_arr_impl<0,2,1>(inp, out, d0, d1, d2*d3);
  return;
}

template<>
void sort_arr_impl<0,3,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  sort_arr_impl<0,2,1>(inp, out, d0, d1*d2, d3);
  return;
}

template<>
void sort_arr_impl<0,3,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  for (size_t j3 = 0; j3 != d3; ++j3)
    for (size_t j2 = 0; j2 != d2; ++j2)
      for (size_t j1 = 0; j1 != d1; ++j1)
        std::copy_n(inp + d0*(j1+d1*(j2+d2*j3)), d0, out + d0*(j3+d3*(j2+d2*j1)));
  return;
}

template<>
void sort_arr_impl<1,0,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  sort_arr_impl<1,0,2>(inp, out, d0, d1, d2*d3);
  return;
}

template<>
void sort_arr_impl<1,2,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  sort_arr_impl<1,0,2>(inp, out, d0, d1*d2, d3);
  return;
}

template<>
void sort_arr_impl<1,2,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  const double alpha = 1.0;
  transpose_arr(alpha, d0, d1*d2*d3, inp, out);
  return;
}

template<>
void sort_arr_impl<2,0,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  sort_arr_impl<1,0,2>(inp, out, d0*d1, d2, d3);
  return;
}

template<>
void sort_arr_impl<2,3,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  transpose_arr(1.0, d0*d1, d2*d3, inp, out);
  return;
}

template<>
void sort_arr_impl<3,0,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3) {
  transpose_arr(1.0, d0*d1*d2, d3, inp, out);
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  std::copy_n(inp, d0*d1*d2*d3*d4, out); 
  return;
}

template<>
void sort_arr_impl<0,1,2,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,2,1>(inp, out, d0*d1*d2, d3, d4);
  return;
}

template<>
void sort_arr_impl<0,1,3,2,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,2,1,3>(inp, out, d0*d1, d2, d3, d4);    
  return;
}

template<>
void sort_arr_impl<0,1,3,4,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,2,3,1>(inp, out, d0*d1, d2, d3, d4); 
  return;
}

template<>
void sort_arr_impl<0,1,4,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,3,1,2>(inp, out, d0*d1, d2, d3, d4);
  return;
}

template<>
void sort_arr_impl<0,1,4,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,3,2,1>(inp, out, d0*d1, d2, d3, d4); 
  return;
}

template<>
void sort_arr_impl<0,2,1,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,2,1,3>(inp, out, d0, d1, d2, d3*d4);
  return;
}

template<>
void sort_arr_impl<0,2,1,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  for (size_t j4 = 0; j4 != d4; ++j4) {
    for (size_t j3 = 0; j3 != d3; ++j3) {
      for (size_t j2 = 0; j2 != d2; ++j2) {
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*j4))), d0, out + d0*(j2+d2*(j1+d1*(j4+d4*j3))));
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,2,3,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,2,1,3>(inp, out, d0, d1, d2*d3, d4);  
  return;
}

template<>
void sort_arr_impl<0,2,3,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,2,1>(inp, out, d0, d1, d2*d3*d4);
  return;
}

template<>
void sort_arr_impl<0,2,4,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  for (size_t j4 = 0; j4 != d4; ++j4) {
    for (size_t j3 = 0; j3 != d3; ++j3) {
      for (size_t j2 = 0; j2 != d2; ++j2) {
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*j4))), d0, out + d0*(j2+d2*(j4+d4*(j1+d1*j3))));
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,2,4,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  for (size_t j4 = 0; j4 != d4; ++j4) {
    for (size_t j3 = 0; j3 != d3; ++j3) {
      for (size_t j2 = 0; j2 != d2; ++j2) {
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*j4))), d0, out + d0*(j2+d2*(j4+d4*(j3+d3*j1))));
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,3,1,2,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,2,1,3>(inp, out, d0, d1*d2, d3, d4);
  return;
}

template<>
void sort_arr_impl<0,3,1,4,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  for (size_t j4 = 0; j4 != d4; ++j4) {
    for (size_t j3 = 0; j3 != d3; ++j3) {
      for (size_t j2 = 0; j2 != d2; ++j2) {
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*j4))), d0, out + d0*(j3+d3*(j1+d1*(j4+d4*j2))));
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,3,2,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  for (size_t j4 = 0; j4 != d4; ++j4) {
    for (size_t j3 = 0; j3 != d3; ++j3) {
      for (size_t j2 = 0; j2 != d2; ++j2) {
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*j4))), d0, out + d0*(j3+d3*(j2+d2*(j1+d1*j4))));
        }
      }
    }
  } 
  return;
}

template<>
void sort_arr_impl<0,3,2,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  for (size_t j4 = 0; j4 != d4; ++j4) {
    for (size_t j3 = 0; j3 != d3; ++j3) {
      for (size_t j2 = 0; j2 != d2; ++j2) {
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*j4))), d0, out + d0*(j3+d3*(j2+d2*(j4+d4*j1))));
        }
      }
    }
  }   
  return;
}

template<>
void sort_arr_impl<0,3,4,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,2,3,1>(inp, out, d0, d1*d2, d3, d4);
  return;
}

template<>
void sort_arr_impl<0,3,4,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  for (size_t j4 = 0; j4 != d4; ++j4) {
    for (size_t j3 = 0; j3 != d3; ++j3) {
      for (size_t j2 = 0; j2 != d2; ++j2) {
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*j4))), d0, out + d0*(j3+d3*(j4+d4*(j2+d2*j1))));
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,4,1,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,2,1>(inp, out, d0, d1*d2*d3, d4);
  return;
}

template<>
void sort_arr_impl<0,4,1,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  for (size_t j4 = 0; j4 != d4; ++j4) {
    for (size_t j3 = 0; j3 != d3; ++j3) {
      for (size_t j2 = 0; j2 != d2; ++j2) {
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*j4))), d0, out + d0*(j4+d4*(j1+d1*(j3+d3*j2))));
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,4,2,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  for (size_t j4 = 0; j4 != d4; ++j4) {
    for (size_t j3 = 0; j3 != d3; ++j3) {
      for (size_t j2 = 0; j2 != d2; ++j2) {
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*j4))), d0, out + d0*(j4+d4*(j2+d2*(j1+d1*j3))));
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<0,4,2,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,3,2,1>(inp, out, d0, d1, d2*d3, d4);
  return;
}

template<>
void sort_arr_impl<0,4,3,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<0,3,2,1>(inp, out, d0, d1*d2, d3, d4);
  return;
}

template<>
void sort_arr_impl<0,4,3,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  for (size_t j4 = 0; j4 != d4; ++j4) {
    for (size_t j3 = 0; j3 != d3; ++j3) {
      for (size_t j2 = 0; j2 != d2; ++j2) {
        for (size_t j1 = 0; j1 != d1; ++ j1) {
          std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*j4))), d0, out + d0*(j4+d4*(j3+d3*(j2+d2*j1))));
        }
      }
    }
  }
  return;
}

template<>
void sort_arr_impl<1,0,2,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0,2>(inp, out, d0, d1, d2*d3*d4);
  return;
}

template<>
void sort_arr_impl<1,0,2,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1>(inp, tmp.get(), d0*d1*d2, d3, d4);
  sort_arr_impl<1,0,2>(tmp.get(), out, d0, d1, d2*d4*d3);
  return;
}

//template<>
//void sort_arr_impl<1,0,3,2,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//  return;
//}

template<>
void sort_arr_impl<1,0,3,4,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,3,1>(inp, tmp.get(), d0*d1, d2, d3, d4);
  sort_arr_impl<1,0,2>(tmp.get(), out, d0, d1, d2*d3*d4);
  return;
}

template<>
void sort_arr_impl<1,0,4,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,3,1,2>(inp, tmp.get(), d0*d1, d2, d3, d4);
  sort_arr_impl<1,0,2>(tmp.get(), out, d0, d1, d2*d3*d4);
  return;
}

//template<>
//void sort_arr_impl<1,0,4,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//  return;
//}

template<>
void sort_arr_impl<1,2,0,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0,2>(inp, out, d0, d1*d2, d3*d4);
  return;
}

template<>
void sort_arr_impl<1,2,0,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1>(inp, tmp.get(), d0*d1*d2, d3, d4);
  sort_arr_impl<1,0,2>(tmp.get(), out, d0, d1*d2, d3*d4);
  return;
}

template<>
void sort_arr_impl<1,2,3,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0,2>(inp, out, d0, d1*d2*d3, d4);
  return;
}

template<>
void sort_arr_impl<1,2,3,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0>(inp, out, d0, d1*d2*d3*d4);
  return;
}

template<>
void sort_arr_impl<1,2,4,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1>(inp, tmp.get(), d0*d1*d2, d3, d4);
  sort_arr_impl<1,0,2>(tmp.get(), out, d0, d1*d2*d4, d3);
  return;
}

template<>
void sort_arr_impl<1,2,4,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1>(inp, tmp.get(), d0*d1*d2, d3, d4);
  sort_arr_impl<1,0>(tmp.get(), out, d0, d1*d2*d3*d4);
  return;
}

//template<>
//void sort_arr_impl<1,3,2,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//  return;
//}

template<>
void sort_arr_impl<1,3,2,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1,3>(inp, tmp.get(), d0*d1, d2, d3, d4);
  sort_arr_impl<1,0>(tmp.get(), out, d0, d1*d2*d3*d4);
  return;
}

template<>
void sort_arr_impl<1,3,4,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1>(inp, tmp.get(), d0*d1, d2, d3*d4);
  sort_arr_impl<1,0,2>(tmp.get(), out, d0, d1*d3*d4, d2);
  return;
}

template<>
void sort_arr_impl<1,3,4,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1>(inp, tmp.get(), d0*d1, d2, d3*d4);
  sort_arr_impl<1,0>(tmp.get(), out, d0, d1*d2*d3*d4);
  return;
}

template<>
void sort_arr_impl<1,4,0,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1>(inp, tmp.get(), d0*d1, d2*d3, d4);
  sort_arr_impl<1,0,2>(tmp.get(), out, d0, d1*d4, d2*d3);
  return;
}

//template<>
//void sort_arr_impl<1,4,0,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//  return;
//}

template<>
void sort_arr_impl<1,4,2,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,3,1,2>(inp, tmp.get(), d0*d1, d2, d3, d4);
  sort_arr_impl<1,0,2>(tmp.get(), out, d0, d1*d2*d4, d3);
  return;
}

template<>
void sort_arr_impl<1,4,2,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,3,1,2>(inp, tmp.get(), d0*d1, d2, d3, d4);
  sort_arr_impl<1,0>(tmp.get(), out, d0, d1*d2*d3*d4);
  return;
}

//template<>
//void sort_arr_impl<1,4,3,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//  return;
//}

template<>
void sort_arr_impl<1,4,3,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,3,2,1>(inp, tmp.get(), d0*d1, d2, d3, d4);
  sort_arr_impl<1,0>(tmp.get(), out, d0, d1*d2*d3*d4);
  return;
}

template<>
void sort_arr_impl<2,0,1,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0,2>(inp, out, d0*d1, d2, d3*d4);
  return;
}

template<>
void sort_arr_impl<2,0,1,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0,3,2>(inp, out, d0*d1, d2, d3, d4);
  return;
}

//template<>
//void sort_arr_impl<2,0,3,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {  
//  return;
//}

template<>
void sort_arr_impl<2,0,3,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<2,0,3,1>(inp, out, d0, d1, d2, d3*d4);
  return;
}

//template<>
//void sort_arr_impl<2,0,4,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<2,0,4,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//  return;
//}

template<>
void sort_arr_impl<2,1,0,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<2,1,0,3>(inp, out, d0, d1, d2, d3*d4);
  return;
}

//template<>
//void sort_arr_impl<2,1,0,4,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//  return;
//}

//template<>
//void sort_arr_impl<2,1,3,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<2,1,3,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<2,1,4,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<2,1,4,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<2,3,0,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0,2>(inp, out, d0*d1, d2*d3, d4);
  return;
}

template<>
void sort_arr_impl<2,3,0,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<1,0>(inp, tmp.get(), d0*d1, d2*d3*d4);
  sort_arr_impl<0,2,1,3>(tmp.get(), out, d2*d3, d4, d0, d1);
  return;
}

//template<>
//void sort_arr_impl<2,3,1,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<2,3,1,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<1,0>(inp, tmp.get(), d0, d1*d2*d3*d4);
  sort_arr_impl<1,0,2>(tmp.get(), out, d1, d2*d3, d0*d4);
  return;
}

template<>
void sort_arr_impl<2,3,4,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0>(inp, out, d0*d1, d2*d3*d4);
  return;
}

template<>
void sort_arr_impl<2,3,4,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<1,0>(inp, tmp.get(), d0, d1*d2*d3*d4);
  sort_arr_impl<1,0,2>(tmp.get(), out, d1, d2*d3*d4, d0);
  return;
}

template<>
void sort_arr_impl<2,4,0,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,3,0,2>(inp, out, d0*d1, d2, d3, d4);
  return;
}

//template<>
//void sort_arr_impl<2,4,0,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<2,4,1,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<2,4,1,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<2,4,3,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1>(inp, tmp.get(), d0*d1*d2, d3,d4);
  sort_arr_impl<1,0>(tmp.get(), out, d0*d1, d2*d3*d4);
  return;
}

//template<>
//void sort_arr_impl<2,4,3,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<3,0,1,2,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0,2>(inp, out, d0*d1*d2, d3, d4);
  return;
}

template<>
void sort_arr_impl<3,0,1,4,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<2,0,3,1>(inp, out, d0*d1, d2, d3, d4);
  return;
}

//template<>
//void sort_arr_impl<3,0,2,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}


//template<>
//void sort_arr_impl<3,0,2,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<3,0,4,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<1,0>(inp, tmp.get(), d0*d1*d2, d3*d4);
  sort_arr_impl<0,2,1,3>(tmp.get(), out, d3, d4, d0, d1*d2);
  return;
}

//template<>
//void sort_arr_impl<3,0,4,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<3,1,0,2,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<3,1,0,4,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<3,1,2,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<2,1,0,3>(inp, out, d0, d1*d2, d3, d4);
  return;
}

template<>
void sort_arr_impl<3,1,2,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<2,1,3,0>(inp, out, d0, d1*d2, d3, d4);
  return;
}

//template<>
//void sort_arr_impl<3,1,4,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<3,1,4,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<3,2,0,1,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<2,1,0,3>(inp, out, d0*d1, d2, d3, d4);
  return;
}

//template<>
//void sort_arr_impl<3,2,0,4,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<3,2,1,0,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<3,2,1,4,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<3,2,4,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<1,0>(inp, tmp.get(), d0*d1*d2, d3*d4);
  sort_arr_impl<0,2,1>(tmp.get(), out, d3, d0*d1*d4, d2);
  return;
}

//template<>
//void sort_arr_impl<3,2,4,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<3,4,0,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0>(inp, out, d0*d1*d2, d3*d4);
  return;
}

template<>
void sort_arr_impl<3,4,0,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<1,0>(inp, tmp.get(), d0*d1*d2, d3*d4);
  sort_arr_impl<0,2,1>(tmp.get(), out, d0*d3*d4, d1, d2);
  return;
}

template<>
void sort_arr_impl<3,4,1,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<3,1,0,2>(inp, out, d0, d1, d2, d3*d4);
  return;
}

template<>
void sort_arr_impl<3,4,1,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<2,1,0>(inp, out, d0, d1*d2, d3*d4);
  return;
}

template<>
void sort_arr_impl<3,4,2,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<2,1,0>(inp, out, d0*d1, d2, d3*d4);
  return;
}

//template<>
//void sort_arr_impl<3,4,2,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<4,0,1,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<1,0>(inp, out, d0*d1*d2*d3, d4);
  return;
}

template<>
void sort_arr_impl<4,0,1,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<1,0>(inp, tmp.get(), d0*d1*d2*d3, d4);
  sort_arr_impl<0,2,1>(tmp.get(), out, d0*d1*d4, d2, d3);
  return;
}

template<>
void sort_arr_impl<4,0,2,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  const size_t ns = d0*d1*d2*d3*d4;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1,3>(inp, tmp.get(), d0, d1, d2, d3*d4);
  sort_arr_impl<1,0>(tmp.get(), out, d0*d1*d2*d3, d4);
  return;
}

//template<>
//void sort_arr_impl<4,0,2,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<4,0,3,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<3,0,2,1>(inp, out, d0, d1*d2, d3, d4);
  return;
}

//template<>
//void sort_arr_impl<4,0,3,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<4,1,0,2,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<3,1,0,2>(inp, out, d0, d1, d2*d3, d4);
  return;
}

//template<>
//void sort_arr_impl<4,1,0,3,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<4,1,2,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<3,1,0,2>(inp, out, d0, d1*d2, d3, d4);
  return;
}

template<>
void sort_arr_impl<4,1,2,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<2,1,0>(inp, out, d0, d1*d2*d3, d4);
  return;
}

//template<>
//void sort_arr_impl<4,1,3,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<4,1,3,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<4,2,0,1,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<3,1,0,2>(inp, out, d0*d1, d2, d3, d4);
  return;
}

//template<>
//void sort_arr_impl<4,2,0,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//  return;
//}

//template<>
//void sort_arr_impl<4,2,1,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

//template<>
//void sort_arr_impl<4,2,1,3,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<4,2,3,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<3,1,2,0>(inp, out, d0*d1, d2, d3, d4);
  return;
}

template<>
void sort_arr_impl<4,2,3,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<3,2,1,0>(inp, out, d0, d1, d2*d3, d4);
  return;
}

template<>
void sort_arr_impl<4,3,0,1,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<2,1,0>(inp, out, d0*d1*d2, d3, d4);
  return;
}

//template<>
//void sort_arr_impl<4,3,0,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//  return;
//}

//template<>
//void sort_arr_impl<4,3,1,0,2>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//  return;
//}

template<>
void sort_arr_impl<4,3,1,2,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<3,2,1,0>(inp, out, d0, d1*d2, d3, d4);
  return;
}

template<>
void sort_arr_impl<4,3,2,0,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
  sort_arr_impl<3,2,1,0>(inp, out, d0*d1, d2, d3, d4);
  return;
}

//template<>
//void sort_arr_impl<4,3,2,1,0>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4) {
//   return;
//}

template<>
void sort_arr_impl<0,1,2,3,4,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5) {
  std::copy_n(inp, d0*d1*d2*d3*d4*d5, out); 
  return;
}

template<>
void sort_arr_impl<0,2,4,1,3,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5) {
  for (size_t j5 = 0; j5 != d5; ++j5)
    for (size_t j4 = 0; j4 != d4; ++j4)
      for (size_t j3 = 0; j3 != d3; ++j3)
        for (size_t j2 = 0; j2 != d2; ++j2)
          for (size_t j1 = 0; j1 != d1; ++j1)
            std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*(j4+d4*j5)))), d0, out + d0*(j2+d2*(j4+d4*(j1+d1*(j3+d3*j5)))));
  return;
}

template<>
void sort_arr_impl<0,3,1,2,4,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5) {
  sort_arr_impl<0,2,1,3>(inp, out, d0, d1*d2, d3, d4*d5);
  return;
}

template<>
void sort_arr_impl<0,5,4,3,2,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5) {
  for (size_t j5 = 0; j5 != d5; ++j5)
    for (size_t j4 = 0; j4 != d4; ++j4)
      for (size_t j3 = 0; j3 != d3; ++j3)
        for (size_t j2 = 0; j2 != d2; ++j2)
          for (size_t j1 = 0; j1 != d1; ++j1)
            std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*(j4+d4*j5)))), d0, out + d0*(j5+d5*(j4+d4*(j3+d3*(j2+d2*j1)))));
  return;
}

template<>
void sort_arr_impl<1,2,4,5,0,3>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5) {
  const size_t ns = d0*d1*d2*d3*d4*d5;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<1,0>(inp, tmp.get(), d0, d1*d2*d3*d4*d5);
  sort_arr_impl<0,2,1>(tmp.get(), out, d1*d2, d3, d0*d4*d5);
  return;
}

template<>
void sort_arr_impl<1,2,5,0,3,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5) {
  const size_t ns = d0*d1*d2*d3*d4*d5;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<1,0>(inp, tmp.get(), d0, d1*d2*d3*d4*d5);
  sort_arr_impl<0,2,1>(tmp.get(), out, d1*d2, d3*d4, d0*d5);
  return;
}

template<>
void sort_arr_impl<2,3,0,1,4,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5) {
  sort_arr_impl<1,0,2>(inp, out, d0*d1, d2*d3, d4*d5);
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4,5,6>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6) {
  std::copy_n(inp, d0*d1*d2*d3*d4*d5*d6, out); 
  return;
}

template<>
void sort_arr_impl<0,1,2,4,6,3,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6) {
  const size_t ns = d0*d1*d2*d3*d4*d5*d6;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,1>(inp, tmp.get(), d0*d1*d2, d3, d4*d5*d6);
  sort_arr_impl<0,2,1>(tmp.get(), out, d0*d1*d2*d4, d5, d3*d6);
  return;
}

template<>
void sort_arr_impl<0,2,3,4,5,1,6>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6) {
  sort_arr_impl<0,2,1,3>(inp, out, d0, d1, d2*d3*d4*d5, d6);
  return;
}

template<>
void sort_arr_impl<0,2,4,6,5,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6) {
  for (size_t j6 = 0; j6 != d6; ++j6)
    for (size_t j5 = 0; j5 != d5; ++j5)
      for (size_t j4 = 0; j4 != d4; ++j4)
        for (size_t j3 = 0; j3 != d3; ++j3)
          for (size_t j2 = 0; j2 != d2; ++j2)
            for (size_t j1 = 0; j1 != d1; ++j1)
            std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*(j4+d4*(j5+d5*j6))))), d0, out + d0*(j2+d2*(j4+d4*(j6+d6*(j5+d5*(j3+d3*j1))))));
  return;
}

template<>
void sort_arr_impl<1,3,4,5,0,2,6>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6) {
  sort_arr_impl<1,3,0,2,4>(inp, out, d0, d1, d2, d3*d4*d5, d6);
  return;
}

template<>
void sort_arr_impl<2,4,6,0,1,3,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6) {
  const size_t ns = d0*d1*d2*d3*d4*d5*d6;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,4,1,3,5>(inp, tmp.get(), d0*d1, d2, d3, d4, d5, d6);
  sort_arr_impl<1,0>(tmp.get(), out, d0*d1*d3*d5, d2*d4*d6);
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4,5,6,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  std::copy_n(inp, d0*d1*d2*d3*d4*d5*d6*d7, out); 
  return;
}

template<>
void sort_arr_impl<0,1,2,3,6,7,4,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  sort_arr_impl<0,2,1>(inp, out, d0*d1*d2*d3, d4*d5, d6*d7);
  return;
}

template<>
void sort_arr_impl<0,1,3,4,6,7,2,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  sort_arr_impl<0,2,4,1,3>(inp, out, d0*d1, d2, d3*d4, d5, d6*d7);
  return;
}

template<>
void sort_arr_impl<0,1,6,2,3,4,5,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  sort_arr_impl<0,2,1,3>(inp, out, d0*d1, d2*d3*d4*d5, d6, d7);
  return;
}

template<>
void sort_arr_impl<0,2,3,4,5,1,6,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  sort_arr_impl<0,2,1,3>(inp, out, d0, d1, d2*d3*d4*d5, d6*d7);
  return;
}

template<>
void sort_arr_impl<0,4,5,3,6,1,2,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  sort_arr_impl<0,3,2,4,1,5>(inp, out, d0, d1*d2, d3, d4*d5, d6, d7);
}

template<>
void sort_arr_impl<2,3,0,1,4,5,6,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  sort_arr_impl<1,0,2>(inp, out, d0*d1, d2*d3, d4*d5*d6*d7);
}

//template<>
//void sort_arr_impl<2,5,3,1,6,0,7,4>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
//  return;
//}

template<>
void sort_arr_impl<3,4,2,5,1,6,0,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  sort_arr_impl<3,2,4,1,5,0,6>(inp, out, d0, d1, d2, d3*d4, d5, d6, d7);
  return;
}

template<>
void sort_arr_impl<5,6,0,1,2,3,4,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  sort_arr_impl<1,0,2>(inp, out, d0*d1*d2*d3*d4, d5*d6, d7);  
}

template<>
void sort_arr_impl<6,4,2,0,1,3,5,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7) {
  sort_arr_impl<5,3,1,0,2,4,6>(inp, out, d0*d1, d2, d3, d4, d5, d6, d7);  
}

template<>
void sort_arr_impl<0,1,2,3,4,5,6,7,8>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8) {
  std::copy_n(inp, d0*d1*d2*d3*d4*d5*d6*d7*d8, out); 
  return;
}

template<>
void sort_arr_impl<8,6,4,2,0,1,3,5,7>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8) {
  sort_arr_impl<7,5,3,1,0,2,4,6>(inp, out, d0*d1, d2, d3, d4, d5, d6, d7, d8);
  return;
}

template<>
void sort_arr_impl<0,2,4,6,8,7,5,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8) {
  for (size_t j8 = 0; j8 != d8; ++j8)
    for (size_t j7 = 0; j7 != d7; ++j7)
      for (size_t j6 = 0; j6 != d6; ++j6)
        for (size_t j5 = 0; j5 != d5; ++j5)
          for (size_t j4 = 0; j4 != d4; ++j4)
            for (size_t j3 = 0; j3 != d3; ++j3)
              for (size_t j2 = 0; j2 != d2; ++j2)
                for (size_t j1 = 0; j1 != d1; ++j1)
                  std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*(j4+d4*(j5+d5*(j6+d6*(j7+d7*j8))))))), d0, out + d0*(j2+d2*(j4+d4*(j6+d6*(j8+d8*(j7+d7*(j5+d5*(j3+d3*j1))))))));
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4,5,6,7,8,9>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9) {
  std::copy_n(inp, d0*d1*d2*d3*d4*d5*d6*d7*d8*d9, out); 
  return;
}

template<>
void sort_arr_impl<0,1,2,3,4,9,8,7,6,5>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9) {
  sort_arr_impl<0,5,4,3,2,1>(inp, out, d0*d1*d2*d3*d4, d5, d6, d7, d8, d9);
  return;
}

template<>
void sort_arr_impl<0,2,4,6,8,1,3,5,7,9>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9) {
  for (size_t j9 = 0; j9 != d9; ++j9)
    for (size_t j8 = 0; j8 != d8; ++j8)
      for (size_t j7 = 0; j7 != d7; ++j7)
        for (size_t j6 = 0; j6 != d6; ++j6)
          for (size_t j5 = 0; j5 != d5; ++j5)
            for (size_t j4 = 0; j4 != d4; ++j4)
              for (size_t j3 = 0; j3 != d3; ++j3)
                for (size_t j2 = 0; j2 != d2; ++j2)
                  for (size_t j1 = 0; j1 != d1; ++j1)
                    std::copy_n(inp + d0*(j1+d1*(j2+d2*(j3+d3*(j4+d4*(j5+d5*(j6+d6*(j7+d7*(j8+d8*j9)))))))), d0, out + d0*(j2+d2*(j4+d4*(j6+d6*(j8+d8*(j1+d1*(j3+d3*(j5+d5*(j7+d7*j9))))))))); 
 return;
}

template<>
void sort_arr_impl<0,2,4,6,8,9,7,5,3,1>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9) {
  sort_arr_impl<0,2,4,6,8,7,5,3,1>(inp, out, d0, d1, d2, d3, d4, d5, d6, d7, d8*d9);
  return;
}
template<>
void sort_arr_impl<8,6,4,2,0,1,3,5,7,9>(const double* inp, double* out, const size_t d0, const size_t d1, const size_t d2, const size_t d3, const size_t d4, const size_t d5, const size_t d6, const size_t d7, const size_t d8, const size_t d9) {
  const size_t ns = d0*d1*d2*d3*d4*d5*d6*d7*d8*d9;
  std::unique_ptr<double[]> tmp(new double[ns]);
  sort_arr_impl<0,2,4,6,8,7,5,3,1>(inp, tmp.get(), d0*d1, d2, d3, d4, d5, d6, d7, d8, d9);
  sort_arr_impl<1,0>(tmp.get(), out, d0*d1*d3*d5*d7*d9, d2*d4*d6*d8);
  return;
}
