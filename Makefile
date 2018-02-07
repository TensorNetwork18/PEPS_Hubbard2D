CXX = mpicxx
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -DMKL_ILP64 -m64 -DNDEBUG

MKLDIR=/opt/intel/compilers_and_libraries_2017/mac/mkl
MKLINC=-I$(MKLDIR)/include
MKLLIB=-L$(MKLDIR)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
#MKLLIB=-L$(MKLDIR)/lib -lmkl_intel_lp64 -lmkl_tbb_thread -lstdc++ -lmkl_core -lpthread -lm -ldl

MPILIB=-DHAVE_MPI_H
TBBLIB=-ltbb -DHAVE_TBB #-DTBB_USE_THREADING_TOOLS -DTBB_USE_PERFORMACE_WARNINGS -DTBB_USE_DEBUG

INC_FLAGS=-I. $(MKLINC) $(GTESTINC)
LIB_FLAGS=$(MKLLIB) $(TBBLIB) $(MPILIB)

SRC = mpi/mpi_interface.cc static.cc
SRC_MATH = math/matrix.cc math/tensor.cc math/vec.cc math/sort_arr_impl.cc
SRC_SparseTensor = sparse_tensor/qnum.cc sparse_tensor/basis.cc sparse_tensor/qtensor.cc sparse_tensor/qtensor_tool.cc sparse_tensor/qtensor_permute.cc sparse_tensor/qtensor_contract.cc sparse_tensor/qtensor_svd.cc
SRC_PEPS = peps/peps.cc peps/etensor.cc peps/mps2.cc peps/tmatrix2.cc peps/compress2.cc peps/gauge.cc
SRC_HUBBARD = hubbard/etensor_pool.cc hubbard/block.cc hubbard/effectiveTN.cc hubbard/minimizeE.cc
SRC_GTEST = main.cc

OBJ = $(SRC:.cpp=.o)
OBJ_MATH  = $(SRC_MATH:.cpp=.o)
OBJ_SparseTensor  = $(SRC_SparseTensor:.cpp=.o)
OBJ_PEPS = $(SRC_PEPS:.cpp=.o)
OBJ_HUBBARD = $(SRC_HUBBARD:.cpp=.o)
OBJ_GTEST = $(SRC_GTEST:.cpp=.o)

.cpp.o  :
	$(CXX) $(CXXFLAGS) $(INC_FLAGS) $(LIB_FLAGS) -c $< -o $@

all : test

test : $(OBJ) $(OBJ_MATH) $(OBJ_SparseTensor) $(OBJ_PEPS) $(OBJ_HUBBARD) $(OBJ_GTEST)
	$(CXX) $(CXXFLAGS) -o Hubbard2D.x $(LIB_FLAGS) $(OBJ) $(OBJ_MATH) $(OBJ_SparseTensor) $(OBJ_PEPS) $(OBJ_HUBBARD) $(OBJ_GTEST)

clean: 
	rm -f Hubbard2D.x *.o
