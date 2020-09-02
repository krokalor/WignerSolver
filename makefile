# Makefile
# Author: KK

CXXFLAGS = -std=c++11 -g -O2 -fopenmp -Wall -Wpedantic -Werror #-Wextra
CXX = gcc # icpc
OBJS = src/WignerFunction.o \
	src/wfIO.o src/wfBoundCond.o src/wignerTools.o \
	src/solveWignerPoisson.o \
	src/poisson1D.o \
	main.o
LDLIBS = -lsuperlu -larmadillo -lopenblas -llapack -lm -lstdc++ # -mkl
#-lsuperlu -L/home/karol/intel/compilers_and_libraries_2019.0.117/linux/mkl/lib -llapack -L/opt/OpenBLAS/lib/ -lm

exec.out: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(OBJS) exec.out
