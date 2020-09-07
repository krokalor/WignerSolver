# Makefile
# Author: KK

CXXFLAGS = -std=c++11 -g -O2 -fopenmp \
	-pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy \
	-Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op \
	-Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast \
	-Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo \
	-Wstrict-null-sentinel -Wswitch-default -Wundef -Werror -Wno-unused #  -Wstrict-overflow=4
CXX = gcc # icpc
OBJS = src/WignerFunction.o \
	src/wfIO.o src/wfBoundCond.o src/wignerTools.o \
	src/solveWignerPoisson.o \
	src/poisson1D.o \
	main.o
LDLIBS = -lsuperlu -larmadillo -lopenblas -llapack -lm -lstdc++ #-m64 -I${MKLROOT}/include # -mkl
#-lsuperlu -L/home/karol/intel/compilers_and_libraries_2019.0.117/linux/mkl/lib -llapack -L/opt/OpenBLAS/lib/ -lm

exec.out: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(OBJS) exec.out
