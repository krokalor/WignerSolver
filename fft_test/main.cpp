#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#define ARMA_USE_SUPERLU 1
#include <armadillo>

using std::cout;
using std::endl;
using namespace arma;


int main() {

    double a = 1, x_max = 4, x = 0;
    size_t N = 100;
    double dx = 2*x_max/(N-1);

    vec X(N, fill::zeros);

    for(size_t i = 0; i < N; ++i) {
        x = -x_max + i*dx;
        X(i) = exp(-a*x*x);
    }

    vec Y = real(fft(X, N));

    X.print();
    cout<<endl<<endl;
    Y.print();

}
