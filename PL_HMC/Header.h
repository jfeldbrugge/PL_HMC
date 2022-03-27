#include <iostream>
#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <array>
#include <algorithm>
#include <random>
#include <time.h>
#include <stdlib.h>
#include <tuple>

const double epsilon = std::numeric_limits<double>::epsilon() ;
const double pi = 3.141592654;
const std::complex<double> I (0,1.0);

#include "Point.h"
#include "Matrix.h"
#include "Flow.h"
#include "Gradient.h"
#include "Leapfrog.h"
#include "Utility.h"

template <size_t dimensions>
std::complex<double> S(const point<dimensions> p);

template <size_t dimensions>
double h(const point<dimensions> p);

template <size_t dimensions>
double H(const point<dimensions> p);

template <size_t dimensions>
double hamiltonian(point<dimensions> z, point<dimensions> p, const double m);

template <size_t dimensions>
point<dimensions> gradS(const point<dimensions> p);

template <size_t dimensions>
matrix<dimensions> hessS(const point<dimensions> p);

template <size_t dimensions>
void flow(const point<dimensions> x,
          point<dimensions> &z,
          matrix<dimensions> &J,
          point<dimensions> &f,
          const double tau, const int N_tau);

template <size_t dimensions>
point<dimensions> force(const point<dimensions> x, const double tau, const int N_tau);

template <size_t dimensions>
std::tuple<point<dimensions>, double> leapfrog(const point<dimensions> x, const point<dimensions> p,
                                               const double Delta_s, const int N_s, const double tau, const int N_tau, const double m);

double uniform() {
    return (double) rand() / RAND_MAX;
}

template<size_t dimensions>
void writeB(std::vector<point<dimensions>> xi, std::string fileName);

template <size_t dimensions>
std::complex<double> expectation(std::complex<double> (*OO)(point<dimensions>), const std::vector<point<dimensions>> &xi, const double tau, const int N_tau);
