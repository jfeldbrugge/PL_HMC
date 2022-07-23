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
#include "Hamiltonian_Monte_Carlo.h"
#include "Utility.h"

double uniform() {
    return (double) rand() / RAND_MAX;
}
