template <size_t dimensions>
double h(const point<dimensions> p) {
    return std::real(I * S(p));
}

template <size_t dimensions>
double H(const point<dimensions> p) {
    return std::imag(I * S(p));
}

template <size_t dimensions>
point<dimensions> gradS(const point<dimensions> p) {
    const double eps = 0.000001;
    std::array<std::complex<double>, dimensions> grad;
    
    for (int i = 0; i < dimensions; i++) {
        std::array<std::complex<double>, dimensions> temp;
        temp[i] = eps;
        point<dimensions> Peps(temp);
        grad[i] = (S(p + Peps) - S(p - Peps)) / (2. * eps);
    }
    
    return point<dimensions>(grad);
}

template <size_t dimensions>
matrix<dimensions> hessS(const point<dimensions> p) {
    const double eps = 0.0001;
    std::array<std::complex<double>, dimensions * dimensions> hessian;
    for (int i = 0; i < dimensions; i++) {
        for (int j = 0; j < dimensions; j++) {
            
            if (i == j) {
                std::array<std::complex<double>, dimensions> temp;
                temp[i] = eps;
                point<dimensions> Peps(temp);
                hessian[i * dimensions + j] = (S(p + Peps) -2. * S(p) + S(p - Peps)) / (eps * eps);
                
            } else {
                std::array<std::complex<double>, dimensions> temp1, temp2;
                temp1[i] = eps; temp2[j] = eps;
                point<dimensions> Peps1(temp1), Peps2(temp2);
                hessian[i * dimensions + j] = (S(p + Peps1 + Peps2) - S(p + Peps1 - Peps2) - S(p - Peps1 + Peps2) + S(p - Peps1 - Peps2)) / (4. * eps * eps);
            }
        }
    }
    return matrix<dimensions>(hessian);
}

template <size_t dimensions>
void flow(const point<dimensions> x, point<dimensions> &z, const double tau, const int N_tau) {
    std::complex<double> eps = tau / N_tau;
 
    std::vector<point<dimensions>> zi(N_tau + 1);
    std::vector<matrix<dimensions>> Ji(N_tau + 1);
    std::vector<point<dimensions>> fi(N_tau + 1);
    
    zi[0] = x;
    for (int i = 1; i < N_tau + 1; i++) {
        zi[i] = zi[i - 1] + gradS(zi[i - 1]).conjugate() * I * eps;
    }
    z = zi[N_tau];
    
//    for (int i = 0; i < zi.size(); i++) {
//        std::cout << "zi[" << i << "] = " << zi[i] << std::endl;
//    }
}

template <size_t dimensions>
void flow(const point<dimensions> x, point<dimensions> &z, matrix<dimensions> &J, const double tau, const int N_tau) {
    std::complex<double> eps = tau / N_tau;
 
    std::vector<point<dimensions>> zi(N_tau + 1);
    std::vector<matrix<dimensions>> Ji(N_tau + 1);
    
    zi[0] = x;
    for (int i = 1; i < N_tau + 1; i++) {
        zi[i] = zi[i - 1] + gradS(zi[i - 1]).conjugate() * I * eps;
    }
    z = zi[N_tau];
    
    Ji[0] = identity<dimensions>();
    for (int i = 1; i < N_tau + 1; i++) {
        Ji[i] = Ji[i - 1] + hessS(zi[i - 1]).conjugate() * Ji[i - 1].conjugate() * I * eps;
    }
    J = Ji[N_tau];
    
//    for (int i = 0; i < zi.size(); i++) {
//        std::cout << "zi[" << i << "] = " << zi[i] << std::endl;
//    }
//    std::cout << std::endl;
//    for (int i = 0; i < Ji.size(); i++) {
//        std::cout << "J[" << i << "] = " << Ji[i] << std::endl;
//    }
//    std::cout << std::endl;
}

template <size_t dimensions>
point<dimensions> force(const point<dimensions> x, const double tau, const int N_tau) {
    std::complex<double> eps = tau / N_tau;
 
    std::vector<point<dimensions>> zi(N_tau + 1);
    std::vector<point<dimensions>> fi(N_tau + 1);
    
    zi[0] = x;
    for (int i = 1; i < N_tau + 1; i++) {
        zi[i] = zi[i - 1] + gradS(zi[i - 1]).conjugate() * I * eps;
    }
    
    fi[N_tau] = (gradS(zi.back()) * I).re();
    for (int i = N_tau - 1; i >= 0; i--) {
        fi[i] = fi[i + 1] - multiply(fi[i + 1].conjugate(), hessS(zi[i])) * I * eps;
    }

    return fi[0].re() * 2;
}

template <size_t dimensions>
double hamiltonian(point<dimensions> z, point<dimensions> p, const matrix<dimensions> M_inv) {
    return std::real(M_inv * p * p / 2.) - h(z);
}

template <size_t dimensions>
double hamiltonian(point<dimensions> x, point<dimensions> p, const matrix<dimensions> M_inv, const double tau, const int N_tau) {
    point<dimensions> z, f;
    flow(x, z, tau, N_tau);
    return hamiltonian(z, p, M_inv);
}
    
template <size_t dimensions>
std::tuple<point<dimensions>, double> leapfrog(const point<dimensions> x, const point<dimensions> p,
              const double Delta_s, const int N_s, const double tau, const int N_tau, const matrix<dimensions> Minv) {
    point<dimensions> X = x, P = p;
    point<dimensions> F = force(X, tau, N_tau);
    for (int i = 1; i < N_s + 1; i++) {
        // Kick
        P = P + F * Delta_s / 2.;
        // Drift
        X = X + Minv * P * Delta_s;
        //Kick
        F = force(X, tau, N_tau);
        P = P + F * Delta_s / 2.;
    }
    return {X, hamiltonian(X, P, Minv, tau, N_tau)};
}

template <size_t dimensions>
std::complex<double> expectation(std::complex<double> (*OO)(point<dimensions>),
                                 const std::vector<point<dimensions>> &xi, const double tau, const int N_tau) {
    std::complex<double> numerator = 0., denominator = 0.;
    for (int i = 0; i < xi.size(); i++) {
        point<dimensions> z;
        matrix<dimensions> J;
        flow(xi[i], z, J, tau, N_tau);
        numerator   = numerator   + OO(z) * std::exp(I * H(z)) * J.determinant();
        denominator = denominator +         std::exp(I * H(z)) * J.determinant();
    }
    return numerator / denominator;
}

template <size_t dimensions>
std::complex<double> variance(std::complex<double> (*OO)(point<dimensions>), const std::complex<double> exp,
                              const std::vector<point<dimensions>> &xi, const double tau, const int N_tau) {
    std::complex<double> numerator = 0., denominator = 0.;
    for (int i = 0; i < xi.size(); i++) {
        point<dimensions> z, f;
        matrix<dimensions> J;
        flow(xi[i], z, J, tau, N_tau);
        numerator   = numerator   + (OO(z) - exp) * std::conj(OO(z) - exp) * std::exp(I * H(z)) * J.determinant();
        denominator = denominator +                                          std::exp(I * H(z)) * J.determinant();
    }
    return numerator / denominator;
}
