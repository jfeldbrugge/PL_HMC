template <size_t dimensions>
void flow(const point<dimensions> x,
          point<dimensions> &z,
          matrix<dimensions> &J,
          point<dimensions> &f,
          const double tau, const int N_tau)
{
    
    std::complex<double> eps = tau / N_tau;
 
    std::vector<point<dimensions>> zi(N_tau + 1);
    std::vector<matrix<dimensions>> Ji(N_tau + 1);
    std::vector<point<dimensions>> fi(N_tau + 1);
    
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
    
    fi[0] = gradS(x) * I;
    for (int i = 1; i < N_tau + 1; i++) {
        fi[i] = fi[i - 1] - multiply(fi[i - 1].conjugate(), hessS(zi[N_tau - i])) * I * eps;
    }
    f = fi[N_tau];
    
//    for (int i = 0; i < zi.size(); i++) {
//        std::cout << "zi[" << i << "] = " << zi[i] << std::endl;
//    }
//    std::cout << std::endl;
//    for (int i = 0; i < Ji.size(); i++) {
//        std::cout << "J[" << i << "] = " << Ji[i] << std::endl;
//    }
//    std::cout << std::endl;
//    for (int i = 0; i < Ji.size(); i++) {
//        std::cout << "f[" << i << "] = " << fi[i] << std::endl;
//    }
}

template <size_t dimensions>
point<dimensions> force(const point<dimensions> x, const double tau, const int N_tau) {
    point<dimensions> z, f;
    matrix<dimensions> J;
    
    flow(x, z, J, f, tau, N_tau);
    return multiply(f, J).re();
}

template <size_t dimensions>
double h(const point<dimensions> p) {
    return std::real(I * S(p));
}

template <size_t dimensions>
double H(const point<dimensions> p) {
    return std::imag(I * S(p));
}
