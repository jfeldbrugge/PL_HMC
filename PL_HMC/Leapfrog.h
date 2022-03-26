template <size_t dimensions>
double Hamiltonian(point<dimensions> z, point<dimensions> p, const double m) {
    
    double p2 = 0;
    for(size_t i = 0; i < dimensions; ++i)
    {
        p2 = p2 + std::real(p.get(i) * std::conj(p.get(i)));
    }
    return p2 / (2. * m) - h(z);
}

template <size_t dimensions>
double Hamiltonian(point<dimensions> x, point<dimensions> p, const double m, const double tau, const int N_tau) {
    point<dimensions> z, f;
    matrix<dimensions> J;
    flow(x, z, J, f, tau, N_tau);
    return Hamiltonian(z, p, m);
}
    
template <size_t dimensions>
std::tuple<point<dimensions>, double> leapfrog(const point<dimensions> x, const point<dimensions> p,
              const double Delta_s, const int N_s, const double tau, const int N_tau, const double m) {
    point<dimensions> X = x, P = p;
    for (int i = 1; i < N_s + 1; i++) {
        P = P + force(X, tau, N_tau) * Delta_s / 2.;
        X = X + P * Delta_s / m;
        P = P + force(X, tau, N_tau) * Delta_s / 2.;
        
//        std::cout << "(X, P, H) = (" << X << ", " << P << ", " << Hamiltonian(X, P, m, tau, N_tau) << ")" << std::endl;
    }
    
    return {X, Hamiltonian(X, P, m, tau, N_tau)};
}
