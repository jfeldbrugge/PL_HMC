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
    
    point<dimensions> F = force(X, tau, N_tau);
    for (int i = 1; i < N_s + 1; i++) {
        P = P + F * Delta_s / 2.;
        X = X + P * Delta_s / m;
        F = force(X, tau, N_tau);
        P = P + F * Delta_s / 2.;
    }

    return {X, Hamiltonian(X, P, m, tau, N_tau)};
}

template <size_t dimensions>
std::complex<double> expectation(std::complex<double> (*OO)(point<dimensions>), const std::vector<point<dimensions>> &xi, const double tau, const int N_tau) {
    std::complex<double> numerator = 0., denominator = 0.;
    for (int i = 0; i < xi.size(); i++) {
        point<dimensions> z, f;
        matrix<dimensions> J;
        flow(xi[i], z, J, f, tau, N_tau);
        numerator = numerator + OO(z) * std::exp(I * H(z)) * J.determinant();
        denominator = denominator + std::exp(I * H(z)) * J.determinant();
    }
    
    return numerator / denominator;
}
